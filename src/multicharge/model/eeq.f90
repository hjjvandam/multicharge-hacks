! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file multicharge/model/eeq.f90
!> Provides implementation of the electronegativity equilibration model (EEQ)

!> Electronegativity equlibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
module multicharge_model_eeq
   use mctc_env, only: error_type, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3
   use mctc_ncoord, only: new_ncoord, cn_count
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   use multicharge_ewald, only: get_alpha
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans, get_rec_trans
   use multicharge_model_cache, only: cache_container, model_cache
   implicit none
   private

   public :: eeq_model, new_eeq_model

   type, extends(model_cache), public :: eeq_cache
   end type eeq_cache

   type, extends(mchrg_model_type) :: eeq_model
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate EN vector derivatives
      procedure :: get_xvec_derivs
   end type eeq_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

subroutine new_eeq_model(self, mol, error, chi, rad, eta, kcnchi, &
   & cutoff, cn_exp, rcov, cn_max)
   !> Electronegativity equilibration model
   type(eeq_model), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Electronegativity
   real(wp), intent(in) :: chi(:)
   !> Exponent gaussian charge
   real(wp), intent(in) :: rad(:)
   !> Chemical hardness
   real(wp), intent(in) :: eta(:)
   !> CN scaling factor for electronegativity
   real(wp), intent(in) :: kcnchi(:)
   !> Cutoff radius for coordination number
   real(wp), intent(in), optional :: cutoff
   !> Steepness of the CN counting function
   real(wp), intent(in), optional :: cn_exp
   !> Covalent radii for CN
   real(wp), intent(in), optional :: rcov(:)
   !> Maximum CN cutoff for CN
   real(wp), intent(in), optional :: cn_max

   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi

   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max)

end subroutine new_eeq_model

subroutine update(self, mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in), optional :: qloc(:)
   real(wp), intent(in), optional :: dcndr(:, :, :)
   real(wp), intent(in), optional :: dcndL(:, :, :)
   real(wp), intent(in), optional :: dqlocdr(:, :, :)
   real(wp), intent(in), optional :: dqlocdL(:, :, :)

   type(eeq_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Refer CN arrays in cache
   ptr%cn = cn
   if (present(dcndr) .and. present(dcndL)) then
      ptr%dcndr = dcndr
      ptr%dcndL = dcndL
   end if

   if (any(mol%periodic)) then
      ! Create WSC
      call new_wignerseitz_cell(ptr%wsc, mol)
      call get_alpha(mol%lattice, ptr%alpha)
   end if

end subroutine update

subroutine get_xvec(self, mol, cache, xvec)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: xvec(:)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, xvec, ptr) private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) / sqrt(ptr%cn(iat) + reg)
      xvec(iat) = -self%chi(izp) + tmp * ptr%cn(iat)
   end do
   xvec(mol%nat + 1) = mol%charge

end subroutine get_xvec

subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out), contiguous :: dxdr(:, :, :)
   real(wp), intent(out), contiguous :: dxdL(:, :, :)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   dxdr(:, :, :) = 0.0_wp
   dxdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, ptr, dxdr, dxdL) &
   !$omp private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) / sqrt(ptr%cn(iat) + reg)
      dxdr(:, :, iat) = 0.5_wp * tmp * ptr%dcndr(:, :, iat) + dxdr(:, :, iat)
      dxdL(:, :, iat) = 0.5_wp * tmp * ptr%dcndL(:, :, iat) + dxdL(:, :, iat)
   end do

end subroutine get_xvec_derivs

subroutine get_coulomb_matrix(self, mol, cache, amat)
   use omp_lib
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: amat(:, :)

   type(eeq_cache), pointer :: ptr

   real(wp) :: t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12

   !DEBUG
   real(wp), allocatable :: amat01(:,:)
   real(wp), allocatable :: amat02(:,:)
   real(wp), allocatable :: amat03(:,:)
   real(wp), allocatable :: amat04(:,:)
   real(wp), allocatable :: amat05(:,:)
   real(wp), allocatable :: amat06(:,:)
   real(wp), allocatable :: amat07(:,:)
   real(wp), allocatable :: amat08(:,:)
   real(wp), allocatable :: amat09(:,:)
   real(wp), allocatable :: amat10(:,:)
   real(wp), allocatable :: amat11(:,:)

   allocate(amat01, source=amat)
   allocate(amat02, source=amat)
   allocate(amat03, source=amat)
   allocate(amat04, source=amat)
   allocate(amat05, source=amat)
   allocate(amat06, source=amat)
   allocate(amat07, source=amat)
   allocate(amat08, source=amat)
   allocate(amat09, source=amat)
   allocate(amat10, source=amat)
   allocate(amat11, source=amat)
   !DEBUG

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, ptr%wsc, ptr%alpha, amat)
   else
      t0 = omp_get_wtime()
      call get_amat_0d(self, mol, amat)
      t1 = omp_get_wtime()
      call get_amat_0ds(self, mol, amat01)
      t2 = omp_get_wtime()
      call get_amat_0d_1(self, mol, amat02)
      t3 = omp_get_wtime()
      call get_amat_0d_1a(self, mol, amat03)
      t4 = omp_get_wtime()
      call get_amat_0d_1c(self, mol, amat04)
      t5 = omp_get_wtime()
      call get_amat_0d_1d(self, mol, amat05)
      t6 = omp_get_wtime()
      call get_amat_0d_2a(self, mol, amat06)
      t7 = omp_get_wtime()
      call get_amat_0d_3a(self, mol, amat07)
      t8 = omp_get_wtime()
      call get_amat_0d_3b(self, mol, amat08)
      t9 = omp_get_wtime()
      call get_amat_0d_4a(self, mol, amat09)
      t10 = omp_get_wtime()
      call get_amat_0d_4c(self, mol, amat10)
      t11 = omp_get_wtime()
      call get_amat_0d_4d(self, mol, amat11)
      t12 = omp_get_wtime()
      write(*,*)"timings (s) resolution: ",omp_get_wtick()
      write(*,'("get_amat_0d   : ",f12.6,"  ",d12.6)')t1-t0,0.0d0
      write(*,'("get_amat_0ds  : ",f12.6,"  ",d12.6)')t2-t1,norm2(amat01-amat)
      write(*,'("get_amat_0d_1 : ",f12.6,"  ",d12.6)')t3-t2,norm2(amat02-amat)
      write(*,'("get_amat_0d_1a: ",f12.6,"  ",d12.6)')t4-t3,norm2(amat03-amat)
      write(*,'("get_amat_0d_1c: ",f12.6,"  ",d12.6)')t5-t4,norm2(amat04-amat)
      write(*,'("get_amat_0d_1d: ",f12.6,"  ",d12.6)')t6-t5,norm2(amat05-amat)
      write(*,'("get_amat_0d_2a: ",f12.6,"  ",d12.6)')t7-t6,norm2(amat06-amat)
      write(*,'("get_amat_0d_3a: ",f12.6,"  ",d12.6)')t8-t7,norm2(amat07-amat)
      write(*,'("get_amat_0d_3b: ",f12.6,"  ",d12.6)')t9-t8,norm2(amat08-amat)
      write(*,'("get_amat_0d_4a: ",f12.6,"  ",d12.6)')t10-t9,norm2(amat09-amat)
      write(*,'("get_amat_0d_4c: ",f12.6,"  ",d12.6)')t11-t10,norm2(amat10-amat)
      write(*,'("get_amat_0d_4d: ",f12.6,"  ",d12.6)')t12-t11,norm2(amat11-amat)
      !write(*,*)"get_amat_0d_4b: ",t10-t9
   end if
   !DEBUG
   deallocate(amat01)
   deallocate(amat02)
   deallocate(amat03)
   deallocate(amat04)
   deallocate(amat05)
   deallocate(amat06)
   deallocate(amat07)
   deallocate(amat08)
   deallocate(amat09)
   deallocate(amat10)
   !DEBUG
end subroutine get_coulomb_matrix

subroutine get_amat_0d(self, mol, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   ! Thread-private array for reduction
   real(wp), allocatable :: amat_local(:, :)

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, amat_local)
   allocate(amat_local, source=amat)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat_local(jat, iat) = amat_local(jat, iat) + tmp
         amat_local(iat, jat) = amat_local(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat_local(iat, iat) = amat_local(iat, iat) + tmp
   end do
   !$omp end do
   !$omp critical (get_amat_0d_)
   amat(:, :) = amat(:, :) + amat_local(:, :)
   !$omp end critical (get_amat_0d_)
   deallocate(amat_local)
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d

subroutine get_amat_0ds(self, mol, amat)
   !
   ! Same as get_amat_0d except that we use static scheduling
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   ! Thread-private array for reduction
   real(wp), allocatable :: amat_local(:, :)

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, amat_local)
   allocate(amat_local, source=amat)
   !$omp do schedule(static)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat_local(jat, iat) = amat_local(jat, iat) + tmp
         amat_local(iat, jat) = amat_local(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat_local(iat, iat) = amat_local(iat, iat) + tmp
   end do
   !$omp end do
   !$omp critical (get_amat_0d_)
   amat(:, :) = amat(:, :) + amat_local(:, :)
   !$omp end critical (get_amat_0d_)
   deallocate(amat_local)
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0ds

subroutine get_amat_0d_1(self, mol, amat)
   !
   ! Same as get_amat_0d except that we use 1 shared 3D buffer with a
   ! separate matrix for each thread. Instead of using critical regions
   ! to accumulate the results we use 3 loops, the outer 2 of which we
   ! parallelise with OpenMP.
   !
   use omp_lib
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   ! Thread-private array for reduction
   real(wp), allocatable :: amat_threads(:, :, :)
   integer :: nsize(3)
   integer :: idth

   amat(:, :) = 0.0_wp

   nsize(1:2) = shape(amat)
   nsize(3) = omp_get_max_threads()
   allocate(amat_threads(nsize(1),nsize(2),nsize(3)))

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self, amat_threads) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, idth)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      idth = omp_get_thread_num()+1
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat_threads(jat, iat, idth) = amat_threads(jat, iat, idth) + tmp
         amat_threads(iat, jat, idth) = amat_threads(iat, jat, idth) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat_threads(iat, iat, idth) = amat_threads(iat, iat, idth) + tmp
   end do
   !$omp end do
   !$omp end parallel
   do idth = 1, omp_get_max_threads()
     !$omp parallel default(none) &
     !$omp shared(amat, amat_threads, idth, mol) &
     !$omp private(iat, jat)
     !$omp do collapse(2) schedule(static)
     do iat = 1, mol%nat
       do jat = 1, mol%nat
         amat(jat, iat) = amat(jat, iat) + amat_threads(jat, iat, idth)
       enddo
     enddo
     !$omp end do
     !$omp end parallel
   enddo
   deallocate(amat_threads)

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_1

subroutine get_amat_0d_1a(self, mol, amat)
   !
   ! Same as get_amat_0d except that we eliminate the buffers for
   ! partial results altogether.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat(jat, iat) = amat(jat, iat) + tmp
         amat(iat, jat) = amat(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat(iat, iat) = amat(iat, iat) + tmp
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_1a

! OpenMP cannot collapse two loops when the loop limit of the second
! depends on the counter value of the first.
!
!subroutine get_amat_0d_1b(self, mol, amat)
!   ! Here we eliminate the intermediate matrix amat_local altogether.
!   ! This approach might not work if calculating the contributions
!   ! to the matrix involves a lot of parallelizable work, but it
!   ! this case it is easy.
!   !
!   ! Building on get_amat_0d_1a we now collapse and parallelize
!   ! over the iat and jat loop.
!   class(eeq_model), intent(in) :: self
!   type(structure_type), intent(in) :: mol
!   real(wp), intent(out) :: amat(:, :)
!
!   integer :: iat, jat, izp, jzp
!   real(wp) :: vec(3), r2, gam, tmp
!
!   amat(:, :) = 0.0_wp
!
!   !$omp parallel default(none) &
!   !$omp shared(amat, mol, self) &
!   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
!   !$omp do collapse(2) schedule(runtime)
!   do iat = 1, mol%nat
!      do jat = 1, iat
!         izp = mol%id(iat)
!         if (iat.ne.jat) then
!           jzp = mol%id(jat)
!           vec = mol%xyz(:, jat) - mol%xyz(:, iat)
!           r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
!           gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
!           tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
!           amat(jat, iat) = amat(jat, iat) + tmp
!           amat(iat, jat) = amat(iat, jat) + tmp
!         else
!           tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
!           amat(iat, iat) = amat(iat, iat) + tmp
!         endif
!      end do
!   end do
!   !$omp end do
!   !$omp end parallel
!
!   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
!   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
!   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
!
!end subroutine get_amat_0d_1b

subroutine get_amat_0d_1c(self, mol, amat)
   !
   ! Same as get_amat_0d_1 except that we replace the 2 loops over the
   ! upper triangle with 1 loop as demonstrated in
   ! `multicharge/snippets/prog_single.f90`.
   ! In addition we use static scheduling.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: izp, jzp
   integer(kind=8) :: iat, jat, idx, ntop
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   ntop = mol%nat*(mol%nat+1)/2
   !$omp parallel default(none) &
   !$omp shared(amat, mol, self, ntop) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, idx)
   !$omp do schedule(static)
   do idx = 1, ntop
         iat = int((1.0d0 + sqrt(1.0d0 + 8.0d0 * real(idx-1, 8))) / 2.0d0)
         jat = idx - iat*(iat-1)/2
         izp = mol%id(iat)
         if (iat.ne.jat) then
           jzp = mol%id(jat)
           vec = mol%xyz(:, jat) - mol%xyz(:, iat)
           r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
           gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
           tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
           amat(jat, iat) = amat(jat, iat) + tmp
           amat(iat, jat) = amat(iat, jat) + tmp
         else
           tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
           amat(iat, iat) = amat(iat, iat) + tmp
         endif
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_1c

subroutine get_amat_0d_1d(self, mol, amat)
   !
   ! Same as get_amat_0d_1 except that we replace the 2 loops over the
   ! upper triangle with 2 loops running over a rectangle as demonstrated
   ! in `multicharge/snippets/prog_mapping.f90`.
   ! In addition we collapse the 2 loops and use static scheduling.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: jcol, irow
   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do collapse(2) schedule(static)
   do jcol = 1, mol%nat
      do irow = 1, (mol%nat+1)/2
         if (irow .ge. jcol) then
           jat = irow
           iat = jcol
           izp = mol%id(iat)
           if (iat.ne.jat) then
             jzp = mol%id(jat)
             vec = mol%xyz(:, jat) - mol%xyz(:, iat)
             r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
             gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
             tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
             amat(jat, iat) = amat(jat, iat) + tmp
             amat(iat, jat) = amat(iat, jat) + tmp
           else
             tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
             amat(iat, iat) = amat(iat, iat) + tmp
           endif
         endif
         if (irow .le. mol%nat/2 .and. irow .le. jcol) then
           jat = mol%nat - irow + 1
           iat = mol%nat - jcol + 1
           izp = mol%id(iat)
           if (iat.ne.jat) then
             jzp = mol%id(jat)
             vec = mol%xyz(:, jat) - mol%xyz(:, iat)
             r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
             gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
             tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
             amat(jat, iat) = amat(jat, iat) + tmp
             amat(iat, jat) = amat(iat, jat) + tmp
           else
             tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
             amat(iat, iat) = amat(iat, iat) + tmp
           endif
         endif
      end do
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_1d

subroutine get_amat_0d_2a(self, mol, amat)
   !
   ! Same as get_amat_0d_1a except that we run the outer loop over the
   ! upper triangle backwards, thus ensuring the largest tasks are
   ! scheduled first to improve load balancing.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do schedule(runtime)
   do iat = mol%nat, 1, -1
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat(jat, iat) = amat(jat, iat) + tmp
         amat(iat, jat) = amat(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat(iat, iat) = amat(iat, iat) + tmp
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_2a

! This doesn't work, see the comments at get_amat_0d_1b for the
! reasons "why".
!
!subroutine get_amat_0d_2b(self, mol, amat)
!   ! Building on get_amat_0d_1a we improve load balancing
!   ! by running the loop over which we parallelise backwards.
!   ! This means that the biggest tasks are started first,
!   ! and the smallest last. In combination with dynamic
!   ! load balancing this should improve the time to solution
!   ! as the load imbalance is proportional to the size of the
!   ! last tasks.
!   class(eeq_model), intent(in) :: self
!   type(structure_type), intent(in) :: mol
!   real(wp), intent(out) :: amat(:, :)
!
!   integer :: iat, jat, izp, jzp
!   real(wp) :: vec(3), r2, gam, tmp
!
!   amat(:, :) = 0.0_wp
!
!   !$omp parallel default(none) &
!   !$omp shared(amat, mol, self) &
!   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
!   !$omp do collapse(2) schedule(runtime)
!   do iat = mol%nat, 1, -1
!      do jat = 1, iat
!         izp = mol%id(iat)
!         if (iat.ne.jat) then
!           jzp = mol%id(jat)
!           vec = mol%xyz(:, jat) - mol%xyz(:, iat)
!           r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
!           gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
!           tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
!           amat(jat, iat) = amat(jat, iat) + tmp
!           amat(iat, jat) = amat(iat, jat) + tmp
!         else
!           tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
!           amat(iat, iat) = amat(iat, iat) + tmp
!         endif
!      end do
!   end do
!   !$omp end do
!   !$omp end parallel
!
!   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
!   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
!   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
!
!end subroutine get_amat_0d_2b

subroutine get_amat_0d_3a(self, mol, amat)
   !
   ! Same as get_amat_0d_1a except that we calculate the entire
   ! matrix explicitly. This avoids the transposed matrix accesses
   ! at the cost of duplicating the compute.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, mol%nat
         if (iat.ne.jat) then
           jzp = mol%id(jat)
           vec = mol%xyz(:, jat) - mol%xyz(:, iat)
           r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
           gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
           tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
           amat(jat, iat) = amat(jat, iat) + tmp
         else
           tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
           amat(iat, iat) = amat(iat, iat) + tmp
         endif
      end do
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_3a

subroutine get_amat_0d_3b(self, mol, amat)
   !
   ! Same as get_amat_0d_3a except that we collapse the 2 loops
   ! and use static scheduling.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do collapse(2) schedule(static)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         izp = mol%id(iat)
         if (iat.ne.jat) then
           jzp = mol%id(jat)
           vec = mol%xyz(:, jat) - mol%xyz(:, iat)
           r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
           gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
           tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
           amat(jat, iat) = amat(jat, iat) + tmp
         else
           tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
           amat(iat, iat) = amat(iat, iat) + tmp
         endif
      end do
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_3b

subroutine get_amat_0d_4a(self, mol, amat)
   !
   ! Same as get_amat_0d_2a except that we first calculate
   ! just the upper triangle. Afterwards we populate the
   ! lower triangle by copying the data.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do schedule(runtime)
   do iat = mol%nat, 1, -1
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
         amat(jat, iat) = amat(jat, iat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat(iat, iat) = amat(iat, iat) + tmp
   end do
   !$omp end do
   !$omp end parallel

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, jat)
   !$omp do schedule(runtime)
   do iat = mol%nat, 1, -1
      do jat = 1, iat - 1
         amat(iat, jat) = amat(jat, iat)
      end do
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_4a

! as feared OpenMP currently cannot collapse loops when the
! loop limit of one loop depends on the iteration number of
! the other (not even for the special case of a triangle).
!subroutine get_amat_0d_4b(self, mol, amat)
!   ! Building on get_amat_0d_1a we improve load balancing
!   ! by running the loop over which we parallelise backwards.
!   ! This means that the biggest tasks are started first,
!   ! and the smallest last. In combination with dynamic
!   ! load balancing this should improve the time to solution
!   ! as the load imbalance is proportional to the size of the
!   ! last tasks.
!   !
!   ! In addition we just calculate the triangle first and then
!   ! copy the results in a dedicated loop.
!   class(eeq_model), intent(in) :: self
!   type(structure_type), intent(in) :: mol
!   real(wp), intent(out) :: amat(:, :)
!
!   integer :: iat, jat, izp, jzp
!   real(wp) :: vec(3), r2, gam, tmp
!
!   amat(:, :) = 0.0_wp
!
!   !$omp parallel default(none) &
!   !$omp shared(amat, mol, self) &
!   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
!   !$omp do collapse(2) schedule(runtime)
!   do iat = mol%nat, 1, -1
!      do jat = 1, iat - 1
!         izp = mol%id(iat)
!         jzp = mol%id(jat)
!         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
!         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
!         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
!         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
!         amat(jat, iat) = amat(jat, iat) + tmp
!      end do
!      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
!      amat(iat, iat) = amat(iat, iat) + tmp
!   end do
!   !$omp end do
!   !$omp end parallel
!
!   !$omp parallel default(none) &
!   !$omp shared(amat, mol, self) &
!   !$omp private(iat, jat)
!   !$omp do collapse(2) schedule(runtime)
!   do iat = mol%nat, 1, -1
!      do jat = 1, iat - 1
!         amat(iat, jat) = amat(jat, iat)
!      end do
!   end do
!   !$omp end do
!   !$omp end parallel
!
!   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
!   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
!   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
!
!end subroutine get_amat_0d_4b

subroutine get_amat_0d_4c(self, mol, amat)
   !
   ! Same as get_amat_0d_4a except that we manually collapse
   ! the 2 loops over the upper triangle as shown in
   ! multicharge/snippets/prog_single.f90.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: izp, jzp
   integer(kind=8) :: iat, jat, ntop, idx
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   ntop = mol%nat*(mol%nat+1)/2
   !$omp parallel default(none) &
   !$omp shared(amat, mol, self, ntop) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
   !$omp do schedule(static)
   do idx = 1, ntop
         iat = int((1.0d0 + sqrt(1.0d0 + 8.0d0 * real(idx-1, 8))) / 2.0d0)
         jat = idx - iat*(iat-1)/2
         izp = mol%id(iat)
         jzp = mol%id(jat)
         if (iat.ne.jat) then
           vec = mol%xyz(:, jat) - mol%xyz(:, iat)
           r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
           gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
           tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
           amat(jat, iat) = amat(jat, iat) + tmp
         else
           tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
           amat(iat, iat) = amat(iat, iat) + tmp
         endif
   end do
   !$omp end do
   !$omp end parallel

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self, ntop) &
   !$omp private(iat, jat, idx)
   !$omp do schedule(static)
   do idx = 1, ntop
         iat = int((1.0d0 + sqrt(1.0d0 + 8.0d0 * real(idx-1, 8))) / 2.0d0)
         jat = idx - iat*(iat-1)/2
         amat(iat, jat) = amat(jat, iat)
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_4c

subroutine get_amat_0d_4d(self, mol, amat)
   !
   ! Same as get_amat_0d_4a except that we loop over
   ! the rectange that represents the upper half of
   ! the matrix, as shown in
   ! multicharge/snippets/prog_mapping.f90.
   !
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, irow, jcol
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, irow, jcol)
   !$omp do collapse(2) schedule(static)
   do jcol = 1, mol%nat
      do irow = 1, (mol%nat+1)/2
         if (irow .ge. jcol) then
           jat = irow
           iat = jcol
           izp = mol%id(iat)
           jzp = mol%id(jat)
           if (iat.ne.jat) then
             vec = mol%xyz(:, jat) - mol%xyz(:, iat)
             r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
             gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
             tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
             amat(jat, iat) = amat(jat, iat) + tmp
           else
             tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
             amat(iat, iat) = amat(iat, iat) + tmp
           endif
         endif
         if (irow .le. mol%nat/2 .and. irow .le. jcol) then
           jat = mol%nat - irow + 1
           iat = mol%nat - jcol + 1
           izp = mol%id(iat)
           jzp = mol%id(jat)
           if (iat.ne.jat) then
             vec = mol%xyz(:, jat) - mol%xyz(:, iat)
             r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
             gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
             tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
             amat(jat, iat) = amat(jat, iat) + tmp
           else
             tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
             amat(iat, iat) = amat(iat, iat) + tmp
           endif
         endif
      end do
   end do
   !$omp end do
   !$omp end parallel

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self) &
   !$omp private(iat, jat, jcol, irow)
   !$omp do collapse(2) schedule(static)
   do jcol = 1, mol%nat
      do irow = 1, (mol%nat+1)/2
         if (irow .ge. jcol) then
           jat = irow
           iat = jcol
           amat(iat, jat) = amat(jat, iat)
         endif
         if (irow .le. mol%nat/2 .and. irow .le. jcol) then
           jat = mol%nat - irow + 1
           iat = mol%nat - jcol + 1
           amat(iat, jat) = amat(jat, iat)
         endif
      end do
   end do
   !$omp end do
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_0d_4d

subroutine get_amat_3d(self, mol, wsc, alpha, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   ! Thread-private array for reduction
   real(wp), allocatable :: amat_local(:, :)

   amat(:, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel default(none) &
   !$omp shared(amat, mol, self, wsc, dtrans, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, gam, wsw, vec, dtmp, rtmp, amat_local)
   allocate(amat_local, source=amat)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat_local(jat, iat) = amat_local(jat, iat) + (dtmp + rtmp) * wsw
            amat_local(iat, jat) = amat_local(iat, jat) + (dtmp + rtmp) * wsw
         end do
      end do

      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
         amat_local(iat, iat) = amat_local(iat, iat) + (dtmp + rtmp) * wsw
      end do

      dtmp = self%eta(izp) + sqrt2pi / self%rad(izp) - 2 * alpha / sqrtpi
      amat_local(iat, iat) = amat_local(iat, iat) + dtmp
   end do
   !$omp end do
   !$omp critical (get_amat_3d_)
   amat(:, :) = amat(:, :) + amat_local(:, :)
   !$omp end critical (get_amat_3d_)
   deallocate(amat_local)
   !$omp end parallel

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = erf(gam * r1) / r1 - erf(alp * r1) / r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, tmp

   amat = 0.0_wp
   fac = 4 * pi / vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      tmp = cos(dot_product(rij, vec)) * fac * exp(-0.25_wp * g2 / (alp * alp)) / g2
      amat = amat + tmp
   end do

end subroutine get_amat_rec_3d

subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_damat_3d(self, mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace)
   else
      call get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   end if
end subroutine get_coulomb_derivs

subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: atrace_local(:, :)
   real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(atrace, dadr, dadL, mol, self, qvec) &
   !$omp private(iat, izp, jat, jzp, gam, r2, vec, dG, dS, dtmp, arg) &
   !$omp private(atrace_local, dadr_local, dadL_local)
   allocate(atrace_local, source=atrace)
   allocate(dadr_local, source=dadr)
   allocate(dadL_local, source=dadL)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         arg = gam * gam * r2
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * sqrt(r2))
         dG = dtmp * vec
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)
         atrace_local(:, iat) = -dG * qvec(jat) + atrace_local(:, iat)
         atrace_local(:, jat) = +dG * qvec(iat) + atrace_local(:, jat)
         dadr_local(:, iat, jat) = -dG * qvec(iat)
         dadr_local(:, jat, iat) = +dG * qvec(jat)
         dadL_local(:, :, jat) = +dS * qvec(iat) + dadL_local(:, :, jat)
         dadL_local(:, :, iat) = +dS * qvec(jat) + dadL_local(:, :, iat)
      end do
   end do
   !$omp end do
   !$omp critical (get_damat_0d_)
   atrace(:, :) = atrace(:, :) + atrace_local(:, :)
   dadr(:, :, :) = dadr(:, :, :) + dadr_local(:, :, :)
   dadL(:, :, :) = dadL(:, :, :) + dadL_local(:, :, :)
   !$omp end critical (get_damat_0d_)
   deallocate(dadL_local, dadr_local, atrace_local)
   !$omp end parallel

end subroutine get_damat_0d

subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vol, gam, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: atrace_local(:, :)
   real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel default(none) &
   !$omp shared(mol, self, wsc, alpha, vol, dtrans, rtrans, qvec) &
   !$omp shared(atrace, dadr, dadL) &
   !$omp private(iat, izp, jat, jzp, img, gam, wsw, vec, dG, dS) &
   !$omp private(dGr, dSr, dGd, dSd, atrace_local, dadr_local, dadL_local)
   allocate(atrace_local, source=atrace)
   allocate(dadr_local, source=dadr)
   allocate(dadL_local, source=dadL)
   !$omp do schedule(runtime)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace_local(:, iat) = -dG * qvec(jat) + atrace_local(:, iat)
         atrace_local(:, jat) = +dG * qvec(iat) + atrace_local(:, jat)
         dadr_local(:, iat, jat) = -dG * qvec(iat) + dadr_local(:, iat, jat)
         dadr_local(:, jat, iat) = +dG * qvec(jat) + dadr_local(:, jat, iat)
         dadL_local(:, :, jat) = +dS * qvec(iat) + dadL_local(:, :, jat)
         dadL_local(:, :, iat) = +dS * qvec(jat) + dadL_local(:, :, iat)
      end do

      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
         dS = dS + (dSd + dSr) * wsw
      end do
      dadL_local(:, :, iat) = +dS * qvec(iat) + dadL_local(:, :, iat)
   end do
   !$omp end do
   !$omp critical (get_damat_3d_)
   atrace(:, :) = atrace(:, :) + atrace_local(:, :)
   dadr(:, :, :) = dadr(:, :, :) + dadr_local(:, :, :)
   dadL(:, :, :) = dadL(:, :, :) + dadL_local(:, :, :)
   !$omp end critical (get_damat_3d_)
   deallocate(dadL_local, dadr_local, atrace_local)
   !$omp end parallel

end subroutine get_damat_3d

subroutine get_damat_dir_3d(rij, gam, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, gam2, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   gam2 = gam * gam
   alp2 = alp * alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1 * r1
      gtmp = +2 * gam * exp(-r2 * gam2) / (sqrtpi * r2) - erf(r1 * gam) / (r2 * r1)
      atmp = -2 * alp * exp(-r2 * alp2) / (sqrtpi * r2) + erf(r1 * alp) / (r2 * r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

subroutine get_damat_rec_3d(rij, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, etmp, dtmp, alp2
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4 * pi / vol
   alp2 = alp * alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      etmp = fac * exp(-0.25_wp * g2 / alp2) / g2
      dtmp = -sin(gv) * etmp
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + etmp * cos(gv) &
                  & * ((2.0_wp / g2 + 0.5_wp / alp2) * spread(vec, 1, 3) * spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_rec_3d

!> Inspect cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the cache
   type(cache_container), target, intent(inout) :: cache
   !> Reference to the cache
   type(eeq_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not. allocated(cache%raw)) then
      block
         type(eeq_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the cache
   type(cache_container), target, intent(inout) :: cache
   !> Reference to the cache
   type(eeq_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(eeq_cache)
      ptr => target
   end select
end subroutine view

end module multicharge_model_eeq
