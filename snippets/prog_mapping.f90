program prog_mapping
  implicit none
  integer :: nat
  integer :: iat, jat
  integer :: irow, jcol
  integer :: idx
  character(len=32) :: arg
  call get_command_argument(1,arg)
  read(arg,*)nat
  write(*,'(2i8)')0,nat
  idx = 0
  do jcol = 1, nat
    do irow = 1, (nat+1)/2
      if (irow .ge. jcol) then
        jat = irow
        iat = jcol
        idx = jat*(jat-1)/2 + iat
        write(*,'(3i8)')idx,iat,jat
      endif
      if (irow .le. nat/2 .and. irow .le. jcol) then
        jat = nat - irow + 1
        iat = nat - jcol + 1
        idx = jat*(jat-1)/2 + iat
        write(*,'(3i8)')idx,iat,jat
      endif
    enddo
  enddo
end program prog_mapping
