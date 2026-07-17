program prog_triangle
  integer :: nat
  integer :: iat, jat
  integer :: idx, idy
  character(len=32) :: arg
  call get_command_argument(1,arg)
  read(arg,*)nat
  write(*,'(2i8)')0,nat
  idx = 0
  do jat = 1, nat
    do iat = 1, jat
      idx = jat*(jat-1)/2+iat
      write(*,'(3i8)')idx,iat,jat
    enddo
  enddo
end program prog_triangle
