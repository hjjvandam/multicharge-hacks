program prog_single
  integer(kind=8) :: nat
  integer(kind=8) :: iat, jat
  integer(kind=8) :: idx
  integer(kind=8) :: ntop
  character(len=32) :: arg
  call get_command_argument(1,arg)
  read(arg,*)nat
  write(*,'(2i8)')0,nat
  ntop = nat*(nat+1)/2
  do idx = 1, ntop
    jat = int((1.0d0 + sqrt(1.0d0 + 8.0d0 * real(idx-1, 8))) / 2.0d0)
    iat = idx - jat*(jat-1)/2
    write(*,'(3i8)')idx,iat,jat
  enddo
end program prog_single
