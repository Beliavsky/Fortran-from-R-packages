program density_example
  use locfit
  implicit none
  real(dp) :: x(10), grid(7), den(7)
  integer :: st(7), i

  x=[-1.6_dp,-1.1_dp,-0.7_dp,-0.25_dp,0.0_dp,0.2_dp,0.55_dp,0.9_dp,1.2_dp,1.7_dp]
  do i=1,size(grid)
    grid(i)=-1.5_dp+0.5_dp*real(i-1,dp)
  end do

  call local_density_1d(x,grid,0.65_dp,den,degree=2,ker=wgaus,link=llog,status=st)
  write(*,'(a)') '       x          density   status'
  do i=1,size(grid)
    write(*,'(2f15.8,2x,i0)') grid(i),den(i),st(i)
  end do
end program density_example
