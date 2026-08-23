program example_kde
  use kernsmooth_mod
  implicit none
  real(dp) :: x(5)=[-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp]
  real(dp) :: grid(101), dens(101)
  integer :: i
  call bkde(x,0.4_dp,grid,dens)
  do i=1,size(grid),10
    write(*,'(2f12.6)') grid(i),dens(i)
  end do
end program
