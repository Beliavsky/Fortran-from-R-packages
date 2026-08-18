program test_bootstrap
  use quantreg, only : dp, rq_bootstrap_xy
  implicit none
  real(dp) :: x(9,1), y(9), b(1,10)
  integer :: info, i
  x=1.0_dp
  do i=1,9
    y(i)=real(i,dp)
  end do
  call rq_bootstrap_xy(x,y,0.5_dp,10,b,info,1234)
  if (info /= 0) error stop 'bootstrap info'
  if (any(b < 1.0_dp) .or. any(b > 9.0_dp)) error stop 'bootstrap range'
  print *, 'test_bootstrap: PASS'
end program
