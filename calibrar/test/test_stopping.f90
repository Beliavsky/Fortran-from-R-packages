program test_stopping
  use calibrar, only : dp, smooth_stop2, smooth_stop3, smooth_stop4, n_stop
  implicit none
  real(dp) :: x(40)
  logical :: b(30)
  x=1.0_dp;b=.true.
  if(.not.smooth_stop2(x,1.0e-8_dp,5)) error stop "smooth_stop2 failed"
  if(.not.smooth_stop3(x,1.0e-8_dp,5)) error stop "smooth_stop3 failed"
  if(.not.smooth_stop4(x,1.0e-8_dp,5)) error stop "smooth_stop4 failed"
  if(.not.n_stop(b,5)) error stop "n_stop failed"
  print *, "PASS test_stopping"
end program test_stopping
