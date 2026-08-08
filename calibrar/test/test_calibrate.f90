program test_calibrate
  use calibrar, only : dp, optim_options, optim_result, calibrate, shifted_quadratic
  implicit none
  type(optim_options) :: op
  type(optim_result) :: res
  real(dp) :: x0(3),lo(3),hi(3)
  integer :: phases(3),reps(2)
  x0=[-2.0_dp,3.0_dp,-1.0_dp];lo=-5.0_dp;hi=5.0_dp;phases=[1,2,1];reps=[1,2]
  op%maxit=250;op%maxfeval=10000;op%reltol=1.0e-7_dp
  call calibrate(x0,shifted_quadratic,phases,res,"BFGS",lo,hi,op,replicates=reps)
  if(maxval(abs(res%par-1.5_dp))>2.0e-4_dp) error stop "phased calibrate failed"
  print *, "PASS test_calibrate"
end program test_calibrate
