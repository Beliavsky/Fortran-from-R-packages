! SPDX-License-Identifier: GPL-3.0-only
program test_means
  use hdshop, only: dp, mean_bs, mean_js, mean_bop19, mean_shrink_result
  implicit none
  real(dp) :: x(3,8), target(3)
  type(mean_shrink_result) :: bs,js,bop
  x = reshape([ &
    0.01_dp,0.005_dp,-0.002_dp, 0.02_dp,-0.01_dp,0.008_dp, &
   -0.01_dp,0.015_dp,0.012_dp, 0.03_dp,0.02_dp,-0.004_dp, &
    0.00_dp,-0.005_dp,0.016_dp, 0.015_dp,0.01_dp,0.006_dp, &
   -0.005_dp,0.012_dp,0.011_dp, 0.025_dp,0.018_dp,0.003_dp], [3,8])
  bs=mean_bs(x);if(.not.bs%ok)error stop 1
  call assert_close(bs%alpha,0.968787239155097_dp,2.0e-13_dp)
  call assert_vec(bs%means,[0.0076584782028458_dp,0.0075804463007335_dp, &
    0.0075219223741493_dp],2.0e-13_dp)
  js=mean_js(x,0.01_dp);if(.not.js%ok)error stop 1
  call assert_close(js%alpha,0.124930993629704_dp,2.0e-13_dp)
  target=0.01_dp;bop=mean_bop19(x,target);if(.not.bop%ok)error stop 1
  call assert_close(bop%alpha,-11.7699933948299_dp,2.0e-11_dp)
  call assert_close(bop%beta,9.65782042396744_dp,2.0e-11_dp)
  print '(a)', 'test_means: PASS'
contains
  subroutine assert_close(a,b,tol)
    real(dp),intent(in)::a,b,tol
    if(abs(a-b)>tol)then;print *,'mismatch:',a,b,abs(a-b);error stop 1;endif
  end subroutine
  subroutine assert_vec(a,b,tol)
    real(dp),intent(in)::a(:),b(:),tol
    if(maxval(abs(a-b))>tol)then;print *,'vector mismatch:',a,b;error stop 1;endif
  end subroutine
end program test_means
