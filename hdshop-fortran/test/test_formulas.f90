! SPDX-License-Identifier: GPL-3.0-only
program test_formulas
  use hdshop, only: dp, sigma_sample_estimator, row_means, alpha_star, &
    alpha_star_gmv, var_alpha_simple
  implicit none
  real(dp) :: x(3,8), b(3), mu(3)
  real(dp), allocatable :: s(:,:)
  x = reshape([ &
    0.01_dp,0.005_dp,-0.002_dp, 0.02_dp,-0.01_dp,0.008_dp, &
   -0.01_dp,0.015_dp,0.012_dp, 0.03_dp,0.02_dp,-0.004_dp, &
    0.00_dp,-0.005_dp,0.016_dp, 0.015_dp,0.01_dp,0.006_dp, &
   -0.005_dp,0.012_dp,0.011_dp, 0.025_dp,0.018_dp,0.003_dp], [3,8])
  s=sigma_sample_estimator(x);mu=row_means(x);b=1.0_dp/3.0_dp
  call assert_close(alpha_star(3.0_dp,mu,s,b,0.375_dp), &
    0.0402697175172202_dp,5.0e-13_dp)
  call assert_close(alpha_star_gmv(s,b,0.375_dp),0.813589537799146_dp,5.0e-13_dp)
  call assert_close(var_alpha_simple(s,b,8),0.179845456767587_dp,5.0e-13_dp)
  print '(a)', 'test_formulas: PASS'
contains
  subroutine assert_close(a,bv,tol)
    real(dp),intent(in)::a,bv,tol
    if(abs(a-bv)>tol)then;print *,'mismatch:',a,bv,abs(a-bv);error stop 1;endif
  end subroutine assert_close
end program test_formulas
