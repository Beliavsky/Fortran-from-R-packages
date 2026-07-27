! SPDX-License-Identifier: GPL-3.0-only
program test_portfolios
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
  use hdshop, only: dp, traditional_portfolio, shrinkage_mv_portfolio, &
    shrinkage_gmv_portfolio, portfolio_result
  implicit none
  real(dp) :: x(3,8), target(3), xwide(5,4), inf
  type(portfolio_result) :: trad,mv,gmv,wide
  x = reshape([ &
    0.01_dp,0.005_dp,-0.002_dp, 0.02_dp,-0.01_dp,0.008_dp, &
   -0.01_dp,0.015_dp,0.012_dp, 0.03_dp,0.02_dp,-0.004_dp, &
    0.00_dp,-0.005_dp,0.016_dp, 0.015_dp,0.01_dp,0.006_dp, &
   -0.005_dp,0.012_dp,0.011_dp, 0.025_dp,0.018_dp,0.003_dp], [3,8])
  target=1.0_dp/3.0_dp
  trad=traditional_portfolio(x,3.0_dp);if(.not.trad%ok)error stop 1
  call assert_vec(trad%weights,[3.82527344391138_dp,0.127490970407655_dp, &
    -2.95276441431900_dp],3.0e-11_dp)
  call assert_close(sum(trad%weights),1.0_dp,1.0e-12_dp)
  mv=shrinkage_mv_portfolio(x,3.0_dp,target,0.05_dp);if(.not.mv%ok)error stop 1
  call assert_close(mv%alpha,-4.94380513227976_dp,5.0e-11_dp)
  if(.not.allocated(mv%weight_intervals))error stop 1
  gmv=shrinkage_gmv_portfolio(x,target,0.05_dp);if(.not.gmv%ok)error stop 1
  call assert_close(gmv%alpha,0.677712128006628_dp,5.0e-12_dp)
  call assert_vec(gmv%weights,[0.258165353225274_dp,0.230255479428008_dp, &
    0.511579167346718_dp],5.0e-12_dp)
  xwide=reshape([0.01_dp,0.02_dp,0.00_dp,-0.01_dp,0.03_dp, &
    0.02_dp,0.01_dp,-0.01_dp,0.00_dp,0.04_dp, &
    0.00_dp,-0.01_dp,0.02_dp,0.01_dp,0.01_dp, &
    0.03_dp,0.00_dp,0.01_dp,0.02_dp,-0.02_dp],[5,4])
  inf=ieee_value(1.0_dp,ieee_positive_inf)
  wide=shrinkage_gmv_portfolio(xwide,[0.2_dp,0.2_dp,0.2_dp,0.2_dp,0.2_dp])
  if(.not.wide%ok .or. abs(sum(wide%weights)-1.0_dp)>1.0e-9_dp)error stop 1
  wide=traditional_portfolio(xwide,inf)
  if(.not.wide%ok .or. abs(sum(wide%weights)-1.0_dp)>1.0e-9_dp)error stop 1
  print '(a)', 'test_portfolios: PASS'
contains
  subroutine assert_close(a,b,tol)
    real(dp),intent(in)::a,b,tol
    if(abs(a-b)>tol)then;print *,'mismatch:',a,b,abs(a-b);error stop 1;endif
  end subroutine
  subroutine assert_vec(a,b,tol)
    real(dp),intent(in)::a(:),b(:),tol
    if(maxval(abs(a-b))>tol)then;print *,'vector mismatch:',a,b;error stop 1;endif
  end subroutine
end program test_portfolios
