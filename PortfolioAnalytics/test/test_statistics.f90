! SPDX-License-Identifier: GPL-3.0-only
program test_statistics
  use portfolio_analytics
  use test_support
  implicit none
  real(dp) :: r(4,2),mu(2),sigma(2,2),w(2),rp(4),rc(2)
  real(dp),allocatable :: m3(:,:),m4(:,:)

  r(1,:)=[0.10_dp,0.00_dp]
  r(2,:)=[0.00_dp,0.10_dp]
  r(3,:)=[-0.10_dp,0.00_dp]
  r(4,:)=[0.00_dp,-0.10_dp]
  w=[0.5_dp,0.5_dp]
  call sample_moments(r,mu,sigma)
  call assert_all_close(mu,[0.0_dp,0.0_dp],1.0e-14_dp,'sample means')
  call assert_close(sigma(1,1),0.02_dp/3.0_dp,1.0e-14_dp,'variance one')
  call assert_close(sigma(1,2),0.0_dp,1.0e-14_dp,'zero covariance')
  call portfolio_returns(r,w,rp)
  call assert_all_close(rp,[0.05_dp,0.05_dp,-0.05_dp,-0.05_dp],1.0e-14_dp,'portfolio returns')
  call assert_close(portfolio_variance(w,sigma),1.0_dp/300.0_dp,1.0e-14_dp,'portfolio variance')
  call assert_close(historical_var(rp,0.25_dp),0.05_dp,1.0e-14_dp,'historical VaR')
  call assert_close(historical_es(rp,0.25_dp),0.05_dp,1.0e-14_dp,'historical ES')
  call assert_close(hhi(w),0.5_dp,1.0e-14_dp,'HHI')
  call assert_close(diversification(w),0.5_dp,1.0e-14_dp,'diversification')
  call assert_close(turnover(w,[1.0_dp,0.0_dp]),0.5_dp,1.0e-14_dp,'turnover')
  call risk_contributions(w,sigma,rc)
  call assert_close(rc(1),rc(2),1.0e-14_dp,'equal risk contributions')
  call assert_true(maximum_drawdown(rp)>0.09_dp,'drawdown positive')
  call assert_true(conditional_second_moment(rp,0.25_dp)>0.0_dp,'CSM positive')
  call assert_true(expected_quadratic_shortfall(rp,0.25_dp)>=0.0_dp,'EQS nonnegative')
  allocate(m3(2,4),m4(2,8))
  call sample_coskewness(r,m3)
  call sample_cokurtosis(r,m4)
  call assert_close(portfolio_skewness(w,m3),0.0_dp,1.0e-14_dp,'portfolio skewness')
  call assert_true(portfolio_kurtosis(w,m4)>0.0_dp,'portfolio fourth moment')
  print '(a)','test_statistics: PASS'
end program test_statistics
