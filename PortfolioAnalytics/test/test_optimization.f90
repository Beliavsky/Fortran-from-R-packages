! SPDX-License-Identifier: GPL-3.0-only
program test_optimization
  use portfolio_analytics
  use test_support
  implicit none
  real(dp) :: r(6,2),mu(2),sigma(2,2),w(2)
  type(portfolio_constraints) :: c
  type(portfolio_options) :: options
  type(portfolio_result) :: result
  logical :: ok

  r(1,:)=[0.01_dp,0.03_dp]
  r(2,:)=[-0.01_dp,0.02_dp]
  r(3,:)=[0.02_dp,-0.01_dp]
  r(4,:)=[0.00_dp,0.04_dp]
  r(5,:)=[0.01_dp,0.00_dp]
  r(6,:)=[-0.02_dp,0.01_dp]
  mu=[0.05_dp,0.08_dp]
  sigma=reshape([0.04_dp,0.0_dp,0.0_dp,0.09_dp],[2,2])
  call initialize_constraints(c,2)
  options%objective=obj_min_variance
  options%optimizer=opt_projected_gradient
  options%max_iterations=500
  options%tolerance=1.0e-9_dp
  call optimize_portfolio(r,c,options,result,mu,sigma)
  call assert_true(result%feasible,'minimum variance feasible')
  call assert_close(result%weights(1),9.0_dp/13.0_dp,5.0e-5_dp,'GMV weight one')
  call assert_close(result%weights(2),4.0_dp/13.0_dp,5.0e-5_dp,'GMV weight two')

  options%objective=obj_max_return
  call optimize_portfolio(r,c,options,result,mu,sigma)
  call assert_true(result%weights(2)>0.9999_dp,'maximum return corner')

  options%objective=obj_quadratic_utility
  options%risk_aversion=4.0_dp
  call optimize_portfolio(r,c,options,result,mu,sigma)
  call assert_true(result%weights(1)>0.4_dp .and. result%weights(1)<0.9_dp,'quadratic utility interior')

  call equal_weight_portfolio(c,w,ok)
  call assert_all_close(w,[0.5_dp,0.5_dp],1.0e-12_dp,'equal weight')
  call inverse_volatility_portfolio(sigma,c,w,ok)
  call assert_close(w(1),0.6_dp,1.0e-12_dp,'inverse volatility weight')
  print '(a)','test_optimization: PASS'
end program test_optimization
