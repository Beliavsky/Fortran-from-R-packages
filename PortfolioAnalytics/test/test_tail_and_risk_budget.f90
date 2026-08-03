! SPDX-License-Identifier: GPL-3.0-only
program test_tail_and_risk_budget
  use portfolio_analytics
  use test_support
  implicit none
  real(dp) :: r(10,2),mu(2),sigma(2,2),budgets(2)
  type(portfolio_constraints) :: c
  type(portfolio_options) :: options
  type(portfolio_result) :: result

  r(:,1)=[0.01_dp,0.01_dp,0.00_dp,-0.01_dp,0.01_dp,0.00_dp,-0.01_dp,0.01_dp,0.00_dp,-0.01_dp]
  r(:,2)=[0.08_dp,0.07_dp,0.06_dp,0.05_dp,0.04_dp,0.03_dp,-0.20_dp,-0.15_dp,-0.10_dp,-0.08_dp]
  call sample_moments(r,mu,sigma)
  call initialize_constraints(c,2)
  options%objective=obj_min_es
  options%optimizer=opt_differential_evolution
  options%max_iterations=80
  options%population_size=30
  options%seed=42
  options%alpha=0.2_dp
  call optimize_portfolio(r,c,options,result,mu,sigma)
  call assert_true(result%feasible,'ES portfolio feasible')
  call assert_true(result%weights(1)>0.95_dp,'ES selects safer asset')

  sigma=reshape([0.04_dp,0.01_dp,0.01_dp,0.04_dp],[2,2])
  mu=0.0_dp
  budgets=[0.5_dp,0.5_dp]
  options%objective=obj_risk_parity
  options%max_iterations=80
  options%seed=7
  call optimize_portfolio(r,c,options,result,mu,sigma,budgets)
  call assert_close(result%weights(1),0.5_dp,2.0e-4_dp,'risk parity weight one')
  call assert_close(result%weights(2),0.5_dp,2.0e-4_dp,'risk parity weight two')
  print '(a)','test_tail_and_risk_budget: PASS'
end program test_tail_and_risk_budget
