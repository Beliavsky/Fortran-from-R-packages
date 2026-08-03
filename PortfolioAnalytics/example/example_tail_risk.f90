! SPDX-License-Identifier: GPL-3.0-only
program example_tail_risk
  use portfolio_analytics
  implicit none
  real(dp) :: returns(10,2)
  type(portfolio_constraints) :: constraints
  type(portfolio_options) :: options
  type(portfolio_result) :: result

  returns(:,1)=[0.01_dp,0.01_dp,0.00_dp,-0.01_dp,0.01_dp,0.00_dp,-0.01_dp,0.01_dp,0.00_dp,-0.01_dp]
  returns(:,2)=[0.08_dp,0.07_dp,0.06_dp,0.05_dp,0.04_dp,0.03_dp,-0.20_dp,-0.15_dp,-0.10_dp,-0.08_dp]
  call initialize_constraints(constraints,2)
  options%objective=obj_min_es
  options%optimizer=opt_differential_evolution
  options%alpha=0.20_dp
  options%max_iterations=100
  options%population_size=30
  call optimize_portfolio(returns,constraints,options,result)
  write(*,'(a,2f10.5)') 'minimum-ES weights:',result%weights
  write(*,'(a,f10.6)') 'expected shortfall:',result%expected_shortfall
end program example_tail_risk
