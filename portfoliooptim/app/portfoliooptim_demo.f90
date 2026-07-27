! SPDX-License-Identifier: GPL-3.0-only
program portfoliooptim_demo
  use portfoliooptim, only : dp, portfolio_result, bdportfolio_optim, &
    portfolio_optim_projection, risk_mad
  implicit none
  real(dp) :: returns(5, 2), probabilities(5), a(2, 2), b(2)
  real(dp) :: lower(2), upper(2), benchmark(2)
  type(portfolio_result) :: benders, projected

  returns(1, :) = [-0.10_dp, 0.02_dp]
  returns(2, :) = [0.04_dp, -0.03_dp]
  returns(3, :) = [0.08_dp, 0.05_dp]
  returns(4, :) = [0.02_dp, 0.01_dp]
  returns(5, :) = [0.06_dp, 0.03_dp]
  probabilities = 0.20_dp
  a(1, :) = 1.0_dp
  a(2, :) = -1.0_dp
  b = [1.0_dp, -1.0_dp]
  lower = 0.0_dp
  upper = 1.0_dp
  benchmark = [0.60_dp, 0.40_dp]

  benders = bdportfolio_optim(returns, probabilities, 0.015_dp, risk_mad, &
    0.80_dp, a, b, lower, upper)
  projected = portfolio_optim_projection(returns, probabilities, 0.015_dp, &
    risk_mad, benchmark, 0.80_dp, a, b, lower, upper)

  print '(a)', 'Benders MAD portfolio'
  print '(a,*(f11.6,1x))', 'weights: ', benders%theta
  print '(a,f11.6)', 'mean:    ', benders%mu
  print '(a,f11.6)', 'MAD:     ', benders%mad
  print '(a)', ''
  print '(a)', 'Closest risk-optimal portfolio to benchmark'
  print '(a,*(f11.6,1x))', 'weights: ', projected%theta
  print '(a,*(f11.6,1x))', 'target:  ', benchmark
end program portfoliooptim_demo
