! SPDX-License-Identifier: GPL-3.0-only
program benchmark_projection
  use portfoliooptim, only : dp, portfolio_result, portfolio_optim_projection, risk_mad
  implicit none
  real(dp) :: returns(4, 2), probabilities(4), a(2, 2), b(2)
  real(dp) :: lower(2), upper(2), benchmark(2)
  type(portfolio_result) :: result

  returns(1, :) = [-0.10_dp, 0.02_dp]
  returns(2, :) = [0.04_dp, -0.03_dp]
  returns(3, :) = [0.08_dp, 0.05_dp]
  returns(4, :) = [0.02_dp, 0.01_dp]
  probabilities = 0.25_dp
  a(1, :) = 1.0_dp
  a(2, :) = -1.0_dp
  b = [1.0_dp, -1.0_dp]
  lower = 0.0_dp
  upper = 1.0_dp
  benchmark = [0.70_dp, 0.30_dp]

  result = portfolio_optim_projection(returns, probabilities, 0.0_dp, &
    risk_mad, benchmark, 0.95_dp, a, b, lower, upper)
  print '(a,2(f11.7,1x))', 'benchmark: ', benchmark
  print '(a,2(f11.7,1x))', 'projected: ', result%theta
  print '(a,f11.7)', 'MAD:       ', result%mad
end program benchmark_projection
