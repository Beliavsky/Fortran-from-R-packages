! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program optimal_allocation
  use r4good_personal_finances
  implicit none
  real(dp) :: expected_returns(3), standard_deviations(3), correlations(3,3)
  type(portfolio_result) :: result

  expected_returns = [0.0472_dp, 0.0504_dp, 0.0275_dp]
  standard_deviations = [0.1588_dp, 0.1718_dp, 0.0562_dp]
  correlations = reshape([1.00_dp,0.87_dp,0.21_dp, 0.87_dp,1.00_dp,0.37_dp, &
    0.21_dp,0.37_dp,1.00_dp], [3,3])

  call optimize_portfolio(0.35_dp, expected_returns, standard_deviations, correlations, result)
  print '(a,3f10.5)', 'Optimal weights: ', result%total
  print '(a,f10.5)', 'Expected return: ', result%expected_return
  print '(a,f10.5)', 'Standard deviation: ', result%standard_deviation
end program optimal_allocation
