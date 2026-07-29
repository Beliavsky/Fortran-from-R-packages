! SPDX-License-Identifier: GPL-3.0-or-later
program drawdown_portfolios
  use frapo
  implicit none

  real(dp) :: prices(8, 3)
  type(portfolio_result) :: maxdd, cdar

  prices(:, 1) = [100.0_dp, 102.0_dp, 101.0_dp, 104.0_dp, &
                  106.0_dp, 105.0_dp, 108.0_dp, 110.0_dp]
  prices(:, 2) = [100.0_dp, 100.5_dp, 101.0_dp, 101.5_dp, &
                  102.0_dp, 102.5_dp, 103.0_dp, 103.5_dp]
  prices(:, 3) = [100.0_dp, 99.0_dp, 98.0_dp, 100.0_dp, &
                  101.0_dp, 103.0_dp, 102.0_dp, 104.0_dp]

  maxdd = pmaxdd(prices, max_drawdown=0.08_dp)
  cdar = pcdar(prices, alpha=0.90_dp, bound=0.08_dp)

  write(*, '(a,*(f10.6,1x))') 'Maximum-DD weights: ', maxdd%weights
  write(*, '(a,f10.6)') 'Realized maximum drawdown: ', maxdd%risk_value
  write(*, '(a,*(f10.6,1x))') 'CDaR-constrained weights: ', cdar%weights
  write(*, '(a,f10.6)') 'Realized CDaR: ', cdar%risk_value
end program drawdown_portfolios
