! SPDX-License-Identifier: GPL-3.0-only
program high_dimensional_portfolio
  use hdshop, only: dp, shrinkage_gmv_portfolio, portfolio_result
  implicit none
  real(dp)::x(5,4),target(5)
  type(portfolio_result)::p
  x=reshape([0.01_dp,0.02_dp,0.00_dp,-0.01_dp,0.03_dp, &
    0.02_dp,0.01_dp,-0.01_dp,0.00_dp,0.04_dp, &
    0.00_dp,-0.01_dp,0.02_dp,0.01_dp,0.01_dp, &
    0.03_dp,0.00_dp,0.01_dp,0.02_dp,-0.02_dp],[5,4])
  target=0.2_dp;p=shrinkage_gmv_portfolio(x,target)
  print '(a,5f10.5)', 'p > n shrinkage weights: ',p%weights
end program high_dimensional_portfolio
