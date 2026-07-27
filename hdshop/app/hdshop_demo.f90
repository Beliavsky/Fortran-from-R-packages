! SPDX-License-Identifier: GPL-3.0-only
program hdshop_demo
  use hdshop, only: dp, shrinkage_gmv_portfolio, portfolio_result
  implicit none
  real(dp)::x(3,8),target(3)
  type(portfolio_result)::portfolio
  x=reshape([0.01_dp,0.005_dp,-0.002_dp,0.02_dp,-0.01_dp,0.008_dp, &
    -0.01_dp,0.015_dp,0.012_dp,0.03_dp,0.02_dp,-0.004_dp, &
    0.0_dp,-0.005_dp,0.016_dp,0.015_dp,0.01_dp,0.006_dp, &
    -0.005_dp,0.012_dp,0.011_dp,0.025_dp,0.018_dp,0.003_dp],[3,8])
  target=1.0_dp/3.0_dp;portfolio=shrinkage_gmv_portfolio(x,target)
  print '(a,3f12.6)', 'shrinkage GMV weights: ',portfolio%weights
  print '(a,f10.6)', 'alpha: ',portfolio%alpha
end program hdshop_demo
