program portfolio_pl
   use quarks
   implicit none
   real(dp) :: asset_returns(4,3), weights(3)
   type(pl_result) :: portfolio

   asset_returns(1,:) = [0.010_dp, 0.004_dp, -0.002_dp]
   asset_returns(2,:) = [-0.006_dp, 0.008_dp, 0.003_dp]
   asset_returns(3,:) = [0.002_dp, -0.004_dp, 0.009_dp]
   asset_returns(4,:) = [0.005_dp, 0.001_dp, -0.003_dp]
   weights = [0.50_dp, 0.30_dp, 0.20_dp]
   portfolio = plop(asset_returns, weights, approximation=1)

   print '(a)', 'portfolio returns:'
   print '(4f12.6)', portfolio%pl
end program portfolio_pl
