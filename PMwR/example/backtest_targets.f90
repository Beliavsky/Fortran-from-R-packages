program backtest_targets
   use pmwr, only : dp, backtest_result, run_backtest
   implicit none
   type(backtest_result) :: result
   real(dp) :: prices(6,2), targets(6,2)

   prices(:,1) = [100.0_dp, 102.0_dp, 105.0_dp, 103.0_dp, 108.0_dp, 110.0_dp]
   prices(:,2) = [50.0_dp, 49.0_dp, 51.0_dp, 53.0_dp, 52.0_dp, 54.0_dp]
   targets = 0.0_dp
   targets(1:2,:) = spread([0.60_dp, 0.40_dp], 1, 2)
   targets(3:4,:) = spread([0.30_dp, 0.70_dp], 1, 2)
   targets(5:6,:) = spread([0.50_dp, 0.50_dp], 1, 2)

   call run_backtest(prices, targets, result, initial_cash=100000.0_dp, &
                     convert_weights=.true., transaction_cost=[0.001_dp], lag=1)
   print '(a,f12.2)', 'Initial wealth: ', result%initial_wealth
   print '(a,f12.2)', 'Final wealth:   ', result%wealth(size(result%wealth))
   print '(a,f12.2)', 'Trading costs:  ', result%cumulative_cost(size(result%cumulative_cost))
   print '(a,i0)', 'Journal rows:   ', result%journal%n
end program backtest_targets
