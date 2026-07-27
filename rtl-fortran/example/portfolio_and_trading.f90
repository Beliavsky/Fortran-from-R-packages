! SPDX-License-Identifier: MIT
program portfolio_and_trading
  use rtl, only: dp, frontier_result, refinery_result, strategy_result, trade_stats_result
  use rtl, only: efficient_frontier_statistics, refinery_lp, moving_average_strategy, trade_stats
  implicit none

  type(frontier_result) :: frontier
  type(refinery_result) :: refinery
  type(strategy_result) :: strategy
  type(trade_stats_result) :: statistics
  real(dp) :: covariance(3, 3), yields(2, 2)
  real(dp) :: open_price(8), close_price(8)

  covariance = reshape([0.04_dp, 0.01_dp, 0.00_dp, &
                        0.01_dp, 0.09_dp, 0.02_dp, &
                        0.00_dp, 0.02_dp, 0.06_dp], shape(covariance))
  frontier = efficient_frontier_statistics(1000, [0.08_dp, 0.12_dp, 0.10_dp], covariance, seed=123)

  yields = reshape([0.5_dp, 0.2_dp, 0.3_dp, 0.6_dp], shape(yields))
  refinery = refinery_lp([40.0_dp, 35.0_dp], [2.0_dp, 3.0_dp], &
    [100.0_dp, 80.0_dp], yields, [100.0_dp, 120.0_dp])

  open_price = [10.0_dp, 10.2_dp, 10.4_dp, 10.6_dp, 10.8_dp, 11.0_dp, 10.8_dp, 10.6_dp]
  close_price = [10.1_dp, 10.3_dp, 10.5_dp, 10.7_dp, 10.9_dp, 10.7_dp, 10.5_dp, 10.4_dp]
  strategy = moving_average_strategy(open_price, close_price, 2, 3)
  statistics = trade_stats(strategy%strategy_return)

  print '(a,3f10.4)', 'Maximum-Sharpe weights: ', &
    frontier%weights(frontier%maximum_sharpe_index, :)
  print '(a,f12.4)', 'Refinery profit: ', refinery%profit
  print '(a,f12.6)', 'Strategy cumulative return: ', statistics%cumulative_return
end program portfolio_and_trading
