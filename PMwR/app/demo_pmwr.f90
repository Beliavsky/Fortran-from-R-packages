program demo_pmwr
   use pmwr, only : dp, backtest_result, nav_summary_result, run_backtest, summarize_nav, &
                    valid_isin, format_quote32
   implicit none
   type(backtest_result) :: bt
   type(nav_summary_result) :: nav
   real(dp) :: prices(5,2), targets(5,2)

   prices(:,1) = [100.0_dp, 103.0_dp, 101.0_dp, 107.0_dp, 109.0_dp]
   prices(:,2) = [80.0_dp, 79.0_dp, 82.0_dp, 81.0_dp, 84.0_dp]
   targets = 0.5_dp

   call run_backtest(prices, targets, bt, initial_cash=10000.0_dp, &
                     convert_weights=.true., transaction_cost=[0.0005_dp], lag=1)
   call summarize_nav(bt%wealth, nav, periods_per_year=252.0_dp)

   print '(a,f12.2)', 'Final wealth: ', bt%wealth(size(bt%wealth))
   print '(a,f8.3,a)', 'Total return: ', 100.0_dp*nav%total_return, '%'
   print '(a,l1)', 'Apple ISIN valid: ', valid_isin('US0378331005')
   print '(a,a)', 'Treasury quote: ', trim(format_quote32(99.5078125_dp))
end program demo_pmwr
