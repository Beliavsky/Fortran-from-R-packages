program returns_and_rebalance
   use pmwr, only : dp, rebalance_result, simple_returns, rebalance_portfolio
   implicit none
   real(dp), allocatable :: r(:)
   type(rebalance_result) :: rb

   call simple_returns([100.0_dp, 104.0_dp, 101.0_dp, 108.0_dp], r)
   print '(a,*(f9.5,1x))', 'Returns: ', r

   call rebalance_portfolio(current=[10.0_dp, 20.0_dp], &
                            target=[0.60_dp, 0.40_dp], &
                            price=[12.0_dp, 5.0_dp], result=rb)
   print '(a,f10.2)', 'Notional: ', rb%notional
   print '(a,*(f10.3,1x))', 'Orders:   ', rb%difference
   print '(a,f10.2)', 'Turnover: ', rb%turnover
end program returns_and_rebalance
