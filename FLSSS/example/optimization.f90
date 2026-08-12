program optimization
  use flsss_api
  implicit none
  type(knapsack_result) :: k
  type(gap_result) :: g
  real(dp) :: cost(2,3),profit(2,3),budget(2)
  k=mm_knapsack(2,[10.0_dp,8.0_dp,7.0_dp,6.0_dp], &
    reshape([2.0_dp,3.0_dp,4.0_dp,1.0_dp, 3.0_dp,2.0_dp,1.0_dp,5.0_dp],[4,2]), &
    [6.0_dp,6.0_dp])
  print '(a,f8.3)', 'knapsack profit: ',k%selection_profit
  cost=reshape([2.0_dp,4.0_dp,3.0_dp,2.0_dp,4.0_dp,1.0_dp],[2,3])
  profit=reshape([8.0_dp,7.0_dp,6.0_dp,9.0_dp,7.0_dp,5.0_dp],[2,3])
  budget=[6.0_dp,5.0_dp]
  g=gap_solve(cost,profit,budget)
  print '(a,f8.3)', 'GAP profit: ',g%total_profit_or_loss
  print '(a,*(i0,1x))','assignment: ',g%assignment
end program optimization
