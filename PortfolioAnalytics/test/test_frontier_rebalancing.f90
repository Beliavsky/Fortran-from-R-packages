! SPDX-License-Identifier: GPL-3.0-only
program test_frontier_rebalancing
  use portfolio_analytics
  use test_support
  implicit none
  real(dp) :: r(20,3),costs(3)
  logical :: rebalance(20)
  type(portfolio_constraints) :: c
  type(portfolio_options) :: options
  type(frontier_result) :: frontier
  type(rebalancing_result) :: backtest
  integer :: t

  do t=1,20
    r(t,1)=0.002_dp+0.004_dp*sin(real(t,dp))
    r(t,2)=0.004_dp+0.006_dp*cos(real(t,dp))
    r(t,3)=0.007_dp+0.010_dp*sin(0.5_dp*real(t,dp))
  end do
  call initialize_constraints(c,3,max_weight=[0.8_dp,0.8_dp,0.8_dp])
  options%objective=obj_min_variance
  options%optimizer=opt_projected_gradient
  options%max_iterations=120
  options%population_size=30
  options%seed=19
  call create_efficient_frontier(r,c,options,4,frontier)
  call assert_true(all(frontier%feasible),'efficient frontier feasibility')
  call assert_true(all(frontier%expected_return(2:4)>=frontier%expected_return(1:3)-1.0e-7_dp), &
                   'frontier return monotonicity')

  rebalance=.false.
  rebalance(5)=.true.
  rebalance(10)=.true.
  rebalance(15)=.true.
  rebalance(20)=.true.
  costs=0.001_dp
  options%objective=obj_max_return
  options%optimizer=opt_projected_gradient
  call optimize_rebalancing(r,rebalance,5,c,options,costs,backtest)
  call assert_true(size(backtest%weights,2)==20,'rebalancing dimensions')
  call assert_true(all(abs(sum(backtest%weights,dim=1)-1.0_dp)<1.0e-8_dp),'rebalanced weight sums')
  call assert_true(all(backtest%wealth>0.0_dp),'positive wealth path')
  call assert_true(sum(backtest%transaction_cost)>=0.0_dp,'nonnegative transaction costs')
  print '(a)','test_frontier_rebalancing: PASS'
end program test_frontier_rebalancing
