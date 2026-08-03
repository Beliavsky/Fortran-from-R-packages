! SPDX-License-Identifier: GPL-3.0-only
program demo_portfolio_analytics
  use portfolio_analytics
  implicit none
  real(dp) :: returns(60,4)
  type(portfolio_constraints) :: constraints
  type(portfolio_options) :: options
  type(portfolio_result) :: result
  type(frontier_result) :: frontier
  integer :: t,i

  do t=1,60
    do i=1,4
      returns(t,i)=0.001_dp*real(i,dp)+0.006_dp*real(i,dp)*sin(0.17_dp*real(t*i,dp))
    end do
  end do
  call initialize_constraints(constraints,4,max_weight=[0.6_dp,0.6_dp,0.6_dp,0.6_dp])
  allocate(constraints%group_a(2,4),constraints%group_lower(2),constraints%group_upper(2))
  constraints%group_a=0.0_dp
  constraints%group_a(1,1:2)=1.0_dp
  constraints%group_a(2,3:4)=1.0_dp
  constraints%group_lower=[0.2_dp,0.2_dp]
  constraints%group_upper=[0.8_dp,0.8_dp]
  options%objective=obj_max_sharpe
  options%optimizer=opt_differential_evolution
  options%max_iterations=120
  options%population_size=40
  call optimize_portfolio(returns,constraints,options,result)
  write(*,'(a,4f9.4)') 'maximum-Sharpe weights:',result%weights
  write(*,'(a,f9.4)') 'Sharpe ratio:',result%sharpe
  options%objective=obj_min_variance
  call create_efficient_frontier(returns,constraints,options,5,frontier)
  write(*,'(a)') 'frontier: return volatility'
  do i=1,5
    write(*,'(2f12.6)') frontier%expected_return(i),frontier%risk(i)
  end do
end program demo_portfolio_analytics
