! SPDX-License-Identifier: GPL-3.0-only
program example_rebalancing
  use portfolio_analytics
  implicit none
  real(dp) :: returns(24,3),costs(3)
  logical :: rebalance(24)
  type(portfolio_constraints) :: constraints
  type(portfolio_options) :: options
  type(rebalancing_result) :: result
  integer :: t

  do t=1,24
    returns(t,1)=0.002_dp+0.005_dp*sin(real(t,dp))
    returns(t,2)=0.004_dp+0.007_dp*cos(real(t,dp))
    returns(t,3)=0.006_dp+0.010_dp*sin(0.4_dp*real(t,dp))
  end do
  rebalance=.false.
  rebalance(6:24:6)=.true.
  costs=0.0005_dp
  call initialize_constraints(constraints,3,max_weight=[0.8_dp,0.8_dp,0.8_dp])
  options%objective=obj_max_return
  call optimize_rebalancing(returns,rebalance,6,constraints,options,costs,result)
  write(*,'(a,f10.5)') 'ending wealth:',result%wealth(24)
  write(*,'(a,3f10.5)') 'final weights:',result%weights(:,24)
end program example_rebalancing
