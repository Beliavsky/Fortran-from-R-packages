! SPDX-License-Identifier: GPL-3.0-only
program example_mean_variance
  use portfolio_analytics
  implicit none
  real(dp) :: returns(8,3)
  type(portfolio_constraints) :: constraints
  type(portfolio_options) :: options
  type(portfolio_result) :: result
  integer :: i

  do i=1,8
    returns(i,1)=0.004_dp+0.008_dp*sin(real(i,dp))
    returns(i,2)=0.006_dp+0.012_dp*cos(real(i,dp))
    returns(i,3)=0.009_dp+0.018_dp*sin(0.5_dp*real(i,dp))
  end do
  call initialize_constraints(constraints,3,max_weight=[0.7_dp,0.7_dp,0.7_dp])
  options%objective=obj_min_variance
  call optimize_portfolio(returns,constraints,options,result)
  write(*,'(a,3f10.5)') 'weights:',result%weights
  write(*,'(a,f10.6)') 'expected return:',result%expected_return
  write(*,'(a,f10.6)') 'volatility:',result%risk
end program example_mean_variance
