! SPDX-License-Identifier: GPL-3.0-only
program test_inference_frontier
  use hdshop, only: dp, test_mvsp, mvsp_test_result, bayesian_frontier, &
    frontier_result, traditional_portfolio, portfolio_result
  implicit none
  real(dp) :: x(3,12), target(3), weights(3,2), rj
  integer :: j
  type(mvsp_test_result) :: tst
  type(frontier_result) :: fr
  type(portfolio_result) :: p
  do j=1,12
    rj=real(j,dp)
    x(1,j)=0.01_dp+0.02_dp*sin(0.7_dp*rj)+0.003_dp*cos(0.2_dp*rj)
    x(2,j)=0.008_dp+0.015_dp*cos(0.5_dp*rj)-0.002_dp*sin(0.3_dp*rj)
    x(3,j)=0.006_dp+0.012_dp*sin(0.9_dp*rj+1.0_dp)+0.004_dp*cos(0.4_dp*rj)
  end do
  target=1.0_dp/3.0_dp;tst=test_mvsp(3.0_dp,x,target)
  if(.not.tst%ok .or. tst%alpha_sd<=0.0_dp .or. tst%p_value<0.0_dp .or. tst%p_value>1.0_dp)error stop 1
  p=traditional_portfolio(x,3.0_dp);weights(:,1)=target;weights(:,2)=p%weights
  fr=bayesian_frontier(x,weights,25)
  if(.not.fr%ok .or. size(fr%frontier_sd)/=25)error stop 1
  if(any(fr%frontier_sd(2:)<fr%frontier_sd(:24)))error stop 1
  if(any(fr%frontier_return(2:)<fr%frontier_return(:24)))error stop 1
  print '(a)', 'test_inference_frontier: PASS'
end program test_inference_frontier
