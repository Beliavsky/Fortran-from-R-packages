! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program test_returns_drawdowns
  use kinds_mod, only: dp
  use returns_mod
  use drawdown_mod
  use test_support_mod
  implicit none
  real(dp) :: prices(4), r0(3), expected(3), dd(4), rr(4,2), w(2)
  real(dp) :: levels(5), converted(2)
  real(dp), allocatable :: r(:)
  type(portfolio_result) :: pr
  type(drawdown_episode), allocatable :: episodes(:)
  integer :: n, ne

  prices=[100.0_dp,110.0_dp,99.0_dp,108.9_dp]
  call calculate_returns(prices,'discrete',r)
  expected=[0.1_dp,-0.1_dp,0.1_dp]
  call assert_vector_close(r,expected,1.0e-13_dp,'discrete returns')
  call calculate_returns(prices,'log',r)
  call assert_close(r(1),log(1.1_dp),1.0e-13_dp,'log return')
  r0=expected
  call assert_close(cumulative_return(r0),0.089_dp,1.0e-13_dp,'cumulative return')
  call assert_close(annualized_return(r0,3.0_dp),0.089_dp,1.0e-13_dp,'annualized return')

  call level_from_returns(r0,100.0_dp,levels)
  call assert_close(levels(4),108.9_dp,1.0e-12_dp,'return levels')
  call convert_return_frequency([0.01_dp,0.02_dp,-0.01_dp,0.03_dp],2,converted,n)
  call assert_true(n==2,'converted count')
  call assert_close(converted(1),0.0302_dp,1.0e-13_dp,'converted first')

  rr=reshape([0.10_dp,0.00_dp,-0.10_dp,0.05_dp, &
              0.00_dp,0.10_dp,0.05_dp,-0.05_dp],[4,2])
  w=[0.6_dp,0.4_dp]
  call portfolio_returns(rr,w,pr)
  call assert_close(pr%returns(1),0.06_dp,1.0e-13_dp,'portfolio first return')
  call assert_close(pr%contributions(1,1),0.06_dp,1.0e-13_dp,'portfolio contribution')
  call assert_close(sum(pr%end_weights(1,:)),1.0_dp,1.0e-13_dp,'portfolio end weights')
  call portfolio_returns(rr,w,pr,rebalance_every=1,transaction_cost=0.001_dp)
  call assert_true(all(pr%turnover>=0.0_dp),'portfolio turnover')

  call drawdown_series([0.10_dp,-0.20_dp,0.05_dp,0.20_dp],dd)
  call assert_close(dd(1),0.0_dp,1.0e-13_dp,'drawdown initial peak')
  call assert_close(dd(2),-0.20_dp,1.0e-13_dp,'drawdown trough')
  call assert_close(max_drawdown([0.10_dp,-0.20_dp,0.05_dp,0.20_dp]),0.20_dp,1.0e-13_dp,'max drawdown')
  call find_drawdowns([0.10_dp,-0.20_dp,0.05_dp,0.20_dp],episodes,ne)
  call assert_true(ne==1,'drawdown episodes')
  call assert_true(episodes(1)%trough_index==2,'drawdown trough index')
  call assert_true(pain_index([0.10_dp,-0.20_dp,0.05_dp,0.20_dp])>0.0_dp,'pain index')
  call assert_true(ulcer_index([0.10_dp,-0.20_dp,0.05_dp,0.20_dp])>0.0_dp,'ulcer index')
  call assert_true(conditional_drawdown_at_risk([0.10_dp,-0.20_dp,0.05_dp,0.20_dp],0.8_dp)>=0.20_dp-1.0e-12_dp,'cdar')
  write(*,'(a)')'Return, portfolio, and drawdown tests passed.'
end program test_returns_drawdowns
