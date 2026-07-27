! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
program test_backtest_monitor
  use fportfolio_kinds, only: dp
  use fportfolio_types
  use fportfolio_optimization, only: initialize_constraints
  use fportfolio_backtest
  use fportfolio_monitor
  implicit none
  real(dp)::data(80,3),index(80),stats(4)
  real(dp),allocatable::macd(:),signal(:),hist(:),position(:),dd(:),net(:),wealth(:),md(:),cd(:)
  real(dp)::turn(20),gross(20)
  integer,allocatable::idx(:),dir(:)
  type(linear_constraints)::con
  type(backtest_result)::bt
  integer::i,fails
  fails=0
  do i=1,80
    data(i,1)=0.0005_dp+0.01_dp*sin(0.17_dp*real(i,dp))
    data(i,2)=0.0003_dp+0.006_dp*cos(0.11_dp*real(i,dp))
    data(i,3)=0.0007_dp+0.014_dp*sin(0.07_dp*real(i,dp)+0.4_dp)
  end do
  call initialize_constraints(3,con)
  call run_backtest(data,20,5,'minvariance',con,bt,transaction_cost=0.001_dp,smoothing_lambda=0.2_dp)
  call check(size(bt%portfolio_returns)==60,'backtest length',fails)
  call check(all(bt%weights>=-1.0e-12_dp),'backtest nonnegative weights',fails)
  call check(maxval(abs(sum(bt%weights,dim=1)-1.0_dp))<1.0e-8_dp,'backtest budget',fails)
  call check(bt%max_drawdown>=0.0_dp,'backtest drawdown',fails)
  call check(all(rolling_sigma(bt%portfolio_returns,10)>=0.0_dp),'rolling sigma',fails)
  call check(all(rolling_var(bt%portfolio_returns,10,0.1_dp)>=0.0_dp),'rolling VaR',fails)
  call check(all(rolling_cvar(bt%portfolio_returns,10,0.1_dp)>=0.0_dp),'rolling CVaR',fails)
  index(1)=100.0_dp
  do i=2,80;index(i)=index(i-1)*(1.0_dp+data(i,1));end do
  call macd_indicator(index,0.8_dp,0.9_dp,0.85_dp,macd,signal,hist,position)
  call check(size(macd)==80 .and. all(abs(abs(position)-1.0_dp)<1.0e-14_dp),'MACD indicator',fails)
  call drawdown_indicator(data(:,1),0.8_dp,0.9_dp,0.2_dp,dd,macd,signal,hist,position)
  call check(all(dd>=0.0_dp) .and. all(position>=0.0_dp),'drawdown indicator',fails)
  call turning_points(index,0.8_dp,idx,dir)
  call check(size(idx)==size(dir),'turning points',fails)
  call rolling_stability(data,20,md,cd)
  call check(all(md>=0.0_dp) .and. all(cd>=0.0_dp),'rolling stability',fails)
  gross=0.001_dp;turn=0.0_dp;turn(1:20:5)=0.5_dp
  call net_performance(gross,turn,0.001_dp,0.0001_dp,net,wealth)
  call check(wealth(20)>0.0_dp .and. all(net<=gross+1.0e-15_dp),'net performance',fails)
  stats=rebalancing_statistics(data(:,1),position)
  call check(stats(1)>=abs(stats(4)),'rebalancing statistics',fails)
  if(fails>0)error stop 'backtest/monitor tests failed'
  print '(a)','Rolling backtest, transaction-cost, smoothing, and monitoring tests passed.'
contains
  subroutine check(ok,name,fails)
    logical,intent(in)::ok
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    if(.not.ok)then;print '(a)','FAIL: '//trim(name);fails=fails+1;end if
  end subroutine check
end program test_backtest_monitor
