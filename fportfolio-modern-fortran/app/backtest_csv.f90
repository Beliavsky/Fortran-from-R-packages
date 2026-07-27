! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
program backtest_csv
  use fportfolio_kinds, only: dp
  use fportfolio_types
  use fportfolio_io, only: read_returns_csv
  use fportfolio_optimization, only: initialize_constraints
  use fportfolio_backtest, only: run_backtest
  implicit none
  character(len=512)::filename,strategy,arg
  real(dp),allocatable::data(:,:)
  character(len=64),allocatable::names(:)
  type(linear_constraints)::con
  type(backtest_result)::result
  integer::window,rebalance
  if(command_argument_count()<1)then
    print '(a)','usage: backtest_csv FILE [strategy] [window] [rebalance]'
    stop 1
  end if
  call get_command_argument(1,filename);strategy='minvariance';window=60;rebalance=20
  if(command_argument_count()>=2)call get_command_argument(2,strategy)
  if(command_argument_count()>=3)then;call get_command_argument(3,arg);read(arg,*)window;end if
  if(command_argument_count()>=4)then;call get_command_argument(4,arg);read(arg,*)rebalance;end if
  call read_returns_csv(trim(filename),data,names)
  call initialize_constraints(size(data,2),con)
  call run_backtest(data,window,rebalance,trim(strategy),con,result,transaction_cost=0.001_dp)
  print '(a,a)','strategy: ',trim(strategy)
  print '(a,es16.8)','total return: ',result%total_return
  print '(a,es16.8)','annualized return: ',result%annualized_return
  print '(a,es16.8)','annualized volatility: ',result%annualized_volatility
  print '(a,es16.8)','Sharpe ratio: ',result%sharpe
  print '(a,es16.8)','maximum drawdown: ',result%max_drawdown
  print '(a,es16.8)','total turnover: ',sum(result%turnover)
end program backtest_csv
