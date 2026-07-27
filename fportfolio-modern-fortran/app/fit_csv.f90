! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
program fit_csv
  use fportfolio_kinds, only: dp
  use fportfolio_types
  use fportfolio_io, only: read_returns_csv, write_weights_csv
  use fportfolio_statistics, only: sample_estimator, shrinkage_estimator, ewma_estimator
  use fportfolio_optimization
  use fportfolio_risk, only: portfolio_risk_report
  implicit none
  character(len=512)::filename,method,arg,outfile
  real(dp),allocatable::data(:,:),mu(:),sigma(:,:)
  character(len=64),allocatable::names(:)
  type(linear_constraints)::con
  type(optimizer_result)::res
  type(risk_report)::report
  real(dp)::alpha,rf,target
  integer::narg,i
  narg=command_argument_count()
  if(narg<1)then
    print '(a)','usage: fit_csv FILE [minvariance|tangency|riskparity|maxdiv|mad|cvar|efficient|maxreturn] [value]'
    stop 1
  end if
  call get_command_argument(1,filename);method='minvariance';if(narg>=2)call get_command_argument(2,method)
  alpha=0.05_dp;rf=0.0_dp;target=0.0_dp
  if(narg>=3)then;call get_command_argument(3,arg);read(arg,*)target;alpha=target;rf=target;end if
  call read_returns_csv(trim(filename),data,names)
  call sample_estimator(data,mu,sigma)
  call initialize_constraints(size(mu),con)
  select case(trim(adjustl(method)))
  case('minvariance')
    call minvariance_portfolio(mu,sigma,con,res)
  case('tangency')
    call tangency_portfolio(mu,sigma,rf,con,res)
  case('riskparity')
    call risk_parity_portfolio(mu,sigma,con,res)
  case('maxdiv')
    call maximum_diversification_portfolio(mu,sigma,con,res)
  case('mad')
    call minimum_mad_portfolio(data,mu,con,res)
  case('cvar')
    if(narg<3)alpha=0.05_dp
    call minimum_cvar_portfolio(data,mu,alpha,con,res)
  case('efficient')
    if(narg<3)target=sum(mu)/real(size(mu),dp)
    call efficient_portfolio(mu,sigma,target,con,res)
  case('maxreturn')
    call maxreturn_portfolio(mu,sigma,con,res)
  case default
    error stop 'unknown portfolio method'
  end select
  if(.not.res%converged)print '(a)','warning: optimizer did not meet its convergence criterion'
  call portfolio_risk_report(data,res%weights,0.05_dp,report)
  print '(a,a)','method: ',trim(method)
  print '(a,i0)','observations: ',size(data,1)
  print '(a,i0)','assets: ',size(data,2)
  print '(a,es16.8)','expected return: ',report%expected_return
  print '(a,es16.8)','volatility: ',report%volatility
  print '(a,es16.8)','historical VaR 5%: ',report%var
  print '(a,es16.8)','historical ES 5%: ',report%es
  print '(a)','weights:'
  do i=1,size(names);print '(2x,a,1x,f12.8)',trim(names(i)),res%weights(i);end do
  outfile='weights.csv';call write_weights_csv(trim(outfile),names,res%weights)
  print '(a)','wrote weights.csv'
end program fit_csv
