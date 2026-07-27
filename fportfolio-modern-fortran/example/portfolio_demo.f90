! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
program portfolio_demo
  use fportfolio_kinds, only: dp
  use fportfolio_types
  use fportfolio_statistics, only: sample_estimator
  use fportfolio_optimization
  implicit none
  real(dp)::data(120,4)
  real(dp),allocatable::mu(:),sigma(:,:)
  type(linear_constraints)::con
  type(optimizer_result)::mv,tan,rp,cvar
  type(frontier_result)::front
  integer::i
  do i=1,120
    data(i,1)=0.0004_dp+0.009_dp*sin(0.11_dp*real(i,dp))
    data(i,2)=0.0003_dp+0.006_dp*cos(0.07_dp*real(i,dp))
    data(i,3)=0.0006_dp+0.013_dp*sin(0.05_dp*real(i,dp)+0.8_dp)
    data(i,4)=0.0002_dp+0.004_dp*cos(0.13_dp*real(i,dp)+0.3_dp)
  end do
  call sample_estimator(data,mu,sigma);call initialize_constraints(4,con)
  call minvariance_portfolio(mu,sigma,con,mv)
  call tangency_portfolio(mu,sigma,0.0_dp,con,tan)
  call risk_parity_portfolio(mu,sigma,con,rp)
  call minimum_cvar_portfolio(data,mu,0.05_dp,con,cvar,max_iter=10000)
  call portfolio_frontier(mu,sigma,con,12,front)
  print '(a,4f10.6)','minimum variance weights: ',mv%weights
  print '(a,4f10.6)','tangency weights:        ',tan%weights
  print '(a,4f10.6)','risk parity weights:     ',rp%weights
  print '(a,4f10.6)','minimum CVaR weights:    ',cvar%weights
  print '(a,i0)','frontier points: ',size(front%risk)
end program portfolio_demo
