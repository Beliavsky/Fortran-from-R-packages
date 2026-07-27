! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
program test_optimization
  use fportfolio_kinds, only: dp
  use fportfolio_types
  use fportfolio_optimization
  use fportfolio_risk, only: covariance_risk_contributions, historical_es
  implicit none
  real(dp)::mu(3),sigma(3,3),expected(3),data(12,3),target
  real(dp),allocatable::m(:),c(:),b(:),feasible_w(:)
  type(linear_constraints)::con,ct
  type(optimizer_result)::res,res2
  type(frontier_result)::front
  integer::fails,status
  fails=0;mu=[0.05_dp,0.08_dp,0.12_dp];sigma=0.0_dp
  sigma(1,1)=0.04_dp;sigma(2,2)=0.09_dp;sigma(3,3)=0.16_dp
  call initialize_constraints(3,con)
  call minvariance_portfolio(mu,sigma,con,res)
  expected=[25.0_dp,1.0_dp/0.09_dp,6.25_dp];expected=expected/sum(expected)
  call check(res%converged,'minimum variance converged',fails)
  call check(maxval(abs(res%weights-expected))<2.0e-5_dp,'minimum variance analytic weights',fails)
  target=0.08_dp;call efficient_portfolio(mu,sigma,target,con,res)
  call check(res%converged .and. abs(dot_product(mu,res%weights)-target)<2.0e-6_dp,'efficient target',fails)
  call maxreturn_portfolio(mu,sigma,con,res)
  call check(res%weights(3)>0.999_dp,'maximum return corner',fails)
  call tangency_portfolio(mu,sigma,0.01_dp,con,res)
  call check(res%converged .and. res%sharpe>0.0_dp,'tangency portfolio',fails)
  call risk_parity_portfolio(mu,sigma,con,res)
  call covariance_risk_contributions(sigma,res%weights,m,c,b)
  call check(maxval(abs(b-1.0_dp/3.0_dp))<2.0e-4_dp,'risk parity equal budgets',fails)
  call maximum_diversification_portfolio(mu,sigma,con,res2)
  call check(maxval(abs(res2%weights-[6.0_dp/13.0_dp,4.0_dp/13.0_dp,3.0_dp/13.0_dp]))<2.0e-4_dp, &
    'maximum diversification diagonal solution',fails)
  call portfolio_frontier(mu,sigma,con,8,front)
  call check(all(front%feasible),'frontier feasible',fails)
  call check(all(front%risk>=0.0_dp),'frontier risks',fails)
  data=reshape([ &
    -0.04_dp,-0.02_dp,0.01_dp,0.03_dp,0.02_dp,-0.01_dp,0.01_dp,0.02_dp,0.00_dp,0.01_dp,-0.02_dp,0.03_dp, &
    -0.01_dp,0.00_dp,0.02_dp,0.01_dp,0.03_dp,0.01_dp,-0.01_dp,0.01_dp,0.02_dp,0.00_dp,0.01_dp,0.02_dp, &
    -0.08_dp,0.04_dp,0.05_dp,-0.03_dp,0.06_dp,-0.02_dp,0.04_dp,-0.01_dp,0.03_dp,0.02_dp,-0.04_dp,0.05_dp],[12,3])
  call minimum_mad_portfolio(data,sum(data,dim=1)/12.0_dp,con,res,max_iter=5000)
  call check(res%converged .and. abs(sum(res%weights)-1.0_dp)<1.0e-8_dp,'minimum MAD',fails)
  call minimum_cvar_portfolio(data,sum(data,dim=1)/12.0_dp,0.25_dp,con,res,max_iter=8000)
  call check(abs(sum(res%weights)-1.0_dp)<1.0e-8_dp .and. all(res%weights>=-1.0e-12_dp), &
    'minimum CVaR feasibility',fails)
  call check(res%risk<=historical_es(data,[1.0_dp/3.0_dp,1.0_dp/3.0_dp,1.0_dp/3.0_dp],0.25_dp)+1.0e-8_dp, &
    'minimum CVaR improves equal weights',fails)
  call cardinality_minvariance_portfolio(mu,sigma,con,2,0.1_dp,res)
  call check(res%converged .and. count(res%weights>1.0e-8_dp)<=2,'cardinality constraint',fails)
  ct=con;allocate(ct%a_ineq(1,3),ct%b_ineq(1));ct%a_ineq(1,:)=[1.0_dp,1.0_dp,0.0_dp];ct%b_ineq=0.7_dp
  call minvariance_portfolio(mu,sigma,ct,res)
  call check(sum(res%weights(:2))<=0.70001_dp,'linear inequality constraint',fails)
  call feasible_portfolio(con,feasible_w,status)
  call check(status==0 .and. abs(sum(feasible_w)-1.0_dp)<1.0e-12_dp,'feasible portfolio',fails)
  if(fails>0)error stop 'optimization tests failed'
  print '(a)','Mean-variance, frontier, ratio, MAD, CVaR, risk-parity, and cardinality tests passed.'
contains
  subroutine check(ok,name,fails)
    logical,intent(in)::ok
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    if(.not.ok)then;print '(a)','FAIL: '//trim(name);fails=fails+1;end if
  end subroutine check
end program test_optimization
