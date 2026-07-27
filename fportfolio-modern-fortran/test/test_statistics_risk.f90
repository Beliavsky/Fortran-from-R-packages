! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
program test_statistics_risk
  use fportfolio_kinds, only: dp
  use fportfolio_statistics
  use fportfolio_risk
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  real(dp)::x(6,3),w(3),alpha,mu0(3),risk
  real(dp),allocatable::mu(:),sigma(:,:),marg(:),comp(:),bud(:),lo(:,:),up(:,:),tailmat(:,:)
  integer::fails
  fails=0
  x=reshape([ &
    -0.02_dp, 0.01_dp, 0.03_dp, 0.00_dp,-0.01_dp, 0.02_dp, &
     0.01_dp, 0.02_dp,-0.01_dp, 0.03_dp, 0.00_dp, 0.01_dp, &
     0.00_dp,-0.01_dp, 0.02_dp, 0.01_dp, 0.04_dp,-0.02_dp],[6,3])
  w=[0.5_dp,0.3_dp,0.2_dp];alpha=0.2_dp
  call sample_estimator(x,mu,sigma)
  mu0=sum(x,dim=1)/6.0_dp
  call check(maxval(abs(mu-mu0))<1.0e-14_dp,'sample mean',fails)
  call check(maxval(abs(sigma-transpose(sigma)))<1.0e-14_dp,'sample covariance symmetry',fails)
  risk=covariance_risk(sigma,w)
  call check(risk>0.0_dp,'covariance risk positive',fails)
  call covariance_risk_contributions(sigma,w,marg,comp,bud)
  call check(abs(sum(comp)-risk)<1.0e-12_dp,'Euler covariance contributions',fails)
  call check(abs(sum(bud)-1.0_dp)<1.0e-12_dp,'risk budgets sum',fails)
  call check(historical_es(x,w,alpha)>=historical_var(x,w,alpha)-1.0e-14_dp,'ES above VaR',fails)
  call check(normal_es(mu,sigma,w,alpha)>=normal_var(mu,sigma,w,alpha),'normal ES above VaR',fails)
  call check(modified_es(x,w,alpha)>=0.0_dp .and. modified_var(x,w,alpha)>=0.0_dp,'modified risk finite',fails)
  call empirical_tail_dependence(x,alpha,lo,up)
  call check(all(lo>=0.0_dp) .and. all(up>=0.0_dp),'tail dependence nonnegative',fails)
  call check(maximum_drawdown(portfolio_returns(x,w))>=0.0_dp,'drawdown nonnegative',fails)
  call historical_es_contributions(x,w,alpha,risk,comp,bud)
  call check(abs(sum(comp)-risk)<1.0e-12_dp,'historical ES contributions',fails)
  call modified_var_contributions(x,w,alpha,risk,marg,comp,bud)
  call check(all(ieee_is_finite(comp)) .and. abs(sum(bud)-1.0_dp)<1.0e-8_dp,'modified VaR budgets',fails)
  call modified_es_contributions(x,w,alpha,risk,marg,comp,bud)
  call check(all(ieee_is_finite(comp)) .and. abs(sum(bud)-1.0_dp)<1.0e-8_dp,'modified ES budgets',fails)
  call normal_margin_tail_dependence(x,'upper',tailmat)
  call check(maxval(abs(tailmat-transpose(tailmat)))<1.0e-14_dp,'normal-margin tail dependence',fails)
  call check(portfolio_max_loss(x,w)>=0.0_dp,'portfolio maximum loss',fails)
  call check(size(geometric_portfolio_returns(x,w))==size(x,1),'geometric portfolio returns',fails)
  call shrinkage_estimator(x,mu,sigma)
  call check(maxval(abs(sigma-transpose(sigma)))<1.0e-14_dp,'shrinkage covariance symmetry',fails)
  call ewma_estimator(x,0.94_dp,mu,sigma)
  call check(maxval(abs(sigma-transpose(sigma)))<1.0e-14_dp,'EWMA covariance symmetry',fails)
  if(fails>0)error stop 'statistics/risk tests failed'
  print '(a)','Statistics, estimators, risk measures, and risk-budget tests passed.'
contains
  subroutine check(ok,name,fails)
    logical,intent(in)::ok
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    if(.not.ok)then;print '(a)','FAIL: '//trim(name);fails=fails+1;end if
  end subroutine check
end program test_statistics_risk
