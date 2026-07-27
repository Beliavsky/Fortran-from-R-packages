! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_risk
  use fportfolio_kinds, only: dp
  use fportfolio_probability, only: normal_pdf, normal_cdf, normal_quantile
  use fportfolio_statistics, only: quantile_value, skewness_value, kurtosis_value, column_means, sample_covariance
  use fportfolio_linalg, only: quadratic_form
  use fportfolio_types, only: risk_report
  implicit none
  private
  public :: portfolio_returns, covariance_risk, historical_var, historical_es, &
            normal_var, normal_es, modified_var, modified_es, maximum_drawdown, &
            drawdown_series, drawdown_at_risk, conditional_drawdown_at_risk, &
            covariance_risk_contributions, normal_var_contributions, normal_es_contributions, &
            diversification_ratio, empirical_tail_dependence, portfolio_risk_report, &
            sample_var, sample_cvar, sample_covar, historical_var_contributions, &
            historical_es_contributions, modified_var_contributions, modified_es_contributions, &
            portfolio_max_loss, geometric_portfolio_returns, cfg_tail_dependence, normal_margin_tail_dependence
contains
  pure function portfolio_returns(data,weights) result(r)
    real(dp), intent(in) :: data(:,:),weights(:)
    real(dp) :: r(size(data,1))
    r=matmul(data,weights)
  end function portfolio_returns

  pure real(dp) function covariance_risk(sigma,weights) result(risk)
    real(dp), intent(in) :: sigma(:,:),weights(:)
    risk=sqrt(max(quadratic_form(weights,sigma),0.0_dp))
  end function covariance_risk

  real(dp) function historical_var(data,weights,alpha) result(v)
    real(dp), intent(in) :: data(:,:),weights(:),alpha
    real(dp), allocatable :: r(:)
    allocate(r(size(data,1)))
    r=portfolio_returns(data,weights)
    v=max(0.0_dp,-quantile_value(r,alpha))
  end function historical_var

  real(dp) function historical_es(data,weights,alpha) result(es)
    real(dp), intent(in) :: data(:,:),weights(:),alpha
    real(dp), allocatable :: r(:)
    real(dp) :: q
    integer :: n
    allocate(r(size(data,1)))
    r=portfolio_returns(data,weights)
    q=quantile_value(r,alpha)
    n=count(r<=q)
    if(n>0)then;es=max(0.0_dp,-sum(r,mask=r<=q)/real(n,dp));else;es=max(0.0_dp,-q);end if
  end function historical_es

  pure real(dp) function normal_var(mu,sigma,weights,alpha) result(v)
    real(dp), intent(in) :: mu(:),sigma(:,:),weights(:),alpha
    real(dp)::m,s,z
    m=dot_product(mu,weights);s=covariance_risk(sigma,weights);z=normal_quantile(alpha)
    v=max(0.0_dp,-m-z*s)
  end function normal_var

  pure real(dp) function normal_es(mu,sigma,weights,alpha) result(es)
    real(dp), intent(in) :: mu(:),sigma(:,:),weights(:),alpha
    real(dp)::m,s,z
    m=dot_product(mu,weights);s=covariance_risk(sigma,weights);z=normal_quantile(alpha)
    es=max(0.0_dp,-m+s*normal_pdf(z)/max(alpha,tiny(1.0_dp)))
  end function normal_es

  real(dp) function modified_var(data,weights,alpha) result(v)
    real(dp), intent(in)::data(:,:),weights(:),alpha
    real(dp),allocatable::r(:)
    real(dp)::m,s,z,sk,ku,zcf
    allocate(r(size(data,1)))
    r=portfolio_returns(data,weights);m=sum(r)/real(size(r),dp)
    s=sqrt(sum((r-m)**2)/real(max(1,size(r)-1),dp));z=normal_quantile(alpha)
    sk=skewness_value(r);ku=kurtosis_value(r)-3.0_dp
    zcf=z+(z*z-1.0_dp)*sk/6.0_dp+(z**3-3.0_dp*z)*ku/24.0_dp-(2.0_dp*z**3-5.0_dp*z)*sk*sk/36.0_dp
    v=max(0.0_dp,-(m+s*zcf))
  end function modified_var

  real(dp) function modified_es(data,weights,alpha) result(es)
    real(dp), intent(in)::data(:,:),weights(:),alpha
    real(dp),allocatable::r(:)
    real(dp)::m,s,sk,ku,p,z,zcf,acc
    integer::i,nq
    allocate(r(size(data,1)))
    r=portfolio_returns(data,weights);m=sum(r)/real(size(r),dp)
    s=sqrt(sum((r-m)**2)/real(max(1,size(r)-1),dp));sk=skewness_value(r);ku=kurtosis_value(r)-3.0_dp
    nq=400;acc=0.0_dp
    do i=1,nq
      p=alpha*(real(i,dp)-0.5_dp)/real(nq,dp);z=normal_quantile(p)
      zcf=z+(z*z-1.0_dp)*sk/6.0_dp+(z**3-3.0_dp*z)*ku/24.0_dp-(2.0_dp*z**3-5.0_dp*z)*sk*sk/36.0_dp
      acc=acc+m+s*zcf
    end do
    es=max(0.0_dp,-acc/real(nq,dp))
  end function modified_es

  function drawdown_series(returns) result(dd)
    real(dp), intent(in)::returns(:)
    real(dp)::dd(size(returns)),wealth,peak
    integer::i
    wealth=1.0_dp;peak=1.0_dp
    do i=1,size(returns)
      wealth=wealth*(1.0_dp+returns(i));peak=max(peak,wealth);dd(i)=1.0_dp-wealth/peak
    end do
  end function drawdown_series

  real(dp) function maximum_drawdown(returns) result(mdd)
    real(dp), intent(in)::returns(:)
    real(dp)::dd(size(returns));dd=drawdown_series(returns);mdd=maxval(dd)
  end function maximum_drawdown

  real(dp) function drawdown_at_risk(returns,alpha) result(dar)
    real(dp), intent(in)::returns(:),alpha
    real(dp)::dd(size(returns));dd=drawdown_series(returns);dar=quantile_value(dd,1.0_dp-alpha)
  end function drawdown_at_risk

  real(dp) function conditional_drawdown_at_risk(returns,alpha) result(cdar)
    real(dp), intent(in)::returns(:),alpha
    real(dp)::dd(size(returns)),q
    integer::n
    dd=drawdown_series(returns);q=quantile_value(dd,1.0_dp-alpha);n=count(dd>=q)
    if(n>0)then;cdar=sum(dd,mask=dd>=q)/real(n,dp);else;cdar=q;end if
  end function conditional_drawdown_at_risk

  subroutine covariance_risk_contributions(sigma,weights,marginal,component,budget)
    real(dp), intent(in)::sigma(:,:),weights(:)
    real(dp),allocatable,intent(out)::marginal(:),component(:),budget(:)
    real(dp)::risk
    integer::n
    n=size(weights);allocate(marginal(n),component(n),budget(n));risk=covariance_risk(sigma,weights)
    if(risk>0.0_dp)then
      marginal=matmul(sigma,weights)/risk;component=weights*marginal
      budget=component/risk
    else
      marginal=0.0_dp;component=0.0_dp;budget=0.0_dp
    end if
  end subroutine covariance_risk_contributions

  subroutine normal_var_contributions(mu,sigma,weights,alpha,marginal,component,budget)
    real(dp),intent(in)::mu(:),sigma(:,:),weights(:),alpha
    real(dp),allocatable,intent(out)::marginal(:),component(:),budget(:)
    real(dp)::s,z,total
    integer::n
    n=size(weights);allocate(marginal(n),component(n),budget(n));s=covariance_risk(sigma,weights);z=normal_quantile(alpha)
    if(s>0.0_dp)then;marginal=-mu-z*matmul(sigma,weights)/s;else;marginal=-mu;end if
    component=weights*marginal;total=sum(component)
    if(abs(total)>tiny(1.0_dp))then;budget=component/total;else;budget=0.0_dp;end if
  end subroutine normal_var_contributions

  subroutine normal_es_contributions(mu,sigma,weights,alpha,marginal,component,budget)
    real(dp),intent(in)::mu(:),sigma(:,:),weights(:),alpha
    real(dp),allocatable,intent(out)::marginal(:),component(:),budget(:)
    real(dp)::s,z,k,total
    integer::n
    n=size(weights);allocate(marginal(n),component(n),budget(n));s=covariance_risk(sigma,weights);z=normal_quantile(alpha)
    k=normal_pdf(z)/max(alpha,tiny(1.0_dp))
    if(s>0.0_dp)then;marginal=-mu+k*matmul(sigma,weights)/s;else;marginal=-mu;end if
    component=weights*marginal;total=sum(component)
    if(abs(total)>tiny(1.0_dp))then;budget=component/total;else;budget=0.0_dp;end if
  end subroutine normal_es_contributions

  pure real(dp) function diversification_ratio(sigma,weights) result(dr)
    real(dp),intent(in)::sigma(:,:),weights(:)
    real(dp)::s(size(weights)),den
    integer::i
    do i=1,size(weights);s(i)=sqrt(max(sigma(i,i),0.0_dp));end do
    den=covariance_risk(sigma,weights)
    if(den>0.0_dp)then;dr=dot_product(weights,s)/den;else;dr=0.0_dp;end if
  end function diversification_ratio

  subroutine empirical_tail_dependence(data,alpha,lower,upper)
    real(dp),intent(in)::data(:,:),alpha
    real(dp),allocatable,intent(out)::lower(:,:),upper(:,:)
    real(dp),allocatable::lo(:),hi(:)
    integer::p,n,i,j
    n=size(data,1);p=size(data,2);allocate(lower(p,p),upper(p,p),lo(p),hi(p))
    do i=1,p;lo(i)=quantile_value(data(:,i),alpha);hi(i)=quantile_value(data(:,i),1.0_dp-alpha);end do
    do i=1,p;do j=1,p
      lower(i,j)=real(count(data(:,i)<=lo(i) .and. data(:,j)<=lo(j)),dp)/max(real(n,dp)*alpha,tiny(1.0_dp))
      upper(i,j)=real(count(data(:,i)>=hi(i) .and. data(:,j)>=hi(j)),dp)/max(real(n,dp)*alpha,tiny(1.0_dp))
    end do;end do
  end subroutine empirical_tail_dependence

  subroutine portfolio_risk_report(data,weights,alpha,report)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    type(risk_report),intent(out)::report
    real(dp)::mu(size(data,2)),sigma(size(data,2),size(data,2)),r(size(data,1))
    mu=column_means(data);sigma=sample_covariance(data);r=portfolio_returns(data,weights)
    report%expected_return=dot_product(mu,weights);report%volatility=covariance_risk(sigma,weights)
    report%var=historical_var(data,weights,alpha);report%es=historical_es(data,weights,alpha)
    report%max_drawdown=maximum_drawdown(r);report%dar=drawdown_at_risk(r,alpha);report%cdar=conditional_drawdown_at_risk(r,alpha)
    call covariance_risk_contributions(sigma,weights,report%marginal,report%component,report%budget)
  end subroutine portfolio_risk_report

  subroutine sample_covar(data,weights,risk,marginal,component,budget)
    real(dp),intent(in)::data(:,:),weights(:)
    real(dp),intent(out)::risk
    real(dp),allocatable,intent(out)::marginal(:),component(:),budget(:)
    real(dp)::sigma(size(data,2),size(data,2))
    sigma=sample_covariance(data);risk=covariance_risk(sigma,weights)
    call covariance_risk_contributions(sigma,weights,marginal,component,budget)
  end subroutine sample_covar

  subroutine sample_var(data,weights,alpha,risk)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    real(dp),intent(out)::risk
    risk=historical_var(data,weights,alpha)
  end subroutine sample_var

  subroutine sample_cvar(data,weights,alpha,risk)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    real(dp),intent(out)::risk
    risk=historical_es(data,weights,alpha)
  end subroutine sample_cvar

  subroutine historical_var_contributions(data,weights,alpha,risk,component,budget)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    real(dp),intent(out)::risk
    real(dp),allocatable,intent(out)::component(:),budget(:)
    real(dp),allocatable::r(:)
    real(dp)::q,total
    integer::i,idx,n
    n=size(data,1);allocate(r(n),component(size(weights)),budget(size(weights)))
    r=portfolio_returns(data,weights);q=quantile_value(r,alpha);idx=1
    do i=2,n
      if(abs(r(i)-q)<abs(r(idx)-q))idx=i
    end do
    component=-weights*data(idx,:);risk=sum(component);risk=max(risk,0.0_dp);total=sum(component)
    if(abs(total)>tiny(1.0_dp))then;budget=component/total;else;budget=0.0_dp;end if
  end subroutine historical_var_contributions

  subroutine historical_es_contributions(data,weights,alpha,risk,component,budget)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    real(dp),intent(out)::risk
    real(dp),allocatable,intent(out)::component(:),budget(:)
    real(dp),allocatable::r(:),tail_mean(:)
    real(dp)::q,total
    integer::n,nt,i
    n=size(data,1);allocate(r(n),component(size(weights)),budget(size(weights)),tail_mean(size(weights)))
    r=portfolio_returns(data,weights);q=quantile_value(r,alpha);nt=count(r<=q);tail_mean=0.0_dp
    if(nt>0)then
      do i=1,n
        if(r(i)<=q)tail_mean=tail_mean+data(i,:)
      end do
      tail_mean=tail_mean/real(nt,dp)
    end if
    component=-weights*tail_mean;risk=sum(component);risk=max(risk,0.0_dp);total=sum(component)
    if(abs(total)>tiny(1.0_dp))then;budget=component/total;else;budget=0.0_dp;end if
  end subroutine historical_es_contributions

  subroutine modified_var_contributions(data,weights,alpha,risk,marginal,component,budget)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    real(dp),intent(out)::risk
    real(dp),allocatable,intent(out)::marginal(:),component(:),budget(:)
    real(dp),allocatable::xc(:,:),r(:),dm2(:),dm3(:),dm4(:),dsk(:),dku(:)
    real(dp)::mu(size(weights)),m,m2,m3,m4,s,sk,ku,z,h,total
    integer::n,p,j
    n=size(data,1);p=size(data,2);allocate(xc(n,p),r(n),dm2(p),dm3(p),dm4(p),dsk(p),dku(p))
    allocate(marginal(p),component(p),budget(p));mu=column_means(data);xc=data-spread(mu,1,n);r=matmul(xc,weights)
    m=dot_product(mu,weights);m2=sum(r*r)/real(n,dp);m3=sum(r**3)/real(n,dp);m4=sum(r**4)/real(n,dp)
    s=sqrt(max(m2,tiny(1.0_dp)));sk=m3/s**3;ku=m4/m2**2-3.0_dp;z=normal_quantile(alpha)
    h=z+(z*z-1.0_dp)*sk/6.0_dp+(z**3-3.0_dp*z)*ku/24.0_dp-(2.0_dp*z**3-5.0_dp*z)*sk*sk/36.0_dp
    do j=1,p
      dm2(j)=2.0_dp*sum(r*xc(:,j))/real(n,dp)
      dm3(j)=3.0_dp*sum(r*r*xc(:,j))/real(n,dp)
      dm4(j)=4.0_dp*sum(r**3*xc(:,j))/real(n,dp)
      dsk(j)=dm3(j)/s**3-1.5_dp*m3*dm2(j)/s**5
      dku(j)=dm4(j)/m2**2-2.0_dp*m4*dm2(j)/m2**3
    end do
    marginal=-mu-h*dm2/(2.0_dp*s)-s*((z*z-1.0_dp)*dsk/6.0_dp+ &
      (z**3-3.0_dp*z)*dku/24.0_dp-(2.0_dp*z**3-5.0_dp*z)*sk*dsk/18.0_dp)
    component=weights*marginal;risk=-(m+s*h);risk=max(risk,0.0_dp);total=sum(component)
    if(abs(total)>tiny(1.0_dp))then;budget=component/total;else;budget=0.0_dp;end if
  end subroutine modified_var_contributions

  subroutine modified_es_contributions(data,weights,alpha,risk,marginal,component,budget)
    real(dp),intent(in)::data(:,:),weights(:),alpha
    real(dp),intent(out)::risk
    real(dp),allocatable,intent(out)::marginal(:),component(:),budget(:)
    real(dp),allocatable::xc(:,:),r(:),dm2(:),dm3(:),dm4(:),dsk(:),dku(:),davg(:)
    real(dp)::mu(size(weights)),m,m2,m3,m4,s,sk,ku,pv,z,zcf,avg,total
    integer::n,p,j,i,nq
    n=size(data,1);p=size(data,2);nq=500
    allocate(xc(n,p),r(n),dm2(p),dm3(p),dm4(p),dsk(p),dku(p),davg(p))
    allocate(marginal(p),component(p),budget(p));mu=column_means(data);xc=data-spread(mu,1,n);r=matmul(xc,weights)
    m=dot_product(mu,weights);m2=sum(r*r)/real(n,dp);m3=sum(r**3)/real(n,dp);m4=sum(r**4)/real(n,dp)
    s=sqrt(max(m2,tiny(1.0_dp)));sk=m3/s**3;ku=m4/m2**2-3.0_dp
    do j=1,p
      dm2(j)=2.0_dp*sum(r*xc(:,j))/real(n,dp);dm3(j)=3.0_dp*sum(r*r*xc(:,j))/real(n,dp)
      dm4(j)=4.0_dp*sum(r**3*xc(:,j))/real(n,dp)
      dsk(j)=dm3(j)/s**3-1.5_dp*m3*dm2(j)/s**5;dku(j)=dm4(j)/m2**2-2.0_dp*m4*dm2(j)/m2**3
    end do
    avg=0.0_dp;davg=0.0_dp
    do i=1,nq
      pv=alpha*(real(i,dp)-0.5_dp)/real(nq,dp);z=normal_quantile(pv)
      zcf=z+(z*z-1.0_dp)*sk/6.0_dp+(z**3-3.0_dp*z)*ku/24.0_dp-(2.0_dp*z**3-5.0_dp*z)*sk*sk/36.0_dp
      avg=avg+zcf
      davg=davg+(z*z-1.0_dp)*dsk/6.0_dp+(z**3-3.0_dp*z)*dku/24.0_dp- &
        (2.0_dp*z**3-5.0_dp*z)*sk*dsk/18.0_dp
    end do
    avg=avg/real(nq,dp);davg=davg/real(nq,dp)
    marginal=-mu-avg*dm2/(2.0_dp*s)-s*davg;component=weights*marginal;risk=-(m+s*avg);risk=max(risk,0.0_dp)
    total=sum(component);if(abs(total)>tiny(1.0_dp))then;budget=component/total;else;budget=0.0_dp;end if
  end subroutine modified_es_contributions

  pure real(dp) function portfolio_max_loss(data,weights) result(loss)
    real(dp),intent(in)::data(:,:),weights(:)
    loss=max(0.0_dp,-minval(matmul(data,weights)))
  end function portfolio_max_loss

  function geometric_portfolio_returns(data,weights) result(r)
    real(dp),intent(in)::data(:,:),weights(:)
    real(dp)::r(size(data,1)),asset_wealth(size(weights)),old_portfolio,new_portfolio
    integer::t
    asset_wealth=1.0_dp;old_portfolio=dot_product(weights,asset_wealth)
    do t=1,size(data,1)
      asset_wealth=asset_wealth*(1.0_dp+data(t,:));new_portfolio=dot_product(weights,asset_wealth)
      r(t)=new_portfolio/old_portfolio-1.0_dp;old_portfolio=new_portfolio
    end do
  end function geometric_portfolio_returns

  real(dp) function cfg_tail_dependence(u,v,tail) result(lambda)
    real(dp),intent(in)::u(:),v(:)
    character(len=*),intent(in),optional::tail
    real(dp)::x,y,acc
    integer::i,n
    n=size(u);acc=0.0_dp
    do i=1,n
      x=max(min(u(i),1.0_dp-1.0e-12_dp),1.0e-12_dp);y=max(min(v(i),1.0_dp-1.0e-12_dp),1.0e-12_dp)
      if(present(tail))then;if(trim(tail)=="lower")then;x=1.0_dp-x;y=1.0_dp-y;end if;end if
      acc=acc+log(sqrt(log(1.0_dp/x)*log(1.0_dp/y))/log(1.0_dp/max(x,y)**2))
    end do
    lambda=2.0_dp-2.0_dp*exp(acc/real(n,dp));lambda=max(0.0_dp,min(1.0_dp,lambda))
  end function cfg_tail_dependence

  subroutine normal_margin_tail_dependence(data,tail,matrix)
    real(dp),intent(in)::data(:,:)
    character(len=*),intent(in)::tail
    real(dp),allocatable,intent(out)::matrix(:,:)
    real(dp),allocatable::u(:,:)
    real(dp)::mu(size(data,2)),sd(size(data,2))
    integer::i,j,n,p
    n=size(data,1);p=size(data,2);allocate(u(n,p),matrix(p,p));mu=column_means(data)
    do j=1,p
      sd(j)=sqrt(sum((data(:,j)-mu(j))**2)/real(max(1,n-1),dp))
      do i=1,n;u(i,j)=normal_cdf((data(i,j)-mu(j))/max(sd(j),tiny(1.0_dp)));end do
    end do
    matrix=0.0_dp
    do i=1,p;matrix(i,i)=1.0_dp;do j=i+1,p;matrix(i,j)=cfg_tail_dependence(u(:,i),u(:,j),tail);matrix(j,i)=matrix(i,j);end do;end do
  end subroutine normal_margin_tail_dependence

end module fportfolio_risk
