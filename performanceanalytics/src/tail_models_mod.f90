! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module tail_models_mod
  use kinds_mod, only: dp
  use probability_mod, only: normal_quantile, normal_cdf
  use statistics_mod, only: mean_value, sd_value, quantile_type7, sort_real
  use comoments_mod, only: covariance_matrix
  use risk_mod, only: historical_var, historical_es, component_var
  use performance_ratios_mod, only: sharpe_ratio
  use rng_mod, only: rng_state, rng_seed, rng_integer, multivariate_normal_draws
  implicit none
  private
  public :: gpd_fit_result, gpd_fit, gpd_log_likelihood, gpd_var_value, gpd_es_value
  public :: monte_carlo_asset_risk, monte_carlo_portfolio_risk
  public :: bootstrap_risk_standard_errors
  public :: lognormal_var, lognormal_es, kernel_portfolio_risk

  type :: gpd_fit_result
    real(dp) :: threshold = 0.0_dp
    real(dp) :: scale = 0.0_dp
    real(dp) :: shape = 0.0_dp
    real(dp) :: exceedance_rate = 0.0_dp
    integer :: observations = 0
    integer :: exceedances = 0
    real(dp) :: log_likelihood = -huge(1.0_dp)
    real(dp) :: var_value = 0.0_dp
    real(dp) :: es_value = 0.0_dp
    real(dp) :: var_standard_error = 0.0_dp
    real(dp) :: es_standard_error = 0.0_dp
    real(dp) :: var_lower = 0.0_dp
    real(dp) :: var_upper = 0.0_dp
    real(dp) :: es_lower = 0.0_dp
    real(dp) :: es_upper = 0.0_dp
    logical :: converged = .false.
    integer :: iterations = 0
  end type gpd_fit_result
contains
  pure real(dp) function gpd_log_likelihood(excess, scale, shape) result(ll)
    real(dp), intent(in) :: excess(:), scale, shape
    real(dp), allocatable :: z(:)
    if (scale <= 0.0_dp .or. size(excess) == 0) then
      ll = -huge(1.0_dp)
      return
    end if
    if (abs(shape) < 1.0e-8_dp) then
      ll = -real(size(excess),dp)*log(scale)-sum(excess)/scale
    else
      allocate(z(size(excess)));z=1.0_dp+shape*excess/scale
      if(any(z<=0.0_dp))then;ll=-huge(1.0_dp);return;end if
      ll=-real(size(excess),dp)*log(scale)-(1.0_dp+1.0_dp/shape)*sum(log(z))
    end if
  end function gpd_log_likelihood

  pure real(dp) function gpd_nll_trans(theta, excess) result(v)
    real(dp), intent(in) :: theta(2), excess(:)
    real(dp) :: scale, shape, ll
    scale=exp(theta(1));shape=0.95_dp*tanh(theta(2))
    ll=gpd_log_likelihood(excess,scale,shape)
    if(ll<=-0.5_dp*huge(1.0_dp))then;v=huge(1.0_dp)/100.0_dp;else;v=-ll;end if
  end function gpd_nll_trans

  subroutine nelder_mead_gpd(excess, start, theta, objective, iterations, converged)
    real(dp),intent(in)::excess(:),start(2)
    real(dp),intent(out)::theta(2),objective
    integer,intent(out)::iterations
    logical,intent(out)::converged
    real(dp)::x(2,3),f(3),cent(2),xr(2),xe(2),xc(2),fr,fe,fc,tmpv,tmpx(2),diameter
    integer::i,j,imax,maxit
    real(dp),parameter::alpha=1.0_dp,gamma=2.0_dp,rho=0.5_dp,sigma=0.5_dp
    maxit=800;x(:,1)=start;x(:,2)=start+[0.15_dp,0.0_dp];x(:,3)=start+[0.0_dp,0.15_dp]
    do i=1,3;f(i)=gpd_nll_trans(x(:,i),excess);end do
    converged=.false.
    do iterations=1,maxit
      do i=1,2
        imax=i
        do j=i+1,3;if(f(j)<f(imax))imax=j;end do
        if(imax/=i)then;tmpv=f(i);f(i)=f(imax);f(imax)=tmpv;tmpx=x(:,i);x(:,i)=x(:,imax);x(:,imax)=tmpx;end if
      end do
      diameter=maxval(abs(x-spread(x(:,1),2,3)))
      if(diameter<1.0e-8_dp .and. maxval(abs(f-f(1)))<1.0e-9_dp)then;converged=.true.;exit;end if
      cent=0.5_dp*(x(:,1)+x(:,2));xr=cent+alpha*(cent-x(:,3));fr=gpd_nll_trans(xr,excess)
      if(fr<f(1))then
        xe=cent+gamma*(xr-cent);fe=gpd_nll_trans(xe,excess)
        if(fe<fr)then;x(:,3)=xe;f(3)=fe;else;x(:,3)=xr;f(3)=fr;end if
      else if(fr<f(2))then
        x(:,3)=xr;f(3)=fr
      else
        if(fr<f(3))then;xc=cent+rho*(xr-cent);else;xc=cent-rho*(cent-x(:,3));end if
        fc=gpd_nll_trans(xc,excess)
        if(fc<min(fr,f(3)))then;x(:,3)=xc;f(3)=fc
        else
          x(:,2)=x(:,1)+sigma*(x(:,2)-x(:,1));x(:,3)=x(:,1)+sigma*(x(:,3)-x(:,1))
          f(2)=gpd_nll_trans(x(:,2),excess);f(3)=gpd_nll_trans(x(:,3),excess)
        end if
      end if
    end do
    imax=minloc(f,dim=1);theta=x(:,imax);objective=f(imax)
  end subroutine nelder_mead_gpd

  pure real(dp) function gpd_var_value(threshold,scale,shape,tail_probability,exceedance_rate) result(v)
    real(dp),intent(in)::threshold,scale,shape,tail_probability,exceedance_rate
    real(dp)::ratio
    ratio=max(tail_probability/max(exceedance_rate,tiny(1.0_dp)),tiny(1.0_dp))
    if(abs(shape)<1.0e-8_dp)then
      v=threshold-scale*log(ratio)
    else
      v=threshold+scale/shape*(ratio**(-shape)-1.0_dp)
    end if
  end function gpd_var_value

  pure real(dp) function gpd_es_value(threshold,scale,shape,var_value) result(v)
    real(dp),intent(in)::threshold,scale,shape,var_value
    if(shape>=1.0_dp)then;v=huge(1.0_dp);else;v=(var_value+scale-shape*threshold)/(1.0_dp-shape);end if
  end function gpd_es_value

  subroutine gpd_hessian(excess,scale,shape,h)
    real(dp),intent(in)::excess(:),scale,shape
    real(dp),intent(out)::h(2,2)
    real(dp)::hs,hx,f0,fpp,fpm,fmp,fmm
    hs=max(1.0e-5_dp*scale,1.0e-7_dp);hx=1.0e-5_dp
    f0=-gpd_log_likelihood(excess,scale,shape)
    h(1,1)=(-gpd_log_likelihood(excess,scale+hs,shape)-2.0_dp*f0- &
      gpd_log_likelihood(excess,max(scale-hs,scale*0.5_dp),shape))/(hs*hs)
    h(2,2)=(-gpd_log_likelihood(excess,scale,shape+hx)-2.0_dp*f0- &
      gpd_log_likelihood(excess,scale,shape-hx))/(hx*hx)
    fpp=-gpd_log_likelihood(excess,scale+hs,shape+hx)
    fpm=-gpd_log_likelihood(excess,scale+hs,shape-hx)
    fmp=-gpd_log_likelihood(excess,max(scale-hs,scale*0.5_dp),shape+hx)
    fmm=-gpd_log_likelihood(excess,max(scale-hs,scale*0.5_dp),shape-hx)
    h(1,2)=(fpp-fpm-fmp+fmm)/(4.0_dp*hs*hx);h(2,1)=h(1,2)
  end subroutine gpd_hessian

  subroutine gpd_fit(returns,p,result,p_threshold,threshold,confidence)
    real(dp),intent(in)::returns(:),p
    type(gpd_fit_result),intent(out)::result
    real(dp),intent(in),optional::p_threshold,threshold,confidence
    real(dp),allocatable::loss(:),excess(:)
    real(dp)::ptr,u,theta0(2),theta(2),obj,shape,scale,h(2,2),det,cov(2,2)
    real(dp)::gv(2),ge(2),hs,hx,vp,vm,ep,em,zcrit,conf
    integer::n,m,it
    logical::conv
    n=size(returns);allocate(loss(n));loss=-100.0_dp*returns
    ptr=0.97_dp;if(present(p_threshold))ptr=p_threshold
    if(present(threshold))then;u=100.0_dp*threshold;else;u=normal_quantile(ptr)*sd_value(loss);end if
    m=count(loss>u);result%observations=n;result%exceedances=m;result%threshold=u/100.0_dp
    if(m<3)then;result%converged=.false.;return;end if
    allocate(excess(m));excess=pack(loss-u,loss>u)
    theta0=[log(max(mean_value(excess),1.0e-3_dp)),atanh(0.3_dp/0.95_dp)]
    call nelder_mead_gpd(excess,theta0,theta,obj,it,conv)
    scale=exp(theta(1));shape=0.95_dp*tanh(theta(2))
    result%scale=scale/100.0_dp;result%shape=shape;result%exceedance_rate=real(m,dp)/real(n,dp)
    result%log_likelihood=-obj;result%iterations=it;result%converged=conv
    result%var_value=gpd_var_value(u,scale,shape,1.0_dp-p,result%exceedance_rate)/100.0_dp
    result%es_value=gpd_es_value(u,scale,shape,100.0_dp*result%var_value)/100.0_dp
    call gpd_hessian(excess,scale,shape,h);det=h(1,1)*h(2,2)-h(1,2)*h(2,1)
    if(det>1.0e-18_dp .and. h(1,1)>0.0_dp .and. h(2,2)>0.0_dp)then
      cov(1,1)=h(2,2)/det;cov(2,2)=h(1,1)/det;cov(1,2)=-h(1,2)/det;cov(2,1)=cov(1,2)
      hs=max(1.0e-5_dp*scale,1.0e-7_dp);hx=1.0e-5_dp
      vp=gpd_var_value(u,scale+hs,shape,1.0_dp-p,result%exceedance_rate)
      vm=gpd_var_value(u,max(scale-hs,scale*0.5_dp),shape,1.0_dp-p,result%exceedance_rate)
      gv(1)=(vp-vm)/(2.0_dp*hs)
      vp=gpd_var_value(u,scale,shape+hx,1.0_dp-p,result%exceedance_rate)
      vm=gpd_var_value(u,scale,shape-hx,1.0_dp-p,result%exceedance_rate);gv(2)=(vp-vm)/(2.0_dp*hx)
      ep=gpd_es_value(u,scale+hs,shape,gpd_var_value(u,scale+hs,shape,1.0_dp-p,result%exceedance_rate))
      em=gpd_es_value(u,max(scale-hs,scale*0.5_dp),shape,gpd_var_value(u,max(scale-hs,scale*0.5_dp),shape,1.0_dp-p,result%exceedance_rate));ge(1)=(ep-em)/(2.0_dp*hs)
      ep=gpd_es_value(u,scale,shape+hx,gpd_var_value(u,scale,shape+hx,1.0_dp-p,result%exceedance_rate))
      em=gpd_es_value(u,scale,shape-hx,gpd_var_value(u,scale,shape-hx,1.0_dp-p,result%exceedance_rate));ge(2)=(ep-em)/(2.0_dp*hx)
      result%var_standard_error=sqrt(max(dot_product(gv,matmul(cov,gv)),0.0_dp))/100.0_dp
      result%es_standard_error=sqrt(max(dot_product(ge,matmul(cov,ge)),0.0_dp))/100.0_dp
    end if
    conf=0.95_dp;if(present(confidence))conf=confidence
    zcrit=normal_quantile(0.5_dp+0.5_dp*conf)
    result%var_lower=max(0.0_dp,result%var_value-zcrit*result%var_standard_error)
    result%var_upper=result%var_value+zcrit*result%var_standard_error
    result%es_lower=max(0.0_dp,result%es_value-zcrit*result%es_standard_error)
    result%es_upper=result%es_value+zcrit*result%es_standard_error
  end subroutine gpd_fit

  subroutine monte_carlo_asset_risk(r,p,nsim,seed,var_values,es_values,ok)
    real(dp),intent(in)::r(:,:),p
    integer,intent(in)::nsim
    integer(kind=8),intent(in)::seed
    real(dp),intent(out)::var_values(:),es_values(:)
    logical,intent(out)::ok
    real(dp),allocatable::mu(:),cov(:,:),draws(:,:)
    integer::j,q
    q=size(r,2);allocate(mu(q),cov(q,q));do j=1,q;mu(j)=mean_value(r(:,j));end do
    call covariance_matrix(r,cov,.true.);call multivariate_normal_draws(mu,cov,nsim,seed,draws,ok)
    if(.not.ok)return
    do j=1,min(q,size(var_values),size(es_values))
      var_values(j)=historical_var(draws(:,j),p);es_values(j)=historical_es(draws(:,j),p)
    end do
  end subroutine monte_carlo_asset_risk

  subroutine monte_carlo_portfolio_risk(r,weights,p,nsim,seed,var_value,es_value,var_components,es_components,ok)
    real(dp),intent(in)::r(:,:),weights(:),p
    integer,intent(in)::nsim
    integer(kind=8),intent(in)::seed
    real(dp),intent(out)::var_value,es_value,var_components(:),es_components(:)
    logical,intent(out)::ok
    real(dp),allocatable::mu(:),cov(:,:),draws(:,:),port(:)
    real(dp)::cutoff
    integer::j,n_tail,q
    q=size(r,2);allocate(mu(q),cov(q,q));do j=1,q;mu(j)=mean_value(r(:,j));end do
    call covariance_matrix(r,cov,.true.);call multivariate_normal_draws(mu,cov,nsim,seed,draws,ok)
    if(.not.ok)return
    allocate(port(nsim));port=matmul(draws,weights);var_value=historical_var(port,p);es_value=historical_es(port,p)
    call component_var(draws,weights,p,'historical',var_components,var_value)
    cutoff=quantile_type7(port,1.0_dp-p);n_tail=count(port<=cutoff)
    es_components=0.0_dp
    if(n_tail>0)then
      do j=1,min(q,size(es_components));es_components(j)=-weights(j)*sum(draws(:,j),mask=port<=cutoff)/real(n_tail,dp);end do
    end if
  end subroutine monte_carlo_portfolio_risk


  pure real(dp) function lognormal_var(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp),allocatable::y(:)
    allocate(y(size(r)))
    y=log(max(1.0_dp+r,tiny(1.0_dp)))
    v=1.0_dp-exp(mean_value(y)-sd_value(y)*normal_quantile(p))
  end function lognormal_var

  pure real(dp) function lognormal_es(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp),allocatable::y(:)
    real(dp)::mu,sigma,alpha,ratio
    allocate(y(size(r)));y=log(max(1.0_dp+r,tiny(1.0_dp)))
    mu=mean_value(y);sigma=sd_value(y);alpha=1.0_dp-p
    ratio=exp(mu+0.5_dp*sigma*sigma)*normal_cdf(normal_quantile(alpha)-sigma)/max(alpha,tiny(1.0_dp))
    v=1.0_dp-ratio
  end function lognormal_es

  subroutine kernel_portfolio_risk(r,weights,p,var_value,es_value,var_components,es_components)
    real(dp),intent(in)::r(:,:),weights(:),p
    real(dp),intent(out)::var_value,es_value,var_components(:),es_components(:)
    real(dp),allocatable::port(:),kw(:)
    real(dp)::bandwidth,cutoff,totalw,rawsum
    integer::n,j,ntail
    n=size(r,1);allocate(port(n),kw(n));port=matmul(r,weights)
    var_value=historical_var(port,p);cutoff=-var_value
    bandwidth=2.575_dp*sd_value(port)/max(real(n,dp)**0.2_dp,1.0_dp)
    if(bandwidth<=tiny(1.0_dp))bandwidth=max(sd_value(port),1.0e-8_dp)
    kw=max(0.0_dp,1.0_dp-abs(cutoff-port)/bandwidth);totalw=sum(kw)
    var_components=0.0_dp
    if(totalw>0.0_dp)then
      do j=1,min(size(weights),size(var_components))
        var_components(j)=-weights(j)*sum(kw*r(:,j))/totalw
      end do
      rawsum=sum(var_components)
      if(abs(rawsum)>tiny(1.0_dp))var_components=var_components*var_value/rawsum
    end if
    ntail=count(port<=cutoff)
    if(ntail>0)then
      es_value=-sum(port,mask=port<=cutoff)/real(ntail,dp)
      do j=1,min(size(weights),size(es_components))
        es_components(j)=-weights(j)*sum(r(:,j),mask=port<=cutoff)/real(ntail,dp)
      end do
    else
      es_value=var_value;es_components=var_components
    end if
  end subroutine kernel_portfolio_risk

  subroutine bootstrap_risk_standard_errors(r,p,nboot,block_length,seed,var_se,es_se,sharpe_se)
    real(dp),intent(in)::r(:),p
    integer,intent(in)::nboot,block_length
    integer(kind=8),intent(in)::seed
    real(dp),intent(out)::var_se,es_se,sharpe_se
    real(dp),allocatable::sample(:),v(:),e(:),s(:)
    type(rng_state)::rng
    integer::n,b,i,start,len,j,idx
    real(dp)::mv,me,ms
    n=size(r);allocate(sample(n),v(nboot),e(nboot),s(nboot));call rng_seed(rng,seed)
    do b=1,nboot
      i=1
      do while(i<=n)
        start=rng_integer(rng,1,n);len=min(max(block_length,1),n-i+1)
        do j=1,len
          idx=1+mod(start+j-2,n)
          sample(i+j-1)=r(idx)
        end do
        i=i+len
      end do
      v(b)=historical_var(sample,p);e(b)=historical_es(sample,p);s(b)=mean_value(sample)/max(sd_value(sample),tiny(1.0_dp))
    end do
    mv=mean_value(v);me=mean_value(e);ms=mean_value(s)
    var_se=sqrt(sum((v-mv)**2)/real(max(nboot-1,1),dp));es_se=sqrt(sum((e-me)**2)/real(max(nboot-1,1),dp));sharpe_se=sqrt(sum((s-ms)**2)/real(max(nboot-1,1),dp))
  end subroutine bootstrap_risk_standard_errors
end module tail_models_mod
