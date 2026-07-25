! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module sv_core
  use sv_kinds, only : dp, pi, tiny_var
  use sv_rng, only : randu,randn,rand_inv_gamma,rand_student_t
  use sv_stats, only : normal_logpdf,student_logpdf_std,clamp,safe_logsumexp
  use sv_types
  implicit none
  private
  public :: omori_prob,omori_mean,omori_var,get_omori_constants
  public :: simulate_sv,sv_complete_loglik,svsample,predict_sv,sv_residuals
  public :: svsample_t,svsample_leverage,svsample_t_leverage,bayesian_regression_update
  public :: rolling_sv_forecast,sv_single_sweep
  real(dp),parameter :: omori_prob(10)=[.00609_dp,.04775_dp,.13057_dp,.20674_dp,.22715_dp,.18842_dp,.12047_dp,.05591_dp,.01575_dp,.00115_dp]
  real(dp),parameter :: omori_mean(10)=[1.92677_dp,1.34744_dp,.73504_dp,.02266_dp,-.85173_dp,-1.97278_dp,-3.46788_dp,-5.55246_dp,-8.68384_dp,-14.65000_dp]
  real(dp),parameter :: omori_var(10)=[.11265_dp,.17788_dp,.26768_dp,.40611_dp,.62699_dp,.98583_dp,1.57469_dp,2.54498_dp,4.16591_dp,7.33342_dp]
contains
  subroutine get_omori_constants(prob,mean,var)
    real(dp),intent(out)::prob(10),mean(10),var(10);prob=omori_prob;mean=omori_mean;var=omori_var
  end subroutine get_omori_constants

  subroutine simulate_sv(n,params,out,x)
    integer,intent(in)::n
    type(sv_params),intent(in)::params
    type(sv_sim_result),intent(out)::out
    real(dp),intent(in),optional::x(:,:)
    real(dp)::eta,eps,xb,scale
    integer::t
    allocate(out%y(n),out%latent(n),out%vol(n),out%tau(n));out%params=params
    out%latent0=params%mu+params%sigma/sqrt(max(tiny_var,1.0_dp-params%phi**2))*randn()
    do t=1,n
      eta=randn()
      if(t==1) then
        out%latent(t)=params%mu+params%phi*(out%latent0-params%mu)+params%sigma*eta
      else
        out%latent(t)=params%mu+params%phi*(out%latent(t-1)-params%mu)+params%sigma*eta
      end if
      if(params%nu<huge(1.0_dp)/10.0_dp) then
        out%tau(t)=rand_inv_gamma(0.5_dp*params%nu,0.5_dp*(params%nu-2.0_dp))
      else
        out%tau(t)=1.0_dp
      end if
      if(t<n) then
        ! eta_t drives h_{t+1}; generate it on the next pass, so use an equivalent correlated pair here
        eps=randn()
      else
        eps=randn()
      end if
      xb=0.0_dp
      if(present(x).and.allocated(params%beta)) xb=dot_product(x(t,:),params%beta)
      scale=exp(0.5_dp*out%latent(t))*sqrt(out%tau(t))
      out%y(t)=xb+scale*eps
      out%vol(t)=scale
    end do
    ! Re-simulate exact leverage-coupled path when rho is nonzero.
    if(abs(params%rho)>1.0e-14_dp) call simulate_sv_leverage(n,params,out,x)
    out%vol0=exp(0.5_dp*out%latent0)
  end subroutine simulate_sv

  subroutine simulate_sv_leverage(n,params,out,x)
    integer,intent(in)::n
    type(sv_params),intent(in)::params
    type(sv_sim_result),intent(inout)::out
    real(dp),intent(in),optional::x(:,:)
    real(dp),allocatable::eta(:),z(:)
    real(dp)::xb
    integer::t
    allocate(eta(n),z(n));do t=1,n;eta(t)=randn();z(t)=randn();end do
    out%latent0=params%mu+params%sigma/sqrt(max(tiny_var,1.0_dp-params%phi**2))*randn()
    out%latent(1)=params%mu+params%phi*(out%latent0-params%mu)+params%sigma*randn()
    do t=1,n-1
      out%latent(t+1)=params%mu+params%phi*(out%latent(t)-params%mu)+params%sigma*eta(t)
    end do
    do t=1,n
      if(params%nu<huge(1.0_dp)/10.0_dp) then
        out%tau(t)=rand_inv_gamma(.5_dp*params%nu,.5_dp*(params%nu-2.0_dp))
      else; out%tau(t)=1.0_dp; end if
      xb=0.0_dp;if(present(x).and.allocated(params%beta))xb=dot_product(x(t,:),params%beta)
      if(t<n) then
        out%y(t)=xb+exp(.5_dp*out%latent(t))*sqrt(out%tau(t))*(params%rho*eta(t)+sqrt(1.0_dp-params%rho**2)*z(t))
      else
        out%y(t)=xb+exp(.5_dp*out%latent(t))*sqrt(out%tau(t))*z(t)
      end if
      out%vol(t)=exp(.5_dp*out%latent(t))*sqrt(out%tau(t))
    end do
  end subroutine simulate_sv_leverage

  real(dp) function sv_complete_loglik(y,h,h0,p,tau,x,include_prior,prior) result(ll)
    real(dp),intent(in)::y(:),h(:),h0
    type(sv_params),intent(in)::p
    real(dp),intent(in),optional::tau(:),x(:,:)
    logical,intent(in),optional::include_prior
    type(sv_prior),intent(in),optional::prior
    real(dp)::tt,xb,scale,eta,sd,pr
    integer::t,n
    logical::ip
    n=size(y);ll=normal_logpdf(h0,p%mu,p%sigma/sqrt(max(tiny_var,1.0_dp-p%phi*p%phi)))
    ll=ll+normal_logpdf(h(1),p%mu+p%phi*(h0-p%mu),p%sigma)
    do t=2,n;ll=ll+normal_logpdf(h(t),p%mu+p%phi*(h(t-1)-p%mu),p%sigma);end do
    do t=1,n
      tt=1.0_dp;if(present(tau))tt=tau(t)
      xb=0.0_dp;if(present(x).and.allocated(p%beta))xb=dot_product(x(t,:),p%beta)
      scale=exp(0.5_dp*h(t))*sqrt(tt)
      if(t<n .and. abs(p%rho)>0.0_dp) then
        eta=(h(t+1)-p%mu-p%phi*(h(t)-p%mu))/p%sigma
        sd=scale*sqrt(max(tiny_var,1.0_dp-p%rho*p%rho))
        ll=ll+normal_logpdf(y(t),xb+scale*p%rho*eta,sd)
      else
        ll=ll+normal_logpdf(y(t),xb,scale)
      end if
      if(p%nu<huge(1.0_dp)/10.0_dp .and. present(tau)) then
        ll=ll+(.5_dp*p%nu)*log(.5_dp*(p%nu-2.0_dp))-log_gamma(.5_dp*p%nu) &
          +(-.5_dp*p%nu-1.0_dp)*log(tt)-.5_dp*(p%nu-2.0_dp)/tt
      end if
    end do
    ip=.false.;if(present(include_prior))ip=include_prior
    if(ip.and.present(prior)) then
      pr=normal_logpdf(p%mu,prior%mu_mean,prior%mu_sd)
      pr=pr+(prior%phi_a-1.0_dp)*log(max(tiny_var,.5_dp*(p%phi+1.0_dp)))+(prior%phi_b-1.0_dp)*log(max(tiny_var,.5_dp*(1.0_dp-p%phi)))
      pr=pr+(prior%sigma_shape-1.0_dp)*log(p%sigma*p%sigma)-prior%sigma_rate*p%sigma*p%sigma+log(2.0_dp*p%sigma)
      pr=pr+normal_logpdf(atanh(clamp(p%rho,-.999999_dp,.999999_dp)),0.0_dp,prior%rho_sd)
      if(p%nu<huge(1.0_dp)/10.0_dp)pr=pr-prior%nu_rate*(p%nu-2.0_dp)
      if(allocated(p%beta))pr=pr-sum(p%beta*p%beta)/(2.0_dp*prior%beta_sd**2)
      ll=ll+pr
    end if
  end function sv_complete_loglik

  subroutine sample_mixture_indicators(y,h,tau,x,p,r,offset)
    real(dp),intent(in)::y(:),h(:),tau(:),offset
    real(dp),intent(in),optional::x(:,:)
    type(sv_params),intent(in)::p
    integer,intent(out)::r(:)
    real(dp)::ystar,xb,lp(10),u,c
    integer::t,j
    do t=1,size(y)
      xb=0.0_dp;if(present(x).and.allocated(p%beta))xb=dot_product(x(t,:),p%beta)
      ystar=log((y(t)-xb)**2/max(tau(t),tiny_var)+offset)-h(t)
      do j=1,10;lp(j)=log(omori_prob(j))+normal_logpdf(ystar,omori_mean(j),sqrt(omori_var(j)));end do
      lp=exp(lp-safe_logsumexp(lp));u=randu();c=0.0_dp;r(t)=10
      do j=1,10;c=c+lp(j);if(u<=c)then;r(t)=j;exit;end if;end do
    end do
  end subroutine sample_mixture_indicators

  subroutine draw_latent_mixture(y,p,tau,x,h,h0,offset)
    real(dp),intent(in)::y(:),tau(:),offset
    type(sv_params),intent(in)::p
    real(dp),intent(in),optional::x(:,:)
    real(dp),intent(inout)::h(:),h0
    integer,allocatable::r(:)
    real(dp),allocatable::diag(:),off(:),cvec(:),ld(:),lo(:),z(:),tmp(:),meanv(:)
    real(dp)::s2,phi2,xb,ystar
    integer::n,t
    n=size(y);allocate(r(n));call sample_mixture_indicators(y,h,tau,x,p,r,offset)
    allocate(diag(n+1),off(n),cvec(n+1),ld(n+1),lo(n),z(n+1),tmp(n+1),meanv(n+1))
    s2=p%sigma**2;phi2=p%phi**2
    diag(1)=(1.0_dp-p%phi**2+phi2)/s2
    cvec(1)=p%mu*(1.0_dp-p%phi**2-p%phi*(1.0_dp-p%phi))/s2
    do t=2,n
      xb=0.0_dp;if(present(x).and.allocated(p%beta))xb=dot_product(x(t-1,:),p%beta)
      ystar=log((y(t-1)-xb)**2/max(tau(t-1),tiny_var)+offset)
      diag(t)=1.0_dp/omori_var(r(t-1))+(1.0_dp+phi2)/s2
      cvec(t)=(ystar-omori_mean(r(t-1)))/omori_var(r(t-1))+p%mu*(1.0_dp-p%phi)**2/s2
    end do
    xb=0.0_dp;if(present(x).and.allocated(p%beta))xb=dot_product(x(n,:),p%beta)
    ystar=log((y(n)-xb)**2/max(tau(n),tiny_var)+offset)
    diag(n+1)=1.0_dp/omori_var(r(n))+1.0_dp/s2
    cvec(n+1)=(ystar-omori_mean(r(n)))/omori_var(r(n))+p%mu*(1.0_dp-p%phi)/s2
    off=-p%phi/s2
    ld(1)=sqrt(max(tiny_var,diag(1)))
    do t=1,n;lo(t)=off(t)/ld(t);ld(t+1)=sqrt(max(tiny_var,diag(t+1)-lo(t)**2));end do
    tmp(1)=cvec(1)/ld(1);do t=2,n+1;tmp(t)=(cvec(t)-lo(t-1)*tmp(t-1))/ld(t);end do
    meanv(n+1)=tmp(n+1)/ld(n+1);do t=n,1,-1;meanv(t)=(tmp(t)-lo(t)*meanv(t+1))/ld(t);end do
    do t=1,n+1;tmp(t)=randn();end do
    z(n+1)=tmp(n+1)/ld(n+1)
    do t=n,1,-1;z(t)=(tmp(t)-lo(t)*z(t+1))/ld(t);end do
    h0=meanv(1)+z(1);h=meanv(2:n+1)+z(2:n+1)
  end subroutine draw_latent_mixture

  subroutine update_latent_rw(y,p,tau,x,h,h0,prior,opt,accepted)
    real(dp),intent(in)::y(:),tau(:)
    type(sv_params),intent(in)::p
    real(dp),intent(in),optional::x(:,:)
    real(dp),intent(inout)::h(:),h0
    type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt
    integer,intent(inout)::accepted
    real(dp)::old,prop,lold,lnew,prior_touch
    integer::t,s
    prior_touch=prior%mu_mean*0.0_dp
    do s=1,opt%latent_sweeps
      old=h0;prop=old+opt%step_latent*randn();lold=sv_complete_loglik(y,h,h0,p,tau,x);h0=prop;lnew=sv_complete_loglik(y,h,h0,p,tau,x)
      if(log(randu())<lnew-lold)then;accepted=accepted+1;else;h0=old;end if
      do t=1,size(h)
        old=h(t);prop=old+opt%step_latent*randn();lold=sv_complete_loglik(y,h,h0,p,tau,x);h(t)=prop;lnew=sv_complete_loglik(y,h,h0,p,tau,x)
        if(log(randu())<lnew-lold)then;accepted=accepted+1;else;h(t)=old;end if
      end do
    end do
  end subroutine update_latent_rw

  subroutine mh_scalar(y,h,h0,p,tau,x,prior,which,step,accepted)
    real(dp),intent(in)::y(:),h(:),h0,tau(:),step
    type(sv_params),intent(inout)::p
    real(dp),intent(in),optional::x(:,:)
    type(sv_prior),intent(in)::prior
    integer,intent(in)::which
    integer,intent(inout)::accepted
    type(sv_params)::q
    real(dp)::oldz,newz,lold,lnew,jold,jnew
    q=p;jold=0.0_dp;jnew=0.0_dp
    select case(which)
    case(1);q%mu=p%mu+step*randn()
    case(2);oldz=atanh(clamp(p%phi,-.999999_dp,.999999_dp));newz=oldz+step*randn();q%phi=tanh(newz);jold=log(1.0_dp-p%phi**2);jnew=log(1.0_dp-q%phi**2)
    case(3);oldz=log(p%sigma);newz=oldz+step*randn();q%sigma=exp(newz);jold=oldz;jnew=newz
    case(4);oldz=atanh(clamp(p%rho,-.999999_dp,.999999_dp));newz=oldz+step*randn();q%rho=tanh(newz);jold=log(1.0_dp-p%rho**2);jnew=log(1.0_dp-q%rho**2)
    case(5);oldz=log(p%nu-2.0_dp);newz=oldz+step*randn();q%nu=2.0_dp+exp(newz);jold=oldz;jnew=newz
    end select
    lold=sv_complete_loglik(y,h,h0,p,tau,x,.true.,prior)+jold
    lnew=sv_complete_loglik(y,h,h0,q,tau,x,.true.,prior)+jnew
    if(log(randu())<lnew-lold)then;p=q;accepted=accepted+1;end if
  end subroutine mh_scalar

  subroutine update_beta_rw(y,h,h0,p,tau,x,prior,opt,accepted)
    real(dp),intent(in)::y(:),h(:),h0,tau(:),x(:,:)
    type(sv_params),intent(inout)::p
    type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt
    integer,intent(inout)::accepted
    integer::j
    type(sv_params)::q
    real(dp)::lold,lnew
    if(.not.allocated(p%beta))return
    do j=1,size(p%beta)
      q=p;q%beta(j)=p%beta(j)+opt%step_beta*randn()
      lold=sv_complete_loglik(y,h,h0,p,tau,x,.true.,prior);lnew=sv_complete_loglik(y,h,h0,q,tau,x,.true.,prior)
      if(log(randu())<lnew-lold)then;p=q;accepted=accepted+1;end if
    end do
  end subroutine update_beta_rw

  subroutine update_tau(y,h,p,tau,x)
    real(dp),intent(in)::y(:),h(:)
    type(sv_params),intent(in)::p
    real(dp),intent(inout)::tau(:)
    real(dp),intent(in),optional::x(:,:)
    real(dp)::xb,u,e2,proposal,shape,scale,eta,mn,sd,lar
    integer::t
    if(p%nu>=huge(1.0_dp)/10.0_dp)return
    shape=.5_dp*(p%nu+1.0_dp)
    do t=1,size(y)
      xb=0.0_dp;if(present(x).and.allocated(p%beta))xb=dot_product(x(t,:),p%beta)
      u=(y(t)-xb)/exp(.5_dp*h(t));e2=u*u;scale=.5_dp*(p%nu-2.0_dp+e2)
      proposal=rand_inv_gamma(shape,scale)
      if(abs(p%rho)>1.0e-14_dp.and.t<size(y))then
        eta=(h(t+1)-p%mu-p%phi*(h(t)-p%mu))/p%sigma
        mn=p%rho*eta;sd=sqrt(max(tiny_var,1.0_dp-p%rho*p%rho))
        lar=(normal_logpdf(u,sqrt(proposal)*mn,sqrt(proposal)*sd)+log_inv_gamma(proposal,.5_dp*p%nu,.5_dp*(p%nu-2.0_dp))-log_inv_gamma(proposal,shape,scale)) &
          -(normal_logpdf(u,sqrt(tau(t))*mn,sqrt(tau(t))*sd)+log_inv_gamma(tau(t),.5_dp*p%nu,.5_dp*(p%nu-2.0_dp))-log_inv_gamma(tau(t),shape,scale))
        if(log(randu())<lar)tau(t)=proposal
      else
        tau(t)=proposal
      end if
    end do
  end subroutine update_tau

  pure real(dp) function log_inv_gamma(x,shape,scale) result(v)
    real(dp),intent(in)::x,shape,scale
    v=shape*log(scale)-log_gamma(shape)-(shape+1.0_dp)*log(x)-scale/x
  end function log_inv_gamma

  subroutine sv_single_sweep(y,p,h,h0,tau,prior,opt,x,accept)
    real(dp),intent(in)::y(:)
    type(sv_params),intent(inout)::p
    real(dp),intent(inout)::h(:),h0,tau(:)
    type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt
    real(dp),intent(in),optional::x(:,:)
    integer,intent(inout)::accept(:)
    if(opt%use_mixture.and..not.opt%sample_rho.and.abs(p%rho)<1.0e-14_dp)then
      call draw_latent_mixture(y,p,tau,x,h,h0,opt%offset)
    else
      call update_latent_rw(y,p,tau,x,h,h0,prior,opt,accept(6))
    end if
    call update_tau(y,h,p,tau,x)
    if(opt%sample_mu)call mh_scalar(y,h,h0,p,tau,x,prior,1,opt%step_mu,accept(1))
    if(opt%sample_phi)call mh_scalar(y,h,h0,p,tau,x,prior,2,opt%step_phi,accept(2))
    if(opt%sample_sigma)call mh_scalar(y,h,h0,p,tau,x,prior,3,opt%step_logsigma,accept(3))
    if(opt%sample_rho)call mh_scalar(y,h,h0,p,tau,x,prior,4,opt%step_rho,accept(4))
    if(opt%sample_nu)call mh_scalar(y,h,h0,p,tau,x,prior,5,opt%step_lognu,accept(5))
    if(opt%sample_beta.and.present(x))call update_beta_rw(y,h,h0,p,tau,x,prior,opt,accept(7))
  end subroutine sv_single_sweep

  subroutine svsample(y,start,prior,opt,out,x)
    real(dp),intent(in)::y(:)
    type(sv_params),intent(in)::start
    type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt
    type(sv_draws),intent(out)::out
    real(dp),intent(in),optional::x(:,:)
    type(sv_params)::p
    real(dp),allocatable::h(:),tau(:)
    real(dp)::h0
    integer::iter,total,keep,pdim,acc(7)
    p=start;pdim=0;if(allocated(p%beta))pdim=size(p%beta)
    allocate(h(size(y)),tau(size(y)));h=log(y*y+opt%offset);h0=sum(h)/real(size(h),dp);tau=1.0_dp
    total=opt%burnin+opt%draws*opt%thin
    out%nobs=size(y);out%ndraws=opt%draws;allocate(out%data(size(y)));out%data=y;if(present(x))then;allocate(out%design(size(x,1),size(x,2)));out%design=x;end if
    allocate(out%para(5,opt%draws),out%latent0(opt%draws),out%accept(7));out%para=0.0_dp
    if(opt%store_latent)allocate(out%latent(size(y),opt%draws))
    if(opt%store_tau)allocate(out%tau(size(y),opt%draws))
    if(pdim>0)allocate(out%beta(pdim,opt%draws))
    keep=0;acc=0
    do iter=1,total
      call sv_single_sweep(y,p,h,h0,tau,prior,opt,x,acc)
      if(iter>opt%burnin.and.mod(iter-opt%burnin,opt%thin)==0)then
        keep=keep+1;out%para(:,keep)=[p%mu,p%phi,p%sigma,p%nu,p%rho];out%latent0(keep)=h0
        if(opt%store_latent)out%latent(:,keep)=h
        if(opt%store_tau)out%tau(:,keep)=tau
        if(pdim>0)out%beta(:,keep)=p%beta
      end if
    end do
    out%accept=real(acc,dp)/real(max(1,total),dp)
  end subroutine svsample

  subroutine predict_sv(draws,steps,each,pred,xnew)
    type(sv_draws),intent(in)::draws
    integer,intent(in)::steps,each
    type(sv_prediction),intent(out)::pred
    real(dp),intent(in),optional::xnew(:,:)
    integer::d,e,s,n,pdim
    real(dp)::hcur,mu,phi,sigma,nu,rho,tau,eps,xb,resilast,taulast,hlast
    n=draws%ndraws;pdim=0;if(allocated(draws%beta))pdim=size(draws%beta,1)
    allocate(pred%y(steps,n,each),pred%latent(steps,n,each),pred%vola(steps,n,each))
    do e=1,each;do d=1,n
      mu=draws%para(1,d);phi=draws%para(2,d);sigma=draws%para(3,d);nu=draws%para(4,d);rho=draws%para(5,d)
      if(allocated(draws%latent))then;hlast=draws%latent(size(draws%latent,1),d);else;hlast=mu;end if
      taulast=1.0_dp;if(allocated(draws%tau))taulast=draws%tau(size(draws%tau,1),d)
      xb=0.0_dp
      if(pdim>0.and.allocated(draws%design))xb=dot_product(draws%design(size(draws%design,1),:),draws%beta(:,d))
      resilast=0.0_dp
      if(allocated(draws%data))resilast=(draws%data(size(draws%data))-xb)/exp(.5_dp*hlast)/sqrt(taulast)
      hcur=mu+phi*(hlast-mu)+sigma*(rho*resilast+sqrt(max(tiny_var,1.0_dp-rho*rho))*randn())
      do s=1,steps
        tau=1.0_dp;if(nu<huge(1.0_dp)/10.0_dp)tau=rand_inv_gamma(.5_dp*nu,.5_dp*(nu-2.0_dp))
        xb=0.0_dp;if(pdim>0.and.present(xnew))xb=dot_product(xnew(s,:),draws%beta(:,d))
        eps=randn();pred%latent(s,d,e)=hcur;pred%vola(s,d,e)=exp(.5_dp*hcur)*sqrt(tau);pred%y(s,d,e)=xb+pred%vola(s,d,e)*eps
        if(s<steps)hcur=mu+phi*(hcur-mu)+sigma*(rho*eps+sqrt(max(tiny_var,1.0_dp-rho*rho))*randn())
      end do
    end do;end do
  end subroutine predict_sv

  subroutine sv_residuals(y,draws,resid,x)
    real(dp),intent(in)::y(:)
    type(sv_draws),intent(in)::draws
    real(dp),intent(out)::resid(:)
    real(dp),intent(in),optional::x(:,:)
    real(dp),allocatable::v(:)
    real(dp)::xb
    integer::t,d
    if(.not.allocated(draws%latent))error stop 'sv_residuals: latent draws required'
    allocate(v(size(y)));v=0.0_dp
    do d=1,draws%ndraws;v=v+exp(draws%latent(:,d));end do;v=v/real(draws%ndraws,dp)
    do t=1,size(y);xb=0.0_dp;if(present(x).and.allocated(draws%beta))xb=dot_product(x(t,:),sum(draws%beta,dim=2)/real(draws%ndraws,dp));resid(t)=(y(t)-xb)/sqrt(v(t));end do
  end subroutine sv_residuals


  subroutine svsample_t(y,start,prior,opt,out,x)
    real(dp),intent(in)::y(:);type(sv_params),intent(in)::start;type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt;type(sv_draws),intent(out)::out;real(dp),intent(in),optional::x(:,:)
    type(sv_params)::p;type(sv_mcmc_options)::o
    p=start;o=opt;if(p%nu>=huge(1.0_dp)/10.0_dp)p%nu=10.0_dp;o%sample_nu=.true.;o%store_tau=.true.;call svsample(y,p,prior,o,out,x)
  end subroutine svsample_t

  subroutine svsample_leverage(y,start,prior,opt,out,x)
    real(dp),intent(in)::y(:);type(sv_params),intent(in)::start;type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt;type(sv_draws),intent(out)::out;real(dp),intent(in),optional::x(:,:)
    type(sv_params)::p;type(sv_mcmc_options)::o
    p=start;o=opt;o%sample_rho=.true.;o%use_mixture=.false.;call svsample(y,p,prior,o,out,x)
  end subroutine svsample_leverage

  subroutine svsample_t_leverage(y,start,prior,opt,out,x)
    real(dp),intent(in)::y(:);type(sv_params),intent(in)::start;type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt;type(sv_draws),intent(out)::out;real(dp),intent(in),optional::x(:,:)
    type(sv_params)::p;type(sv_mcmc_options)::o
    p=start;o=opt;if(p%nu>=huge(1.0_dp)/10.0_dp)p%nu=10.0_dp;o%sample_nu=.true.;o%sample_rho=.true.;o%use_mixture=.false.;o%store_tau=.true.;call svsample(y,p,prior,o,out,x)
  end subroutine svsample_t_leverage

  subroutine bayesian_regression_update(y,x,prior_mean,prior_precision,beta)
    use sv_linalg, only : inverse_spd,mvn_draw
    real(dp),intent(in)::y(:),x(:,:),prior_mean(:),prior_precision(:,:)
    real(dp),intent(out)::beta(:)
    real(dp),allocatable::precision(:,:),cov(:,:),mean(:)
    integer::info
    allocate(precision(size(beta),size(beta)),mean(size(beta)));precision=matmul(transpose(x),x)+prior_precision
    call inverse_spd(precision,cov,info);if(info/=0)error stop 'bayesian_regression_update: precision failure'
    mean=matmul(cov,matmul(transpose(x),y)+matmul(prior_precision,prior_mean));call mvn_draw(mean,cov,beta,info)
    if(info/=0)error stop 'bayesian_regression_update: draw failure'
  end subroutine bayesian_regression_update
  subroutine rolling_sv_forecast(y,window,start,prior,opt,variance_forecast)
    real(dp),intent(in)::y(:)
    integer,intent(in)::window
    type(sv_params),intent(in)::start
    type(sv_prior),intent(in)::prior
    type(sv_mcmc_options),intent(in)::opt
    real(dp),allocatable,intent(out)::variance_forecast(:)
    type(sv_draws)::dr
    integer::t,n,d
    real(dp)::v
    n=size(y)-window;allocate(variance_forecast(n))
    do t=1,n
      call svsample(y(t:t+window-1),start,prior,opt,dr)
      v=0.0_dp
      do d=1,dr%ndraws
        v=v+exp(dr%para(1,d)+dr%para(2,d)*(dr%latent(window,d)-dr%para(1,d))+.5_dp*dr%para(3,d)**2)
      end do
      variance_forecast(t)=v/real(dr%ndraws,dp)
    end do
  end subroutine rolling_sv_forecast
end module sv_core
