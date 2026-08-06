! SPDX-License-Identifier: GPL-2.0-only
module tsgarch_simulation_module
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ghyp_kinds, only : dp, i8
  use tsd_distributions, only : rng_state, seed_rng, rdist
  use tsgarch_types
  use tsgarch_model, only : filter_garch, unconditional_variance, distribution_power_moment, effective_omega
  implicit none
  private
  public :: simulate_garch, simulate_conditional
contains

  function simulate_garch(spec,par,n,paths,burn,seed,vreg) result(out)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    integer,intent(in)::n
    integer,intent(in),optional::paths,burn
    integer(i8),intent(in),optional::seed
    real(dp),intent(in),optional::vreg(:,:)
    type(garch_simulation)::out
    integer::np,nb,m,total,j,t,start_keep
    integer(i8)::base_seed
    real(dp),allocatable::series(:),variance(:),residual(:),stdres(:),permanent(:),transitory(:),power_sigma(:),log_variance(:),z(:)
    type(rng_state)::rng
    real(dp)::uv,omega_eff,vterm,kappa
    logical::ok
    np=1
    if(present(paths))np=max(1,paths)
    nb=200
    if(present(burn))nb=max(0,burn)
    if(n<1)then
    out%message='simulation length must be positive'
    return
    end if
    if(present(vreg))then
      if(size(vreg,1)/=n+nb.or.size(vreg,2)/=size(par%xi))then
      out%message='variance regressor dimensions do not conform'
      return
      end if
    end if
    m=max(1,max(spec%p,spec%q))
    total=m+nb+n
    start_keep=m+nb+1
    allocate(out%series(n,np),out%sigma(n,np),out%innovations(n,np))
    allocate(out%permanent_component(n,np),out%transitory_component(n,np))
    out%series=0.0_dp
    out%sigma=0.0_dp
    out%innovations=0.0_dp
    out%permanent_component=0.0_dp
    out%transitory_component=0.0_dp
    base_seed=104729_i8
    if(present(seed))base_seed=seed
    uv=unconditional_variance(spec,par)
    if(.not.ieee_is_finite(uv).or.uv<=0.0_dp)uv=max(abs(par%omega),1.0_dp)
    omega_eff=par%omega
    kappa=distribution_power_moment(spec%distribution,par%dist,0.0_dp,1.0_dp,0.0_dp,1)
    do j=1,np
      call seed_rng(rng,base_seed+int(104729*j,i8))
      allocate(series(total),variance(total),residual(total),stdres(total), &
        permanent(total),transitory(total),power_sigma(total),log_variance(total))
      series=par%mu
      variance=uv
      residual=0.0_dp
      stdres=0.0_dp
      permanent=uv
      transitory=0.0_dp
      power_sigma=max(uv,1.0e-14_dp)**(par%delta/2.0_dp)
      log_variance=log(max(uv,1.0e-14_dp))
      z=rdist(spec%distribution,m,rng,par%dist)
      do t=1,m
        stdres(t)=z(t)
        residual(t)=sqrt(uv)*z(t)
        series(t)=par%mu+residual(t)
      end do
      do t=m+1,total
        vterm=0.0_dp
        if(present(vreg).and.size(par%xi)>0)vterm=dot_product(vreg(t-m,:),par%xi)
        call variance_step(t,spec,par,omega_eff,vterm,kappa,variance, &
          residual,stdres,permanent,transitory,power_sigma,log_variance,ok)
        if(.not.ok)then
        out%message='variance recursion became nonpositive or nonfinite'
        return
        end if
        z=rdist(spec%distribution,1,rng,par%dist)
        stdres(t)=z(1)
        residual(t)=sqrt(variance(t))*z(1)
        series(t)=par%mu+residual(t)
        if(t>=start_keep)then
          out%series(t-start_keep+1,j)=series(t)
          out%sigma(t-start_keep+1,j)=sqrt(variance(t))
          out%innovations(t-start_keep+1,j)=stdres(t)
          out%permanent_component(t-start_keep+1,j)=permanent(t)
          out%transitory_component(t-start_keep+1,j)=transitory(t)
        end if
      end do
      deallocate(series,variance,residual,stdres,permanent,transitory,power_sigma,log_variance)
    end do
    out%status=tsg_success
    out%message='ok'
  end function simulate_garch

  function simulate_conditional(history,spec,par,horizon,paths,seed,vreg_history,vreg_future) result(out)
    real(dp),intent(in)::history(:)
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    integer,intent(in)::horizon
    integer,intent(in),optional::paths
    integer(i8),intent(in),optional::seed
    real(dp),intent(in),optional::vreg_history(:,:),vreg_future(:,:)
    type(garch_simulation)::out
    integer::np,j,t,n0,total
    integer(i8)::base_seed
    real(dp),allocatable::series(:),variance(:),residual(:),stdres(:),permanent(:),transitory(:),power_sigma(:),log_variance(:),z(:)
    type(garch_filter_result)::filtered
    type(rng_state)::rng
    real(dp)::omega_eff,vterm,kappa
    logical::ok
    np=1000
    if(present(paths))np=max(1,paths)
    if(horizon<1.or.size(history)<max(2,max(spec%p,spec%q)+1))then
    out%message='invalid history or horizon'
    return
    end if
    if(present(vreg_history).neqv.present(vreg_future))then
    out%message='both historical and future regressors are required'
    return
    end if
    if(present(vreg_history))then
      if(size(vreg_history,1)/=size(history).or.size(vreg_future,1)/=horizon.or.&
         size(vreg_history,2)/=size(par%xi).or.size(vreg_future,2)/=size(par%xi))then
        out%message='variance regressor dimensions do not conform'
        return
      end if
      filtered=filter_garch(history,spec,par,vreg_history)
      omega_eff=effective_omega(history,spec,par,vreg_history)
    else
      filtered=filter_garch(history,spec,par)
      omega_eff=effective_omega(history,spec,par)
    end if
    if(filtered%status/=tsg_success)then
    out%message=filtered%message(1:len(out%message))
    return
    end if
    n0=size(history)
    total=n0+horizon
    kappa=distribution_power_moment(spec%distribution,par%dist,0.0_dp,1.0_dp,0.0_dp,1)
    allocate(out%series(horizon,np),out%sigma(horizon,np),out%innovations(horizon,np))
    allocate(out%permanent_component(horizon,np),out%transitory_component(horizon,np))
    out%series=0.0_dp
    out%sigma=0.0_dp
    out%innovations=0.0_dp
    out%permanent_component=0.0_dp
    out%transitory_component=0.0_dp
    base_seed=32452843_i8
    if(present(seed))base_seed=seed
    do j=1,np
      call seed_rng(rng,base_seed+int(49999*j,i8))
      allocate(series(total),variance(total),residual(total),stdres(total), &
        permanent(total),transitory(total),power_sigma(total),log_variance(total))
      series(1:n0)=history
      variance(1:n0)=filtered%variance
      residual(1:n0)=filtered%residuals
      stdres(1:n0)=filtered%standardized_residuals
      permanent(1:n0)=filtered%permanent_component
      transitory(1:n0)=filtered%transitory_component
      power_sigma(1:n0)=max(variance(1:n0),1.0e-14_dp)**(par%delta/2.0_dp)
      log_variance(1:n0)=log(max(variance(1:n0),1.0e-14_dp))
      do t=n0+1,total
        vterm=0.0_dp
        if(present(vreg_future).and.size(par%xi)>0)vterm=dot_product(vreg_future(t-n0,:),par%xi)
        call variance_step(t,spec,par,omega_eff,vterm,kappa,variance, &
          residual,stdres,permanent,transitory,power_sigma,log_variance,ok)
        if(.not.ok)then
        out%message='variance recursion became nonpositive or nonfinite'
        return
        end if
        z=rdist(spec%distribution,1,rng,par%dist)
        stdres(t)=z(1)
        residual(t)=sqrt(variance(t))*z(1)
        series(t)=par%mu+residual(t)
        out%series(t-n0,j)=series(t)
        out%sigma(t-n0,j)=sqrt(variance(t))
        out%innovations(t-n0,j)=stdres(t)
        out%permanent_component(t-n0,j)=permanent(t)
        out%transitory_component(t-n0,j)=transitory(t)
      end do
      deallocate(series,variance,residual,stdres,permanent,transitory,power_sigma,log_variance)
    end do
    out%status=tsg_success
    out%message='ok'
  end function simulate_conditional

  subroutine variance_step(t,spec,par,omega_eff,vterm,kappa,variance, &
    residual,stdres,permanent,transitory,power_sigma,log_variance,ok)
    integer,intent(in)::t
    type(garch_spec),intent(in)::spec
    type(garch_parameters),intent(in)::par
    real(dp),intent(in)::omega_eff,vterm,kappa
    real(dp),intent(inout)::variance(:),residual(:),stdres(:),permanent(:),transitory(:),power_sigma(:),log_variance(:)
    logical,intent(out)::ok
    integer::j
    real(dp)::base,powv,a,beta_value
    base=omega_eff+vterm
    if(spec%multiplicative)base=exp(base)
    select case(trim(spec%model))
    case('garch','igarch')
      variance(t)=base
      do j=1,min(spec%p,t-1)
      variance(t)=variance(t)+par%alpha(j)*residual(t-j)**2
      end do
      do j=1,min(spec%q,t-1)
      variance(t)=variance(t)+par%beta(j)*variance(t-j)
      end do
    case('ewma')
      beta_value=min(max(par%beta(1),1.0e-8_dp),0.999999_dp)
      variance(t)=(1.0_dp-beta_value)*residual(t-1)**2+beta_value*variance(t-1)
    case('gjrgarch')
      variance(t)=base
      do j=1,min(spec%p,t-1)
        variance(t)=variance(t)+par%alpha(j)*residual(t-j)**2
        if(residual(t-j)<=0.0_dp)variance(t)=variance(t)+par%gamma(j)*residual(t-j)**2
      end do
      do j=1,min(spec%q,t-1)
      variance(t)=variance(t)+par%beta(j)*variance(t-j)
      end do
    case('egarch')
      log_variance(t)=omega_eff+vterm
      do j=1,min(spec%p,t-1)
      log_variance(t)=log_variance(t)+par%alpha(j)*stdres(t-j)+par%gamma(j)*(abs(stdres(t-j))-kappa)
      end do
      do j=1,min(spec%q,t-1)
      log_variance(t)=log_variance(t)+par%beta(j)*log_variance(t-j)
      end do
      variance(t)=exp(min(max(log_variance(t),-700.0_dp),700.0_dp))
    case('aparch')
      power_sigma(t)=base
      do j=1,min(spec%p,t-1)
      power_sigma(t)=power_sigma(t)+par%alpha(j)*(abs(residual(t-j))-par%gamma(j)*residual(t-j))**par%delta
      end do
      do j=1,min(spec%q,t-1)
      power_sigma(t)=power_sigma(t)+par%beta(j)*power_sigma(t-j)
      end do
      variance(t)=max(power_sigma(t),1.0e-14_dp)**(2.0_dp/par%delta)
    case('fgarch')
      power_sigma(t)=base
      do j=1,min(spec%p,t-1)
        powv=max(variance(t-j),1.0e-14_dp)**(par%delta/2.0_dp)
        a=abs(stdres(t-j)-par%eta(j))-par%gamma(j)*(stdres(t-j)-par%eta(j))
        power_sigma(t)=power_sigma(t)+par%alpha(j)*powv*a**par%delta
      end do
      do j=1,min(spec%q,t-1)
      power_sigma(t)=power_sigma(t)+par%beta(j)*power_sigma(t-j)
      end do
      variance(t)=max(power_sigma(t),1.0e-14_dp)**(2.0_dp/par%delta)
    case('cgarch')
      permanent(t)=base+par%rho*permanent(t-1)+par%phi*(residual(t-1)**2-variance(t-1))
      transitory(t)=0.0_dp
      do j=1,min(spec%p,t-1)
      transitory(t)=transitory(t)+par%alpha(j)*transitory(t-j)+par%alpha(j)*(residual(t-j)**2-variance(t-j))
      end do
      do j=1,min(spec%q,t-1)
      transitory(t)=transitory(t)+par%beta(j)*transitory(t-j)
      end do
      variance(t)=permanent(t)+transitory(t)
    case default
      variance(t)=-1.0_dp
    end select
    ok=ieee_is_finite(variance(t)).and.variance(t)>1.0e-14_dp.and.variance(t)<1.0e100_dp
  end subroutine variance_step
end module tsgarch_simulation_module
