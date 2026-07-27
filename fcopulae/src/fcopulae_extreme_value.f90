! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae_extreme_value
  use fcopulae_kinds, only : dp, i8
  use fcopulae_rng, only : rng_state, seed_rng, uniform_rng
  use fcopulae_special, only : normal_cdf
  use fcopulae_utils, only : clamp01, numerical_first, numerical_second, kendall_tau_sample
  use fcopulae_integration, only : integration_result, integrate_1d, adapt_integrate2d
  use fcopulae_optimizer, only : optimizer_result, nelder_mead
  use fcopulae_archimedean, only : copula_fit_result
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: ev_param_count, ev_default_param, ev_bounds, ev_check
  public :: ev_dependence, ev_dependence_derivative, ev_cdf, ev_density, ev_rng
  public :: ev_tau, ev_rho, ev_tail_coeff, ev_fit
contains
  pure function ev_param_count(type_name) result(k)
    character(len=*),intent(in)::type_name
    integer::k
    select case(trim(adjustl(type_name)))
    case('tawn');k=3
    case('bb5');k=2
    case default;k=1
    end select
  end function ev_param_count

  subroutine ev_default_param(type_name,param)
    character(len=*),intent(in)::type_name
    real(dp),allocatable,intent(out)::param(:)
    select case(trim(adjustl(type_name)))
    case('gumbel','galambos','husler.reiss');allocate(param(1));param=2.0_dp
    case('tawn');allocate(param(3));param=[0.5_dp,0.5_dp,2.0_dp]
    case('bb5');allocate(param(2));param=[1.0_dp,2.0_dp]
    case default;allocate(param(1));param=1.0_dp
    end select
  end subroutine ev_default_param

  subroutine ev_bounds(type_name,lower,upper,bound)
    character(len=*),intent(in)::type_name
    real(dp),allocatable,intent(out)::lower(:),upper(:)
    real(dp),intent(in),optional::bound
    real(dp)::b
    b=20.0_dp;if(present(bound))b=bound
    select case(trim(adjustl(type_name)))
    case('gumbel');allocate(lower(1),upper(1));lower=1.0_dp;upper=b
    case('galambos','husler.reiss');allocate(lower(1),upper(1));lower=1.0e-4_dp;upper=b
    case('tawn');allocate(lower(3),upper(3));lower=[0.0_dp,0.0_dp,1.0_dp];upper=[1.0_dp,1.0_dp,b]
    case('bb5');allocate(lower(2),upper(2));lower=[1.0e-4_dp,1.0_dp];upper=[b,b]
    case default;allocate(lower(1),upper(1));lower=0.0_dp;upper=b
    end select
  end subroutine ev_bounds

  logical function ev_check(param,type_name)
    real(dp),intent(in)::param(:)
    character(len=*),intent(in)::type_name
    real(dp),allocatable::lo(:),hi(:)
    call ev_bounds(type_name,lo,hi,max(20.0_dp,maxval(abs(param))+1.0_dp))
    ev_check=size(param)==size(lo).and.all(param>=lo).and.all(param<=hi)
  end function ev_check

  function ev_dependence(x,param,type_name) result(a)
    real(dp),intent(in)::x,param(:)
    character(len=*),intent(in)::type_name
    real(dp)::a,z,alpha,beta,r,delta,theta,p,q
    z=max(1.0e-12_dp,min(1.0_dp-1.0e-12_dp,x))
    select case(trim(adjustl(type_name)))
    case('gumbel')
      alpha=param(1);if(abs(alpha-1.0_dp)<1e-12_dp)then;a=1.0_dp;else;a=(z**alpha+(1.0_dp-z)**alpha)**(1.0_dp/alpha);end if
    case('galambos')
      alpha=param(1);if(alpha<=1.0e-10_dp)then;a=1.0_dp;else;a=1.0_dp-(z**(-alpha)+(1.0_dp-z)**(-alpha))**(-1.0_dp/alpha);end if
    case('husler.reiss')
      alpha=param(1);if(alpha<=1.0e-10_dp)then;a=1.0_dp;else
        p=1.0_dp/alpha+0.5_dp*alpha*log(z/(1.0_dp-z));q=1.0_dp/alpha-0.5_dp*alpha*log(z/(1.0_dp-z))
        a=z*normal_cdf(p)+(1.0_dp-z)*normal_cdf(q)
      end if
    case('tawn')
      alpha=param(1);beta=param(2);r=param(3)
      if(alpha<=1.0e-12_dp.or.beta<=1.0e-12_dp.or.abs(r-1.0_dp)<1.0e-12_dp)then;a=1.0_dp
      else;a=1.0_dp-beta+(beta-alpha)*z+((alpha*z)**r+(beta*(1.0_dp-z))**r)**(1.0_dp/r);end if
    case('bb5')
      delta=param(1);theta=param(2)
      if(abs(theta-1.0_dp)<1.0e-12_dp)then;a=1.0_dp-(z**(-delta)+(1.0_dp-z)**(-delta))**(-1.0_dp/delta)
      else;a=(z**theta+(1.0_dp-z)**theta-(z**(-delta*theta)+(1.0_dp-z)**(-delta*theta))**(-1.0_dp/delta))**(1.0_dp/theta);end if
    case('gumbelII');alpha=param(1);a=alpha*z*z-alpha*z+1.0_dp
    case('pi','Cperp');a=1.0_dp
    case('m','Cplus');a=max(z,1.0_dp-z)
    case default;a=1.0_dp
    end select
    a=max(max(z,1.0_dp-z),min(1.0_dp,a))
  end function ev_dependence

  function ev_dependence_derivative(x,param,type_name,order) result(d)
    real(dp),intent(in)::x,param(:)
    character(len=*),intent(in)::type_name
    integer,intent(in)::order
    real(dp)::d
    if(order==1)then;d=numerical_first(local,x,1.0e-8_dp,1.0_dp-1.0e-8_dp)
    else;d=numerical_second(local,x,1.0e-8_dp,1.0_dp-1.0e-8_dp);end if
  contains
    function local(z) result(v);real(dp),intent(in)::z;real(dp)::v;v=ev_dependence(z,param,type_name);end function local
  end function ev_dependence_derivative

  function ev_cdf(u,v,param,type_name) result(c)
    real(dp),intent(in)::u,v,param(:)
    character(len=*),intent(in)::type_name
    real(dp)::c,lu,lv,s,x,a
    if(u<=0.0_dp.or.v<=0.0_dp)then;c=0.0_dp;return;end if
    if(u>=1.0_dp)then;c=clamp01(v);return;end if
    if(v>=1.0_dp)then;c=clamp01(u);return;end if
    lu=log(u);lv=log(v);s=lu+lv;x=lu/s;a=ev_dependence(x,param,type_name)
    c=exp(s*a);c=max(max(0.0_dp,u+v-1.0_dp),min(min(u,v),c))
  end function ev_cdf

  function ev_density(u,v,param,type_name) result(c)
    real(dp),intent(in)::u,v,param(:)
    character(len=*),intent(in)::type_name
    real(dp)::c,lu,lv,s,x,y,a,a1,a2,pref
    if(u<=0.0_dp.or.u>=1.0_dp.or.v<=0.0_dp.or.v>=1.0_dp)then;c=0.0_dp;return;end if
    lu=log(u);lv=log(v);s=lu+lv;x=lu/s;y=lv/s
    a=ev_dependence(x,param,type_name)
    a1=ev_dependence_derivative(x,param,type_name,1)
    a2=ev_dependence_derivative(x,param,type_name,2)
    pref=ev_cdf(u,v,param,type_name)/(u*v)
    c=pref*((-x*y/s)*a2+(a+y*a1)*(a-x*a1));c=max(0.0_dp,c)
    if(.not.ieee_is_finite(c))c=0.0_dp
  end function ev_density

  function conditional_v_cdf(v,u,param,type_name) result(p)
    real(dp),intent(in)::v,u,param(:)
    character(len=*),intent(in)::type_name
    real(dp)::p,h,um,up
    h=2.0e-5_dp;um=max(1.0e-9_dp,u-h);up=min(1.0_dp-1.0e-9_dp,u+h)
    p=(ev_cdf(up,v,param,type_name)-ev_cdf(um,v,param,type_name))/(up-um)
    p=max(0.0_dp,min(1.0_dp,p))
  end function conditional_v_cdf

  subroutine ev_rng(n,param,type_name,x,seed)
    integer,intent(in)::n
    real(dp),intent(in)::param(:)
    character(len=*),intent(in)::type_name
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    type(rng_state)::state
    integer(i8)::s0
    integer::i,j
    real(dp)::u,w,lo,hi,mid
    s0=246813579_i8;if(present(seed))s0=seed;call seed_rng(state,s0);allocate(x(n,2))
    do i=1,n
      u=uniform_rng(state);w=uniform_rng(state);lo=1.0e-10_dp;hi=1.0_dp-1.0e-10_dp
      do j=1,80
        mid=0.5_dp*(lo+hi)
        if(conditional_v_cdf(mid,u,param,type_name)<w)then;lo=mid;else;hi=mid;end if
      end do
      x(i,:)=[u,0.5_dp*(lo+hi)]
    end do
  end subroutine ev_rng

  function ev_tau(param,type_name) result(tau)
    real(dp),intent(in)::param(:)
    character(len=*),intent(in)::type_name
    real(dp)::tau
    type(integration_result)::res
    res=integrate_1d(fun,1.0e-5_dp,1.0_dp-1.0e-5_dp,1.0e-6_dp,20);tau=max(-1.0_dp,min(1.0_dp,res%value))
  contains
    function fun(x) result(z)
      real(dp),intent(in)::x;real(dp)::z,a
      a=ev_dependence(x,param,type_name);z=x*(1.0_dp-x)*ev_dependence_derivative(x,param,type_name,2)/a
    end function fun
  end function ev_tau

  function ev_rho(param,type_name) result(rho)
    real(dp),intent(in)::param(:)
    character(len=*),intent(in)::type_name
    real(dp)::rho
    type(integration_result)::res
    res=integrate_1d(fun,0.0_dp,1.0_dp,1.0e-8_dp,20);rho=max(-1.0_dp,min(1.0_dp,res%value))
  contains
    function fun(x) result(z)
      real(dp),intent(in)::x;real(dp)::z,a
      a=ev_dependence(x,param,type_name);z=12.0_dp/(a+1.0_dp)**2-3.0_dp
    end function fun
  end function ev_rho

  subroutine ev_tail_coeff(param,type_name,lower,upper)
    real(dp),intent(in)::param(:)
    character(len=*),intent(in)::type_name
    real(dp),intent(out)::lower,upper
    lower=0.0_dp;upper=2.0_dp-2.0_dp*ev_dependence(0.5_dp,param,type_name)
    upper=max(0.0_dp,min(1.0_dp,upper))
  end subroutine ev_tail_coeff

  function ev_fit(u,v,type_name,start,bound,max_iter) result(fit)
    real(dp),intent(in)::u(:),v(:)
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::start(:),bound
    integer,intent(in),optional::max_iter
    type(copula_fit_result)::fit
    type(optimizer_result)::opt
    real(dp),allocatable::x0(:),lo(:),hi(:),step(:)
    integer::k,n,niter
    niter=3500;if(present(max_iter))niter=max(1,max_iter)
    call ev_default_param(type_name,x0);call ev_bounds(type_name,lo,hi,bound)
    if(present(start))then;if(size(start)==size(x0))x0=start;end if
    k=size(x0);allocate(step(k));step=0.1_dp*max(1.0_dp,abs(x0));x0=max(lo,min(hi,x0))
    opt=nelder_mead(obj,x0,step,1.0e-7_dp,niter,lo,hi)
    allocate(fit%param(k))
    fit%param=opt%x
    fit%loglik=-opt%value
    n=min(size(u),size(v))
    fit%aic=2.0_dp*real(k,dp)-2.0_dp*fit%loglik
    fit%bic=real(k,dp)*log(real(max(1,n),dp))-2.0_dp*fit%loglik
    fit%iterations=opt%iterations;fit%evaluations=opt%evaluations;fit%converged=opt%converged
  contains
    function obj(p) result(val)
      real(dp),intent(in)::p(:);real(dp)::val,d
      integer::i,m
      val=0.0_dp;m=min(size(u),size(v))
      do i=1,m
        d=ev_density(clamp01(u(i)),clamp01(v(i)),p,type_name)
        if(d<=1.0e-300_dp.or..not.ieee_is_finite(d))then;val=val+690.0_dp;else;val=val-log(d);end if
      end do
    end function obj
  end function ev_fit
end module fcopulae_extreme_value
