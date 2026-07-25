! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae_archimedean
  use fcopulae_kinds, only : dp, i8, pi
  use fcopulae_rng, only : rng_state, seed_rng, uniform_rng
  use fcopulae_utils, only : clamp01, numerical_first, numerical_second, safe_log
  use fcopulae_integration, only : integration_result, integrate_1d, adapt_integrate2d
  use fcopulae_optimizer, only : optimizer_result, nelder_mead
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: archm_default_alpha, archm_range, archm_check
  public :: archm_phi, archm_inv_phi, archm_phi0, archm_phi_derivative
  public :: archm_inv_phi_derivative, archm_k, archm_inv_k
  public :: archm_cdf, archm_density, archm_rng
  public :: archm_tau, archm_rho, archm_tail_coeff, archm_fit, copula_fit_result

  type :: copula_fit_result
    real(dp), allocatable :: param(:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    real(dp) :: bic = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
  end type copula_fit_result
contains
  pure function archm_default_alpha(type_id) result(a)
    integer,intent(in)::type_id
    real(dp)::a
    real(dp),parameter::vals(22)=[1.0_dp,2.0_dp,0.5_dp,2.0_dp,1.0_dp,2.0_dp,0.5_dp,2.0_dp, &
      0.5_dp,0.5_dp,0.2_dp,2.0_dp,1.0_dp,2.0_dp,2.0_dp,1.0_dp,0.5_dp,3.0_dp,1.0_dp,1.0_dp,2.0_dp,0.5_dp]
    if(type_id>=1.and.type_id<=22)then;a=vals(type_id);else;a=1.0_dp;end if
  end function archm_default_alpha

  pure subroutine archm_range(type_id,lower,upper,bound)
    integer,intent(in)::type_id
    real(dp),intent(out)::lower,upper
    real(dp),intent(in),optional::bound
    real(dp)::b
    b=20.0_dp;if(present(bound))b=bound
    select case(type_id)
    case(1);lower=-1.0_dp;upper=b
    case(2);lower=1.0_dp;upper=b
    case(3);lower=-1.0_dp;upper=1.0_dp
    case(4);lower=1.0_dp;upper=b
    case(5);lower=-b;upper=b
    case(6);lower=1.0_dp;upper=b
    case(7);lower=0.0_dp;upper=1.0_dp
    case(8);lower=1.0_dp;upper=b
    case(9,10);lower=0.0_dp;upper=1.0_dp
    case(11);lower=0.0_dp;upper=0.5_dp
    case(12);lower=1.0_dp;upper=b
    case(13);lower=0.0_dp;upper=b
    case(14,15);lower=1.0_dp;upper=b
    case(16);lower=0.0_dp;upper=b
    case(17);lower=-b;upper=b
    case(18);lower=2.0_dp;upper=b
    case(19,20);lower=0.0_dp;upper=b
    case(21);lower=1.0_dp;upper=b
    case(22);lower=0.0_dp;upper=1.0_dp
    case default;lower=0.0_dp;upper=0.0_dp
    end select
  end subroutine archm_range

  pure logical function archm_check(alpha,type_id)
    real(dp),intent(in)::alpha
    integer,intent(in)::type_id
    real(dp)::lo,hi
    call archm_range(type_id,lo,hi,max(20.0_dp,abs(alpha)+1.0_dp))
    archm_check=type_id>=1.and.type_id<=22.and.alpha>=lo.and.alpha<=hi
  end function archm_check

  pure function archm_phi0(alpha,type_id) result(p0)
    real(dp),intent(in)::alpha
    integer,intent(in)::type_id
    real(dp)::p0
    select case(type_id)
    case(1);if(alpha<0.0_dp)then;p0=-1.0_dp/alpha;else;p0=huge(1.0_dp);end if
    case(2);p0=1.0_dp
    case(3:6);p0=huge(1.0_dp)
    case(7);if(abs(alpha)<1e-14_dp)then;p0=1.0_dp;else;p0=-log(1.0_dp-alpha);end if
    case(8);p0=1.0_dp
    case(9,10);p0=huge(1.0_dp)
    case(11);if(abs(alpha)<1e-14_dp)then;p0=huge(1.0_dp);else;p0=log(2.0_dp);end if
    case(12:14);p0=huge(1.0_dp)
    case(15);p0=1.0_dp
    case(16);if(abs(alpha)<1e-14_dp)then;p0=1.0_dp;else;p0=huge(1.0_dp);end if
    case(17);p0=huge(1.0_dp)
    case(18);p0=exp(-alpha)
    case(19,20);p0=huge(1.0_dp)
    case(21);p0=1.0_dp
    case(22);if(abs(alpha)<1e-14_dp)then;p0=huge(1.0_dp);else;p0=0.5_dp*pi;end if
    case default;p0=huge(1.0_dp)
    end select
  end function archm_phi0

  pure function archm_phi(x,alpha,type_id) result(f)
    real(dp),intent(in)::x,alpha
    integer,intent(in)::type_id
    real(dp)::f,z
    z=clamp01(x,1.0e-14_dp)
    select case(type_id)
    case(1)
      if(abs(alpha+1.0_dp)<1e-12_dp)then;f=1.0_dp-z
      else if(abs(alpha)<1e-12_dp)then;f=-log(z)
      else if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp/z-1.0_dp
      else;f=(z**(-alpha)-1.0_dp)/alpha;end if
    case(2)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp-z;else;f=(1.0_dp-z)**alpha;end if
    case(3)
      if(abs(alpha)<1e-12_dp)then;f=-log(z)
      else if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp/z-1.0_dp
      else;f=log((1.0_dp-alpha*(1.0_dp-z))/z);end if
    case(4)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=-log(z);else;f=(-log(z))**alpha;end if
    case(5)
      if(abs(alpha)<1e-10_dp)then;f=-log(z)
      else;f=-log((exp(-alpha*z)-1.0_dp)/(exp(-alpha)-1.0_dp));end if
    case(6)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=-log(z);else;f=-log(1.0_dp-(1.0_dp-z)**alpha);end if
    case(7)
      if(abs(alpha)<1e-12_dp)then;f=1.0_dp-z
      else if(abs(alpha-1.0_dp)<1e-12_dp)then;f=-log(z)
      else;f=-log(alpha*z+1.0_dp-alpha);end if
    case(8)
      if(abs(alpha)<1e-12_dp)then;f=-log(z);else;f=(1.0_dp-z)/(1.0_dp+z*(alpha-1.0_dp));end if
    case(9)
      if(abs(alpha)<1e-12_dp)then;f=-log(z);else;f=log(1.0_dp-alpha*log(z));end if
    case(10)
      if(abs(alpha)<1e-12_dp)then;f=-log(z);else;f=log(2.0_dp*z**(-alpha)-1.0_dp);end if
    case(11)
      if(abs(alpha)<1e-12_dp)then;f=-log(z);else;f=log(2.0_dp-z**alpha);end if
    case(12)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp/z-1.0_dp;else;f=(1.0_dp/z-1.0_dp)**alpha;end if
    case(13)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=-log(z);else;f=(1.0_dp-log(z))**alpha-1.0_dp;end if
    case(14)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp/z-1.0_dp;else;f=(z**(-1.0_dp/alpha)-1.0_dp)**alpha;end if
    case(15)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp-z;else;f=(1.0_dp-z**(1.0_dp/alpha))**alpha;end if
    case(16)
      if(abs(alpha)<1e-12_dp)then;f=1.0_dp-z;else;f=(alpha/z+1.0_dp)*(1.0_dp-z);end if
    case(17)
      if(abs(alpha+1.0_dp)<1e-12_dp)then;f=-log(z)
      else;f=-log(((1.0_dp+z)**(-alpha)-1.0_dp)/(2.0_dp**(-alpha)-1.0_dp));end if
    case(18);f=exp(alpha/(z-1.0_dp))
    case(19)
      if(abs(alpha)<1e-12_dp)then;f=1.0_dp/z-1.0_dp;else;f=exp(alpha/z)-exp(alpha);end if
    case(20)
      if(abs(alpha)<1e-12_dp)then;f=-log(z);else;f=exp(z**(-alpha))-exp(1.0_dp);end if
    case(21)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;f=1.0_dp-z;else;f=1.0_dp-(1.0_dp-(1.0_dp-z)**alpha)**(1.0_dp/alpha);end if
    case(22)
      if(abs(alpha)<1e-12_dp)then;f=-log(z);else;f=asin(max(-1.0_dp,min(1.0_dp,1.0_dp-z**alpha)));end if
    case default;f=huge(1.0_dp)
    end select
    if(x<=0.0_dp)f=archm_phi0(alpha,type_id)
    if(x>=1.0_dp)f=0.0_dp
  end function archm_phi

  pure function archm_inv_phi(x,alpha,type_id) result(u)
    real(dp),intent(in)::x,alpha
    integer,intent(in)::type_id
    real(dp)::u,z,p0
    z=max(0.0_dp,x);p0=archm_phi0(alpha,type_id)
    if(z<=0.0_dp)then;u=1.0_dp;return;end if
    if(p0<huge(1.0_dp)/2.0_dp .and. z>=p0)then;u=0.0_dp;return;end if
    select case(type_id)
    case(1)
      if(abs(alpha+1.0_dp)<1e-12_dp)then;u=1.0_dp-z
      else if(abs(alpha)<1e-12_dp)then;u=exp(-z)
      else if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp/(1.0_dp+z)
      else;u=exp(-log(max(1.0e-300_dp,1.0_dp+alpha*z))/alpha);end if
    case(2)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp-z;else;u=1.0_dp-z**(1.0_dp/alpha);end if
    case(3)
      if(abs(alpha)<1e-12_dp)then;u=exp(-z)
      else if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp/(1.0_dp+z)
      else;u=(1.0_dp-alpha)/(exp(z)-alpha);end if
    case(4)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=exp(-z);else;u=exp(-z**(1.0_dp/alpha));end if
    case(5)
      if(abs(alpha)<1e-10_dp)then;u=exp(-z);else;u=-log(1.0_dp+exp(-z)*(exp(-alpha)-1.0_dp))/alpha;end if
    case(6)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=exp(-z);else;u=1.0_dp-(1.0_dp-exp(-z))**(1.0_dp/alpha);end if
    case(7)
      if(abs(alpha)<1e-12_dp)then;u=1.0_dp-z
      else if(abs(alpha-1.0_dp)<1e-12_dp)then;u=exp(-z)
      else;u=(1.0_dp-exp(z)+alpha*exp(z))/(alpha*exp(z));end if
    case(8)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp-z;else;u=(1.0_dp-z)/((alpha-1.0_dp)*z+1.0_dp);end if
    case(9)
      if(abs(alpha)<1e-12_dp)then;u=exp(-z);else;u=exp((1.0_dp-exp(z))/alpha);end if
    case(10)
      if(abs(alpha)<1e-12_dp)then;u=exp(-z);else;u=((1.0_dp+exp(z))/2.0_dp)**(-1.0_dp/alpha);end if
    case(11)
      if(abs(alpha)<1e-12_dp)then;u=exp(-z);else;u=max(0.0_dp,2.0_dp-exp(z))**(1.0_dp/alpha);end if
    case(12)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp/(1.0_dp+z);else;u=1.0_dp/(1.0_dp+z**(1.0_dp/alpha));end if
    case(13)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=exp(-z);else;u=exp(1.0_dp-(1.0_dp+z)**(1.0_dp/alpha));end if
    case(14)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp/(1.0_dp+z);else;u=(1.0_dp+z**(1.0_dp/alpha))**(-alpha);end if
    case(15)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp-z;else;u=max(0.0_dp,1.0_dp-z**(1.0_dp/alpha))**alpha;end if
    case(16)
      if(abs(alpha)<1e-12_dp)then;u=1.0_dp-z;else;u=(1.0_dp-alpha-z)/2.0_dp+sqrt((1.0_dp-alpha-z)**2/4.0_dp+alpha);end if
    case(17)
      if(abs(alpha+1.0_dp)<1e-12_dp)then;u=exp(-z)
      else;u=(exp(-z)*(2.0_dp**(-alpha)-1.0_dp)+1.0_dp)**(-1.0_dp/alpha)-1.0_dp;end if
    case(18)
      if(z<=0.0_dp)then;u=1.0_dp;else;u=1.0_dp+alpha/log(z);end if
    case(19)
      if(abs(alpha)<1e-12_dp)then;u=1.0_dp/(1.0_dp+z);else;u=alpha/log(z+exp(alpha));end if
    case(20)
      if(abs(alpha)<1e-12_dp)then;u=exp(-z);else;u=exp(-log(log(z+exp(1.0_dp)))/alpha);end if
    case(21)
      if(abs(alpha-1.0_dp)<1e-12_dp)then;u=1.0_dp-z;else;u=1.0_dp-(1.0_dp-(1.0_dp-z)**alpha)**(1.0_dp/alpha);end if
    case(22)
      if(abs(alpha)<1e-12_dp)then;u=exp(-z);else;u=max(0.0_dp,1.0_dp-sin(z))**(1.0_dp/alpha);end if
    case default;u=0.0_dp
    end select
    u=min(1.0_dp,max(0.0_dp,u))
  end function archm_inv_phi

  function archm_phi_derivative(x,alpha,type_id,order) result(d)
    real(dp),intent(in)::x,alpha
    integer,intent(in)::type_id,order
    real(dp)::d
    if(order==1)then;d=numerical_first(local,x,1.0e-10_dp,1.0_dp-1.0e-10_dp)
    else;d=numerical_second(local,x,1.0e-10_dp,1.0_dp-1.0e-10_dp);end if
  contains
    function local(z) result(f);real(dp),intent(in)::z;real(dp)::f;f=archm_phi(z,alpha,type_id);end function local
  end function archm_phi_derivative

  function archm_inv_phi_derivative(x,alpha,type_id,order) result(d)
    real(dp),intent(in)::x,alpha
    integer,intent(in)::type_id,order
    real(dp)::d,hi
    hi=archm_phi0(alpha,type_id);if(hi>1.0e6_dp)hi=max(100.0_dp,2.0_dp*x+10.0_dp)
    if(order==1)then;d=numerical_first(local,x,0.0_dp,hi)
    else;d=numerical_second(local,x,0.0_dp,hi);end if
  contains
    function local(z) result(f);real(dp),intent(in)::z;real(dp)::f;f=archm_inv_phi(z,alpha,type_id);end function local
  end function archm_inv_phi_derivative

  function archm_cdf(u,v,alpha,type_id) result(c)
    real(dp),intent(in)::u,v,alpha
    integer,intent(in)::type_id
    real(dp)::c,s
    if(u<=0.0_dp.or.v<=0.0_dp)then;c=0.0_dp;return;end if
    if(u>=1.0_dp)then;c=clamp01(v);return;end if
    if(v>=1.0_dp)then;c=clamp01(u);return;end if
    s=archm_phi(u,alpha,type_id)+archm_phi(v,alpha,type_id)
    c=archm_inv_phi(s,alpha,type_id)
    c=max(max(0.0_dp,u+v-1.0_dp),min(min(u,v),c))
  end function archm_cdf

  function archm_density(u,v,alpha,type_id) result(c)
    real(dp),intent(in)::u,v,alpha
    integer,intent(in)::type_id
    real(dp)::c,s,p1,p2,q2
    if(u<=0.0_dp.or.u>=1.0_dp.or.v<=0.0_dp.or.v>=1.0_dp)then;c=0.0_dp;return;end if
    s=archm_phi(u,alpha,type_id)+archm_phi(v,alpha,type_id)
    p1=archm_phi_derivative(u,alpha,type_id,1);p2=archm_phi_derivative(v,alpha,type_id,1)
    q2=archm_inv_phi_derivative(s,alpha,type_id,2)
    c=max(0.0_dp,q2*p1*p2)
    if(.not.ieee_is_finite(c))c=0.0_dp
  end function archm_density

  function archm_k(t,alpha,type_id) result(k)
    real(dp),intent(in)::t,alpha
    integer,intent(in)::type_id
    real(dp)::k,d
    if(t<=0.0_dp)then;k=0.0_dp;return;end if
    if(t>=1.0_dp)then;k=1.0_dp;return;end if
    d=archm_phi_derivative(t,alpha,type_id,1)
    if(abs(d)<1.0e-14_dp)then;k=t;else;k=t-archm_phi(t,alpha,type_id)/d;end if
    k=min(1.0_dp,max(0.0_dp,k))
  end function archm_k

  function archm_inv_k(p,alpha,type_id) result(t)
    real(dp),intent(in)::p,alpha
    integer,intent(in)::type_id
    real(dp)::t,lo,hi,mid
    integer::i
    if(p<=0.0_dp)then;t=0.0_dp;return;end if
    if(p>=1.0_dp)then;t=1.0_dp;return;end if
    lo=1.0e-10_dp;hi=1.0_dp-1.0e-10_dp
    do i=1,100
      mid=0.5_dp*(lo+hi)
      if(archm_k(mid,alpha,type_id)<p)then;lo=mid;else;hi=mid;end if
    end do
    t=0.5_dp*(lo+hi)
  end function archm_inv_k

  subroutine archm_rng(n,alpha,type_id,x,seed)
    integer,intent(in)::n,type_id
    real(dp),intent(in)::alpha
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    type(rng_state)::state
    integer::i
    integer(i8)::s
    real(dp)::a,b,t,p
    s=190734863_i8;if(present(seed))s=seed;call seed_rng(state,s);allocate(x(n,2))
    do i=1,n
      a=uniform_rng(state);b=uniform_rng(state);t=archm_inv_k(b,alpha,type_id);p=archm_phi(t,alpha,type_id)
      x(i,1)=archm_inv_phi(a*p,alpha,type_id);x(i,2)=archm_inv_phi((1.0_dp-a)*p,alpha,type_id)
    end do
  end subroutine archm_rng

  function archm_tau(alpha,type_id) result(tau)
    real(dp),intent(in)::alpha
    integer,intent(in)::type_id
    real(dp)::tau
    type(integration_result)::res
    res=integrate_1d(fun,1.0e-6_dp,1.0_dp-1.0e-6_dp,1.0e-7_dp)
    tau=1.0_dp+4.0_dp*res%value
    tau=max(-1.0_dp,min(1.0_dp,tau))
  contains
    function fun(t) result(v)
      real(dp),intent(in)::t;real(dp)::v,d
      d=archm_phi_derivative(t,alpha,type_id,1)
      if(abs(d)<1e-14_dp)then;v=0.0_dp;else;v=archm_phi(t,alpha,type_id)/d;end if
    end function fun
  end function archm_tau

  function archm_rho(alpha,type_id) result(rho)
    real(dp),intent(in)::alpha
    integer,intent(in)::type_id
    real(dp)::rho
    type(integration_result)::res
    res=adapt_integrate2d(fun,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],2.0e-5_dp,8)
    rho=12.0_dp*res%value-3.0_dp
    rho=max(-1.0_dp,min(1.0_dp,rho))
  contains
    function fun(u,v) result(z)
      real(dp),intent(in)::u,v;real(dp)::z
      z=archm_cdf(u,v,alpha,type_id)
    end function fun
  end function archm_rho

  subroutine archm_tail_coeff(alpha,type_id,lower,upper)
    real(dp),intent(in)::alpha
    integer,intent(in)::type_id
    real(dp),intent(out)::lower,upper
    real(dp),parameter::e=1.0e-7_dp
    lower=archm_cdf(e,e,alpha,type_id)/e
    upper=(1.0_dp-2.0_dp*(1.0_dp-e)+archm_cdf(1.0_dp-e,1.0_dp-e,alpha,type_id))/e
    lower=max(0.0_dp,min(1.0_dp,lower));upper=max(0.0_dp,min(1.0_dp,upper))
  end subroutine archm_tail_coeff

  function archm_fit(u,v,type_id,start,bound,max_iter) result(fit)
    real(dp),intent(in)::u(:),v(:)
    integer,intent(in)::type_id
    real(dp),intent(in),optional::start,bound
    integer,intent(in),optional::max_iter
    type(copula_fit_result)::fit
    type(optimizer_result)::opt
    real(dp)::x0(1),lo(1),hi(1),st(1),b
    integer::n,niter
    niter=2500;if(present(max_iter))niter=max(1,max_iter)
    b=20.0_dp;if(present(bound))b=bound;call archm_range(type_id,lo(1),hi(1),b)
    x0(1)=archm_default_alpha(type_id);if(present(start))x0(1)=start
    x0=max(lo,min(hi,x0));st=max(0.02_dp,0.1_dp*(hi-lo))
    opt=nelder_mead(obj,x0,st,1.0e-8_dp,niter,lo,hi)
    allocate(fit%param(1));fit%param=opt%x;fit%loglik=-opt%value
    n=min(size(u),size(v));fit%aic=2.0_dp-2.0_dp*fit%loglik;fit%bic=log(real(max(1,n),dp))-2.0_dp*fit%loglik
    fit%iterations=opt%iterations;fit%evaluations=opt%evaluations;fit%converged=opt%converged
  contains
    function obj(x) result(val)
      real(dp),intent(in)::x(:);real(dp)::val,d
      integer::i,m
      val=0.0_dp;m=min(size(u),size(v))
      do i=1,m
        d=archm_density(clamp01(u(i)),clamp01(v(i)),x(1),type_id)
        if(d<=1.0e-300_dp.or..not.ieee_is_finite(d))then;val=val+690.0_dp;else;val=val-log(d);end if
      end do
    end function obj
  end function archm_fit
end module fcopulae_archimedean
