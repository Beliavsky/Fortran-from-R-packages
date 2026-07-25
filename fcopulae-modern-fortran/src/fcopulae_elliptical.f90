! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fCopulae authors and translation contributors.
! This file is part of fcopulae-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fcopulae_elliptical
  use fcopulae_kinds, only : dp, i8, pi
  use fcopulae_rng, only : rng_state, seed_rng, uniform_rng, gamma_rng
  use fcopulae_special, only : normal_pdf, normal_cdf, normal_quantile, student_t_pdf, student_t_cdf, student_t_quantile
  use fcopulae_distributions, only : dnorm2d, pnorm2d, dt2d, pt2d
  use fcopulae_integration, only : integration_result, integrate_1d, adapt_integrate2d
  use fcopulae_optimizer, only : optimizer_result, nelder_mead
  use fcopulae_utils, only : clamp01, kendall_tau_sample
  use fcopulae_archimedean, only : copula_fit_result
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: elliptical_generator, elliptical_lambda, elliptical_marginal_pdf
  public :: elliptical_marginal_cdf, elliptical_marginal_quantile
  public :: elliptical_joint_pdf, elliptical_copula_density, elliptical_copula_cdf
  public :: elliptical_rng, elliptical_tau, elliptical_rho, elliptical_tail_coeff
  public :: elliptical_fit
contains
  pure function canonical_type(type_name) result(t)
    character(len=*),intent(in)::type_name
    character(len=16)::t
    t=trim(adjustl(type_name))
    if(t=='normal')t='norm'
    if(t=='student')t='t'
  end function canonical_type

  elemental function elliptical_generator(q,type_name,param1,param2) result(g)
    real(dp),intent(in)::q
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::g,nu,r,s,z
    character(len=16)::t
    t=canonical_type(type_name);z=max(0.0_dp,q)
    select case(trim(t))
    case('norm');g=exp(-0.5_dp*z)
    case('cauchy');g=(1.0_dp+z)**(-1.5_dp)
    case('t');nu=4.0_dp;if(present(param1))nu=param1;g=(1.0_dp+z/nu)**(-0.5_dp*(nu+2.0_dp))
    case('logistic');g=exp(-0.5_dp*z)/(1.0_dp+exp(-0.5_dp*z))**2
    case('laplace');g=exp(-sqrt(z))
    case('kotz');r=1.0_dp;if(present(param1))r=param1;g=exp(-0.5_dp*r*z)
    case('epower');r=1.0_dp;s=1.0_dp;if(present(param1))r=param1;if(present(param2))s=param2;g=exp(-r*(0.5_dp*z)**s)
    case default;g=0.0_dp
    end select
  end function elliptical_generator

  pure function elliptical_lambda(type_name,param1,param2) result(lambda)
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::lambda,r,s
    character(len=16)::t
    t=canonical_type(type_name)
    select case(trim(t))
    case('norm','cauchy','t');lambda=1.0_dp/(2.0_dp*pi)
    case('logistic');lambda=1.0_dp/pi
    case('laplace');lambda=1.0_dp/(2.0_dp*pi)
    case('kotz');r=1.0_dp;if(present(param1))r=param1;lambda=r/(2.0_dp*pi)
    case('epower')
      r=1.0_dp; s=1.0_dp
      if(present(param1)) r=param1
      if(present(param2)) s=param2
      lambda=s*r**(1.0_dp/s)/(2.0_dp*pi*gamma(1.0_dp/s))
    case default;lambda=0.0_dp
    end select
  end function elliptical_lambda

  function elliptical_marginal_pdf(x,type_name,param1,param2) result(f)
    real(dp),intent(in)::x
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::f,lambda,nu,r,ymax,h,y,sumv
    integer::i,n
    character(len=16)::t
    t=canonical_type(type_name)
    select case(trim(t))
    case('norm');f=normal_pdf(x);return
    case('cauchy');f=student_t_pdf(x,1.0_dp);return
    case('t');nu=4.0_dp;if(present(param1))nu=param1;f=student_t_pdf(x,nu);return
    case('kotz')
      r=1.0_dp;if(present(param1))r=param1;f=sqrt(r)*normal_pdf(sqrt(r)*x);return
    end select
    lambda=elliptical_lambda(t,param1,param2)
    ymax=integration_radius(t,param1,param2)
    n=400;h=ymax/real(n,dp);sumv=0.0_dp
    do i=0,n
      y=real(i,dp)*h
      if(i==0.or.i==n)then
        sumv=sumv+elliptical_generator(x*x+y*y,t,param1,param2)
      else if(mod(i,2)==0)then
        sumv=sumv+2.0_dp*elliptical_generator(x*x+y*y,t,param1,param2)
      else
        sumv=sumv+4.0_dp*elliptical_generator(x*x+y*y,t,param1,param2)
      end if
    end do
    f=max(0.0_dp,2.0_dp*lambda*h*sumv/3.0_dp)
  end function elliptical_marginal_pdf

  function elliptical_marginal_cdf(x,type_name,param1,param2) result(p)
    real(dp),intent(in)::x
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::p,lambda,q,nu,rmax,h,r,sumv,ang
    integer::i,n
    character(len=16)::t
    t=canonical_type(type_name)
    select case(trim(t))
    case('norm');p=normal_cdf(x);return
    case('cauchy');p=student_t_cdf(x,1.0_dp);return
    case('t');nu=4.0_dp;if(present(param1))nu=param1;p=student_t_cdf(x,nu);return
    end select
    if(abs(x)<1.0e-13_dp)then;p=0.5_dp;return;end if
    q=abs(x);lambda=elliptical_lambda(t,param1,param2)
    rmax=max(q+1.0_dp,integration_radius(t,param1,param2))
    n=300;h=(rmax-q)/real(n,dp);sumv=0.0_dp
    do i=0,n
      r=q+real(i,dp)*h
      if(r<=q)then;ang=0.0_dp;else;ang=acos(min(1.0_dp,q/r));end if
      if(i==0.or.i==n)then
        sumv=sumv+ang*r*elliptical_generator(r*r,t,param1,param2)
      else if(mod(i,2)==0)then
        sumv=sumv+2.0_dp*ang*r*elliptical_generator(r*r,t,param1,param2)
      else
        sumv=sumv+4.0_dp*ang*r*elliptical_generator(r*r,t,param1,param2)
      end if
    end do
    p=2.0_dp*lambda*h*sumv/3.0_dp
    if(x>0.0_dp)p=1.0_dp-p
    p=max(0.0_dp,min(1.0_dp,p))
  end function elliptical_marginal_cdf

  pure function integration_radius(type_name,param1,param2) result(rmax)
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::rmax,r,s
    select case(trim(canonical_type(type_name)))
    case('logistic');rmax=12.0_dp
    case('laplace');rmax=35.0_dp
    case('kotz');r=1.0_dp;if(present(param1))r=param1;rmax=12.0_dp/sqrt(max(r,0.02_dp))
    case('epower')
      r=1.0_dp;s=1.0_dp;if(present(param1))r=param1;if(present(param2))s=param2
      rmax=sqrt(2.0_dp)*(35.0_dp/max(r,0.02_dp))**(0.5_dp/max(s,0.1_dp))
      rmax=min(100.0_dp,max(12.0_dp,rmax))
    case default;rmax=35.0_dp
    end select
  end function integration_radius

  function elliptical_marginal_quantile(p,type_name,param1,param2) result(x)
    real(dp),intent(in)::p
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::x,lo,hi,mid,nu
    integer::i
    character(len=16)::t
    t=canonical_type(type_name)
    if(p<=0.0_dp)then;x=-huge(1.0_dp);return;end if
    if(p>=1.0_dp)then;x=huge(1.0_dp);return;end if
    select case(trim(t))
    case('norm');x=normal_quantile(p);return
    case('cauchy');x=student_t_quantile(p,1.0_dp);return
    case('t');nu=4.0_dp;if(present(param1))nu=param1;x=student_t_quantile(p,nu);return
    end select
    lo=-1.0_dp;hi=1.0_dp
    do while(elliptical_marginal_cdf(lo,t,param1,param2)>p);lo=2.0_dp*lo;end do
    do while(elliptical_marginal_cdf(hi,t,param1,param2)<p);hi=2.0_dp*hi;end do
    do i=1,65
      mid=0.5_dp*(lo+hi)
      if(elliptical_marginal_cdf(mid,t,param1,param2)<p)then;lo=mid;else;hi=mid;end if
    end do
    x=0.5_dp*(lo+hi)
  end function elliptical_marginal_quantile

  elemental function elliptical_joint_pdf(x,y,rho,type_name,param1,param2) result(f)
    real(dp),intent(in)::x,y,rho
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::f,q,den
    den=1.0_dp-rho*rho
    if(den<=0.0_dp)then;f=0.0_dp;return;end if
    q=(x*x-2.0_dp*rho*x*y+y*y)/den
    f=elliptical_lambda(type_name,param1,param2)*elliptical_generator(q,type_name,param1,param2)/sqrt(den)
  end function elliptical_joint_pdf

  function elliptical_copula_density(u,v,rho,type_name,param1,param2) result(c)
    real(dp),intent(in)::u,v,rho
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::c,x,y,fx,fy
    x=elliptical_marginal_quantile(clamp01(u),type_name,param1,param2)
    y=elliptical_marginal_quantile(clamp01(v),type_name,param1,param2)
    fx=elliptical_marginal_pdf(x,type_name,param1,param2);fy=elliptical_marginal_pdf(y,type_name,param1,param2)
    if(fx<=0.0_dp.or.fy<=0.0_dp)then;c=0.0_dp;else;c=elliptical_joint_pdf(x,y,rho,type_name,param1,param2)/(fx*fy);end if
    if(.not.ieee_is_finite(c))c=0.0_dp
  end function elliptical_copula_density

  function elliptical_copula_cdf(u,v,rho,type_name,param1,param2) result(c)
    real(dp),intent(in)::u,v,rho
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::c,x,y,nu,rmax
    type(integration_result)::res
    character(len=16)::t
    if(u<=0.0_dp.or.v<=0.0_dp)then;c=0.0_dp;return;end if
    if(u>=1.0_dp)then;c=clamp01(v);return;end if
    if(v>=1.0_dp)then;c=clamp01(u);return;end if
    t=canonical_type(type_name);x=elliptical_marginal_quantile(u,t,param1,param2);y=elliptical_marginal_quantile(v,t,param1,param2)
    select case(trim(t))
    case('norm');c=pnorm2d(x,y,rho,1.0e-9_dp);return
    case('cauchy');c=pt2d(x,y,rho,1.0_dp,1.0e-8_dp);return
    case('t');nu=4.0_dp;if(present(param1))nu=param1;c=pt2d(x,y,rho,nu,1.0e-8_dp);return
    end select
    rmax=integration_radius(t,param1,param2)
    res=adapt_integrate2d(fun,[-rmax,-rmax],[x,y],2.0e-6_dp,8)
    c=max(0.0_dp,min(min(u,v),res%value))
  contains
    function fun(a,b) result(z)
      real(dp),intent(in)::a,b;real(dp)::z
      z=elliptical_joint_pdf(a,b,rho,t,param1,param2)
    end function fun
  end function elliptical_copula_cdf

  subroutine elliptical_rng(n,rho,type_name,x,seed,param1,param2)
    integer,intent(in)::n
    real(dp),intent(in)::rho
    character(len=*),intent(in)::type_name
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    real(dp),intent(in),optional::param1,param2
    type(rng_state)::state
    integer(i8)::s0
    integer::i
    real(dp)::theta,q,radius,z1,z2,u,nu,r,s,w
    real(dp),allocatable::latent(:,:),ranks(:)
    character(len=16)::t
    t=canonical_type(type_name);s0=543210987_i8
    if(present(seed))s0=seed
    call seed_rng(state,s0)
    allocate(x(n,2),latent(n,2))
    do i=1,n
      u=uniform_rng(state)
      select case(trim(t))
      case('norm');q=-2.0_dp*log(u)
      case('cauchy');q=(u**(-2.0_dp)-1.0_dp)
      case('t');nu=4.0_dp;if(present(param1))nu=param1;q=nu*(u**(-2.0_dp/nu)-1.0_dp)
      case('logistic');q=4.0_dp*atanh(min(1.0_dp-1.0e-14_dp,u))
      case('laplace');radius=gamma_rng(state,2.0_dp,1.0_dp);q=radius*radius
      case('kotz');r=1.0_dp;if(present(param1))r=param1;q=-2.0_dp*log(u)/r
      case('epower')
        r=1.0_dp; s=1.0_dp
        if(present(param1)) r=param1
        if(present(param2)) s=param2
        w=gamma_rng(state,1.0_dp/s,1.0_dp)
        q=2.0_dp*(w/r)**(1.0_dp/s)
      case default;q=-2.0_dp*log(u)
      end select
      radius=sqrt(max(0.0_dp,q));theta=2.0_dp*pi*uniform_rng(state)
      z1=radius*cos(theta);z2=radius*sin(theta)
      latent(i,1)=z1
      latent(i,2)=rho*z1+sqrt(max(0.0_dp,1.0_dp-rho*rho))*z2
    end do
    if(trim(t)=='norm'.or.trim(t)=='cauchy'.or.trim(t)=='t')then
      do i=1,n
        x(i,1)=elliptical_marginal_cdf(latent(i,1),t,param1,param2)
        x(i,2)=elliptical_marginal_cdf(latent(i,2),t,param1,param2)
      end do
    else
      call rank_uniform(latent(:,1),ranks);x(:,1)=ranks
      call rank_uniform(latent(:,2),ranks);x(:,2)=ranks
    end if
  end subroutine elliptical_rng

  subroutine rank_uniform(v,u)
    real(dp),intent(in)::v(:)
    real(dp),allocatable,intent(out)::u(:)
    integer::i,j,n,less
    n=size(v);allocate(u(n))
    do i=1,n
      less=0
      do j=1,n
        if(v(j)<v(i))less=less+1
      end do
      u(i)=real(less+1,dp)/real(n+1,dp)
    end do
  end subroutine rank_uniform

  pure function elliptical_tau(rho) result(tau)
    real(dp),intent(in)::rho
    real(dp)::tau
    tau=2.0_dp*asin(max(-1.0_dp,min(1.0_dp,rho)))/pi
  end function elliptical_tau

  function elliptical_rho(rho,type_name,param1,param2) result(rs)
    real(dp),intent(in)::rho
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::rs
    type(integration_result)::res
    character(len=16)::t
    t=canonical_type(type_name)
    if(trim(t)=='norm')then;rs=6.0_dp*asin(rho/2.0_dp)/pi;return;end if
    res=adapt_integrate2d(fun,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],5.0e-5_dp,7)
    rs=12.0_dp*res%value-3.0_dp
  contains
    function fun(u,v) result(z)
      real(dp),intent(in)::u,v;real(dp)::z
      z=elliptical_copula_cdf(u,v,rho,t,param1,param2)
    end function fun
  end function elliptical_rho

  subroutine elliptical_tail_coeff(rho,type_name,lower,upper,param1,param2)
    real(dp),intent(in)::rho
    character(len=*),intent(in)::type_name
    real(dp),intent(out)::lower,upper
    real(dp),intent(in),optional::param1,param2
    real(dp)::nu,arg,e
    character(len=16)::t
    t=canonical_type(type_name)
    select case(trim(t))
    case('norm');lower=0.0_dp;upper=0.0_dp
    case('cauchy','t')
      nu=1.0_dp;if(trim(t)=='t')then;nu=4.0_dp;if(present(param1))nu=param1;end if
      arg=sqrt((nu+1.0_dp)*(1.0_dp-rho)/(1.0_dp+rho));lower=2.0_dp*(1.0_dp-student_t_cdf(arg,nu+1.0_dp));upper=lower
    case default
      e=1.0e-5_dp;lower=elliptical_copula_cdf(e,e,rho,t,param1,param2)/e
      upper=(1.0_dp-2.0_dp*(1.0_dp-e)+elliptical_copula_cdf(1.0_dp-e,1.0_dp-e,rho,t,param1,param2))/e
      lower=max(0.0_dp,min(1.0_dp,lower));upper=max(0.0_dp,min(1.0_dp,upper))
    end select
  end subroutine elliptical_tail_coeff

  function elliptical_fit(u,v,type_name,start_rho,start_param1,start_param2,max_iter) result(fit)
    real(dp),intent(in)::u(:),v(:)
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::start_rho,start_param1,start_param2
    integer,intent(in),optional::max_iter
    type(copula_fit_result)::fit
    type(optimizer_result)::opt
    real(dp),allocatable::x0(:),lo(:),hi(:),step(:)
    integer::k,n,niter
    character(len=16)::t
    niter=3000;if(present(max_iter))niter=max(1,max_iter)
    t=canonical_type(type_name)
    select case(trim(t))
    case('t');k=2
    case('kotz');k=2
    case('epower');k=3
    case default;k=1
    end select
    allocate(x0(k),lo(k),hi(k),step(k));lo(1)=-0.995_dp;hi(1)=0.995_dp;x0(1)=sin(0.5_dp*pi*kendall_tau_sample(u,v));step(1)=0.08_dp
    if(present(start_rho))x0(1)=start_rho
    if(k>=2)then
      select case(trim(t));case('t');lo(2)=1.01_dp;hi(2)=100.0_dp;x0(2)=4.0_dp
      case default;lo(2)=0.05_dp;hi(2)=10.0_dp;x0(2)=1.0_dp;end select
      if(present(start_param1))x0(2)=start_param1;step(2)=0.2_dp
    end if
    if(k==3)then;lo(3)=0.1_dp;hi(3)=5.0_dp;x0(3)=1.0_dp;if(present(start_param2))x0(3)=start_param2;step(3)=0.1_dp;end if
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
      real(dp),intent(in)::p(:);real(dp)::val,d,a1,a2
      integer::i,m
      a1=1.0_dp;a2=1.0_dp;if(size(p)>=2)a1=p(2);if(size(p)>=3)a2=p(3)
      val=0.0_dp;m=min(size(u),size(v))
      do i=1,m
        d=elliptical_copula_density(clamp01(u(i)),clamp01(v(i)),p(1),t,a1,a2)
        if(d<=1.0e-300_dp.or..not.ieee_is_finite(d))then;val=val+690.0_dp;else;val=val-log(d);end if
      end do
    end function obj
  end function elliptical_fit
end module fcopulae_elliptical
