! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_distributions
  use fmultivar_kinds, only : dp, i8, pi, log_two_pi
  use fmultivar_rng, only : rng_state, seed_rng, normal_rng, chi_square_rng
  use fmultivar_special, only : normal_cdf, normal_quantile, student_t_cdf, &
    student_t_quantile
  use fmultivar_integration, only : integration_result, integrate_1d
  use fmultivar_linalg, only : cholesky_lower, inverse_spd, logdet_spd
  implicit none
  private
  public :: dnorm2d, pnorm2d, rnorm2d, dt2d, pt2d, rt2d
  public :: dcauchy2d, pcauchy2d, rcauchy2d, elliptical2d_density
  public :: mvnorm_logpdf, mvnorm_pdf, mvnorm_rng, mvnorm_rect_prob
  public :: mvnorm_equicoordinate_quantile
  public :: mvt_logpdf, mvt_pdf, mvt_rng, mvt_rect_prob
  public :: mvt_equicoordinate_quantile
contains
  elemental function dnorm2d(x,y,rho) result(f)
    real(dp),intent(in)::x,y,rho
    real(dp)::f,q,den
    den=1.0_dp-rho*rho
    if(den<=0.0_dp)then;f=0.0_dp;return;end if
    q=(x*x-2.0_dp*rho*x*y+y*y)/den
    f=exp(-0.5_dp*q)/(2.0_dp*pi*sqrt(den))
  end function dnorm2d

  function pnorm2d(x,y,rho,tol) result(p)
    real(dp),intent(in)::x,y,rho
    real(dp),intent(in),optional::tol
    real(dp)::p,px,eps_use
    type(integration_result)::res
    eps_use=1.0e-10_dp;if(present(tol))eps_use=tol
    if(abs(x)<=tiny(1.0_dp).and.abs(y)<=tiny(1.0_dp))then
      p=0.25_dp+asin(rho)/(2.0_dp*pi);return
    end if
    if(rho>=1.0_dp-1.0e-12_dp)then
      p=normal_cdf(min(x,y));return
    else if(rho<=-1.0_dp+1.0e-12_dp)then
      p=max(0.0_dp,normal_cdf(x)-normal_cdf(-y));return
    else if(abs(rho)<=epsilon(1.0_dp))then
      p=normal_cdf(x)*normal_cdf(y);return
    end if
    px=normal_cdf(x)
    if(px<=tiny(1.0_dp))then;p=0.0_dp;return;end if
    if(px>=1.0_dp-epsilon(1.0_dp))then;p=normal_cdf(y);return;end if
    res=integrate_1d(integrand,0.0_dp,px,eps_use)
    p=max(0.0_dp,min(1.0_dp,res%value))
  contains
    function integrand(u) result(v)
      real(dp),intent(in)::u
      real(dp)::v,z
      if(u<=tiny(1.0_dp))then;v=0.0_dp;return;end if
      z=normal_quantile(min(1.0_dp-epsilon(1.0_dp),u))
      v=normal_cdf((y-rho*z)/sqrt(1.0_dp-rho*rho))
    end function integrand
  end function pnorm2d

  subroutine rnorm2d(n,rho,x,seed)
    integer,intent(in)::n
    real(dp),intent(in)::rho
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    type(rng_state)::state
    integer::i
    integer(i8)::seed_use
    seed_use=1357911_i8
    if(present(seed))seed_use=seed
    call seed_rng(state,seed_use)
    allocate(x(n,2))
    do i=1,n
      x(i,1)=normal_rng(state)
      x(i,2)=rho*x(i,1)+sqrt(max(0.0_dp,1.0_dp-rho*rho))*normal_rng(state)
    end do
  end subroutine rnorm2d

  elemental function dt2d(x,y,rho,nu) result(f)
    real(dp),intent(in)::x,y,rho,nu
    real(dp)::f,q,den
    if(nu>1.0e12_dp)then;f=dnorm2d(x,y,rho);return;end if
    den=1.0_dp-rho*rho
    if(den<=0.0_dp .or. nu<=0.0_dp)then;f=0.0_dp;return;end if
    q=(x*x-2.0_dp*rho*x*y+y*y)/den
    f=(1.0_dp+q/nu)**(-0.5_dp*(nu+2.0_dp))/(2.0_dp*pi*sqrt(den))
  end function dt2d

  function pt2d(x,y,rho,nu,tol) result(p)
    real(dp),intent(in)::x,y,rho,nu
    real(dp),intent(in),optional::tol
    real(dp)::p,px,eps_use
    type(integration_result)::res
    eps_use=1.0e-9_dp;if(present(tol))eps_use=tol
    if(nu>1.0e12_dp)then;p=pnorm2d(x,y,rho,eps_use);return;end if
    if(rho>=1.0_dp-1.0e-12_dp)then
      p=student_t_cdf(min(x,y),nu);return
    else if(rho<=-1.0_dp+1.0e-12_dp)then
      p=max(0.0_dp,student_t_cdf(x,nu)-student_t_cdf(-y,nu));return
    else if(abs(rho)<=epsilon(1.0_dp))then
      p=student_t_cdf(x,nu)*student_t_cdf(y,nu);return
    end if
    px=student_t_cdf(x,nu)
    if(px<=tiny(1.0_dp))then;p=0.0_dp;return;end if
    if(px>=1.0_dp-epsilon(1.0_dp))then;p=student_t_cdf(y,nu);return;end if
    res=integrate_1d(integrand,0.0_dp,px,eps_use)
    p=max(0.0_dp,min(1.0_dp,res%value))
  contains
    function integrand(u) result(v)
      real(dp),intent(in)::u
      real(dp)::v,z,arg
      if(u<=tiny(1.0_dp))then;v=0.0_dp;return;end if
      z=student_t_quantile(min(1.0_dp-epsilon(1.0_dp),u),nu)
      arg=(y-rho*z)*sqrt((nu+1.0_dp)/(nu+z*z))/sqrt(1.0_dp-rho*rho)
      v=student_t_cdf(arg,nu+1.0_dp)
    end function integrand
  end function pt2d

  subroutine rt2d(n,rho,nu,x,seed)
    integer,intent(in)::n
    real(dp),intent(in)::rho,nu
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    type(rng_state)::state
    integer::i
    integer(i8)::seed_use
    real(dp)::z1,z2,s
    seed_use=2468021_i8
    if(present(seed))seed_use=seed
    call seed_rng(state,seed_use)
    allocate(x(n,2))
    do i=1,n
      z1=normal_rng(state);z2=rho*z1+sqrt(max(0.0_dp,1.0_dp-rho*rho))*normal_rng(state)
      s=sqrt(chi_square_rng(state,nu)/nu)
      x(i,:)=[z1/s,z2/s]
    end do
  end subroutine rt2d

  elemental function dcauchy2d(x,y,rho) result(f)
    real(dp),intent(in)::x,y,rho
    real(dp)::f
    f=dt2d(x,y,rho,1.0_dp)
  end function dcauchy2d

  function pcauchy2d(x,y,rho,tol) result(p)
    real(dp),intent(in)::x,y,rho
    real(dp),intent(in),optional::tol
    real(dp)::p
    if(present(tol))then;p=pt2d(x,y,rho,1.0_dp,tol);else;p=pt2d(x,y,rho,1.0_dp);end if
  end function pcauchy2d

  subroutine rcauchy2d(n,rho,x,seed)
    integer,intent(in)::n
    real(dp),intent(in)::rho
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    if(present(seed))then;call rt2d(n,rho,1.0_dp,x,seed);else;call rt2d(n,rho,1.0_dp,x);end if
  end subroutine rcauchy2d

  elemental function elliptical2d_density(x,y,rho,type_name,param1,param2) result(f)
    real(dp),intent(in)::x,y,rho
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::f,q,g,lambda,nu,r,s,den
    den=1.0_dp-rho*rho
    if(den<=0.0_dp)then;f=0.0_dp;return;end if
    q=(x*x-2.0_dp*rho*x*y+y*y)/den
    select case(trim(adjustl(type_name)))
    case('norm','normal')
      g=exp(-0.5_dp*q);lambda=1.0_dp/(2.0_dp*pi)
    case('cauchy')
      g=(1.0_dp+q)**(-1.5_dp);lambda=1.0_dp/(2.0_dp*pi)
    case('t','student')
      nu=4.0_dp;if(present(param1))nu=param1
      g=(1.0_dp+q/nu)**(-0.5_dp*(nu+2.0_dp));lambda=1.0_dp/(2.0_dp*pi)
    case('logistic')
      g=exp(-0.5_dp*q)/(1.0_dp+exp(-0.5_dp*q))**2;lambda=1.0_dp/pi
    case('laplace')
      r=sqrt(2.0_dp);s=0.5_dp
      g=exp(-r*(0.5_dp*q)**s)
      lambda=s*r**(1.0_dp/s)/(2.0_dp*pi*gamma(1.0_dp/s))
    case('kotz')
      r=sqrt(2.0_dp);if(present(param1))r=param1
      g=exp(-0.5_dp*r*q);lambda=r/(2.0_dp*pi)
    case('epower')
      r=sqrt(2.0_dp);s=0.5_dp
      if(present(param1))r=param1;if(present(param2))s=param2
      g=exp(-r*(0.5_dp*q)**s)
      lambda=s*r**(1.0_dp/s)/(2.0_dp*pi*gamma(1.0_dp/s))
    case default
      f=0.0_dp;return
    end select
    f=lambda*g/sqrt(den)
  end function elliptical2d_density

  function mvnorm_logpdf(x,mean,cov,ok) result(logf)
    real(dp),intent(in)::x(:),mean(:),cov(:,:)
    logical,intent(out),optional::ok
    real(dp)::logf,ld,q
    real(dp),allocatable::inv(:,:),d(:)
    logical::good
    call inverse_spd(cov,inv,good)
    if(.not.good)then;logf=-huge(1.0_dp);if(present(ok))ok=.false.;return;end if
    ld=logdet_spd(cov,good);allocate(d(size(x)));d=x-mean
    q=dot_product(d,matmul(inv,d))
    logf=-0.5_dp*(real(size(x),dp)*log_two_pi+ld+q)
    if(present(ok))ok=.true.
  end function mvnorm_logpdf

  function mvnorm_pdf(x,mean,cov,ok) result(f)
    real(dp),intent(in)::x(:),mean(:),cov(:,:)
    logical,intent(out),optional::ok
    real(dp)::f,lf
    logical::good
    lf=mvnorm_logpdf(x,mean,cov,good);f=merge(exp(lf),0.0_dp,good)
    if(present(ok))ok=good
  end function mvnorm_pdf

  subroutine mvnorm_rng(n,mean,cov,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),cov(:,:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::l(:,:),z(:)
    type(rng_state)::state
    logical::good
    integer::i,j,p
    integer(i8)::seed_use
    p=size(mean);call cholesky_lower(cov,l,good);allocate(x(n,p));x=0.0_dp
    if(.not.good)then;if(present(ok))ok=.false.;return;end if
    seed_use=975318642_i8;if(present(seed))seed_use=seed
    call seed_rng(state,seed_use);allocate(z(p))
    do i=1,n
      do j=1,p;z(j)=normal_rng(state);end do
      x(i,:)=mean+matmul(l,z)
    end do
    if(present(ok))ok=.true.
  end subroutine mvnorm_rng

  subroutine mvt_rng(n,mean,scale,nu,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),scale(:,:),nu
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::l(:,:),z(:)
    type(rng_state)::state
    logical::good
    integer::i,j,p
    integer(i8)::seed_use
    real(dp)::s
    p=size(mean);call cholesky_lower(scale,l,good);allocate(x(n,p));x=0.0_dp
    if(.not.good)then;if(present(ok))ok=.false.;return;end if
    seed_use=864209753_i8;if(present(seed))seed_use=seed
    call seed_rng(state,seed_use);allocate(z(p))
    do i=1,n
      do j=1,p;z(j)=normal_rng(state);end do
      s=sqrt(chi_square_rng(state,nu)/nu)
      x(i,:)=mean+matmul(l,z)/s
    end do
    if(present(ok))ok=.true.
  end subroutine mvt_rng

  function mvt_logpdf(x,mean,scale,nu,ok) result(logf)
    real(dp),intent(in)::x(:),mean(:),scale(:,:),nu
    logical,intent(out),optional::ok
    real(dp)::logf,ld,q,pd
    real(dp),allocatable::inv(:,:),d(:)
    logical::good
    call inverse_spd(scale,inv,good)
    if(.not.good .or. nu<=0.0_dp)then
      logf=-huge(1.0_dp);if(present(ok))ok=.false.;return
    end if
    ld=logdet_spd(scale,good);allocate(d(size(x)));d=x-mean;pd=real(size(x),dp)
    q=dot_product(d,matmul(inv,d))
    logf=log_gamma(0.5_dp*(nu+pd))-log_gamma(0.5_dp*nu)-0.5_dp*(pd*log(nu*pi)+ld) &
      -0.5_dp*(nu+pd)*log(1.0_dp+q/nu)
    if(present(ok))ok=.true.
  end function mvt_logpdf

  function mvt_pdf(x,mean,scale,nu,ok) result(f)
    real(dp),intent(in)::x(:),mean(:),scale(:,:),nu
    logical,intent(out),optional::ok
    real(dp)::f,lf
    logical::good
    lf=mvt_logpdf(x,mean,scale,nu,good);f=merge(exp(lf),0.0_dp,good)
    if(present(ok))ok=good
  end function mvt_pdf

  subroutine mvnorm_rect_prob(lower,upper,mean,cov,prob,se,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mean(:),cov(:,:)
    real(dp),intent(out)::prob,se
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::draw(:,:)
    integer::n,i,count
    integer(i8)::seed_use
    logical::good
    n=100000;if(present(nsim))n=nsim
    seed_use=712367821_i8;if(present(seed))seed_use=seed
    call mvnorm_rng(n,mean,cov,draw,seed_use,good)
    if(.not.good)then;prob=0.0_dp;se=huge(1.0_dp);if(present(ok))ok=.false.;return;end if
    count=0
    do i=1,n;if(all(draw(i,:)>=lower).and.all(draw(i,:)<=upper))count=count+1;end do
    prob=real(count,dp)/real(n,dp);se=sqrt(max(prob*(1.0_dp-prob),0.0_dp)/real(n,dp))
    if(present(ok))ok=.true.
  end subroutine mvnorm_rect_prob

  subroutine mvt_rect_prob(lower,upper,mean,scale,nu,prob,se,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mean(:),scale(:,:),nu
    real(dp),intent(out)::prob,se
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::draw(:,:)
    integer::n,i,count
    integer(i8)::seed_use
    logical::good
    n=100000;if(present(nsim))n=nsim
    seed_use=832761245_i8;if(present(seed))seed_use=seed
    call mvt_rng(n,mean,scale,nu,draw,seed_use,good)
    if(.not.good)then;prob=0.0_dp;se=huge(1.0_dp);if(present(ok))ok=.false.;return;end if
    count=0
    do i=1,n;if(all(draw(i,:)>=lower).and.all(draw(i,:)<=upper))count=count+1;end do
    prob=real(count,dp)/real(n,dp);se=sqrt(max(prob*(1.0_dp-prob),0.0_dp)/real(n,dp))
    if(present(ok))ok=.true.
  end subroutine mvt_rect_prob

  function mvnorm_equicoordinate_quantile(prob,mean,cov,nsim) result(q)
    real(dp),intent(in)::prob,mean(:),cov(:,:)
    integer,intent(in),optional::nsim
    real(dp)::q,lo,hi,mid,p,se
    real(dp),allocatable::lower(:),upper(:)
    integer::i,n
    n=50000;if(present(nsim))n=nsim
    allocate(lower(size(mean)),upper(size(mean)));lower=-huge(1.0_dp)
    lo=minval(mean-10.0_dp*sqrt([(cov(i,i),i=1,size(mean))]))
    hi=maxval(mean+10.0_dp*sqrt([(cov(i,i),i=1,size(mean))]))
    do i=1,45
      mid=0.5_dp*(lo+hi);upper=mid
      call mvnorm_rect_prob(lower,upper,mean,cov,p,se,n,1234567_i8)
      if(p<prob)then;lo=mid;else;hi=mid;end if
    end do
    q=0.5_dp*(lo+hi)
  end function mvnorm_equicoordinate_quantile

  function mvt_equicoordinate_quantile(prob,mean,scale,nu,nsim) result(q)
    real(dp),intent(in)::prob,mean(:),scale(:,:),nu
    integer,intent(in),optional::nsim
    real(dp)::q,lo,hi,mid,p,se
    real(dp),allocatable::lower(:),upper(:)
    integer::i,n
    n=50000;if(present(nsim))n=nsim
    allocate(lower(size(mean)),upper(size(mean)));lower=-huge(1.0_dp)
    lo=minval(mean-50.0_dp*sqrt([(scale(i,i),i=1,size(mean))]))
    hi=maxval(mean+50.0_dp*sqrt([(scale(i,i),i=1,size(mean))]))
    do i=1,45
      mid=0.5_dp*(lo+hi);upper=mid
      call mvt_rect_prob(lower,upper,mean,scale,nu,p,se,n,7654321_i8)
      if(p<prob)then;lo=mid;else;hi=mid;end if
    end do
    q=0.5_dp*(lo+hi)
  end function mvt_equicoordinate_quantile
end module fmultivar_distributions
