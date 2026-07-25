! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_skew
  use fmultivar_kinds, only : dp, i8
  use fmultivar_rng, only : rng_state, seed_rng, normal_rng, chi_square_rng
  use fmultivar_special, only : normal_cdf, student_t_cdf
  use fmultivar_linalg, only : cholesky_lower, inverse_spd, sample_mean_cov, inverse_general
  use fmultivar_distributions, only : mvnorm_pdf, mvt_pdf
  use fmultivar_optimizer, only : optimizer_result, nelder_mead
  implicit none
  private
  public :: skew_fit_result
  public :: mvsnorm_pdf, mvsnorm_rng, mvsnorm_rect_prob
  public :: mvst_pdf, mvst_rng, mvst_rect_prob
  public :: fit_multivariate_normal, fit_skew_normal, fit_skew_t, fit_skew_cauchy, mv_fit

  type :: skew_fit_result
    real(dp), allocatable :: location(:)
    real(dp), allocatable :: omega(:,:)
    real(dp), allocatable :: alpha(:)
    real(dp) :: nu = huge(1.0_dp)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp), allocatable :: hessian(:,:)
    real(dp), allocatable :: covariance(:,:)
    integer :: iterations = 0
    logical :: converged = .false.
  end type skew_fit_result
contains
  function mvsnorm_pdf(x,mu,omega,alpha,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:)
    logical,intent(out),optional::ok
    real(dp)::f,zarg
    real(dp),allocatable::scale(:)
    logical::good
    integer::i,p
    p=size(x);allocate(scale(p))
    do i=1,p;scale(i)=sqrt(max(omega(i,i),tiny(1.0_dp)));end do
    zarg=dot_product(alpha,(x-mu)/scale)
    f=2.0_dp*mvnorm_pdf(x,mu,omega,good)*normal_cdf(zarg)
    if(present(ok))ok=good
  end function mvsnorm_pdf

  function mvst_pdf(x,mu,omega,alpha,nu,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:),nu
    logical,intent(out),optional::ok
    real(dp)::f,zarg,q,arg
    real(dp),allocatable::scale(:),inv(:,:),d(:)
    logical::good
    integer::i,p
    p=size(x);allocate(scale(p),d(p))
    do i=1,p;scale(i)=sqrt(max(omega(i,i),tiny(1.0_dp)));end do
    d=x-mu
    call inverse_spd(omega,inv,good)
    if(.not.good)then;f=0.0_dp;if(present(ok))ok=.false.;return;end if
    q=dot_product(d,matmul(inv,d));zarg=dot_product(alpha,d/scale)
    arg=zarg*sqrt((nu+real(p,dp))/(nu+q))
    f=2.0_dp*mvt_pdf(x,mu,omega,nu,good)*student_t_cdf(arg,nu+real(p,dp))
    if(present(ok))ok=good
  end function mvst_pdf

  subroutine mvsnorm_rng(n,mu,omega,alpha,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::cor(:,:),scale(:),delta(:),rescov(:,:),l(:,:),z(:),u(:)
    real(dp)::den,z0
    type(rng_state)::state
    integer::p,i,j
    logical::good
    integer(i8)::seed_use
    p=size(mu);allocate(cor(p,p),scale(p),delta(p),rescov(p,p),z(p),u(p),x(n,p))
    do i=1,p;scale(i)=sqrt(max(omega(i,i),tiny(1.0_dp)));end do
    do i=1,p;do j=1,p;cor(i,j)=omega(i,j)/(scale(i)*scale(j));end do;end do
    den=sqrt(1.0_dp+dot_product(alpha,matmul(cor,alpha)))
    delta=matmul(cor,alpha)/den
    rescov=cor-matmul(reshape(delta,[p,1]),reshape(delta,[1,p]))
    do i=1,p;rescov(i,i)=rescov(i,i)+1.0e-12_dp;end do
    call cholesky_lower(rescov,l,good)
    if(.not.good)then;x=0.0_dp;if(present(ok))ok=.false.;return;end if
    seed_use=314159265_i8;if(present(seed))seed_use=seed
    call seed_rng(state,seed_use)
    do i=1,n
      z0=normal_rng(state)
      do j=1,p;z(j)=normal_rng(state);end do
      u=delta*abs(z0)+matmul(l,z)
      x(i,:)=mu+scale*u
    end do
    if(present(ok))ok=.true.
  end subroutine mvsnorm_rng

  subroutine mvst_rng(n,mu,omega,alpha,nu,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:),nu
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::cor(:,:),scale(:),delta(:),rescov(:,:),l(:,:),z(:),u(:)
    real(dp)::den,z0,s
    type(rng_state)::state
    integer::p,i,j
    logical::good
    integer(i8)::seed_use
    p=size(mu);allocate(cor(p,p),scale(p),delta(p),rescov(p,p),z(p),u(p),x(n,p))
    do i=1,p;scale(i)=sqrt(max(omega(i,i),tiny(1.0_dp)));end do
    do i=1,p;do j=1,p;cor(i,j)=omega(i,j)/(scale(i)*scale(j));end do;end do
    den=sqrt(1.0_dp+dot_product(alpha,matmul(cor,alpha)))
    delta=matmul(cor,alpha)/den
    rescov=cor-matmul(reshape(delta,[p,1]),reshape(delta,[1,p]))
    do i=1,p;rescov(i,i)=rescov(i,i)+1.0e-12_dp;end do
    call cholesky_lower(rescov,l,good)
    if(.not.good)then;x=0.0_dp;if(present(ok))ok=.false.;return;end if
    seed_use=271828182_i8;if(present(seed))seed_use=seed
    call seed_rng(state,seed_use)
    do i=1,n
      z0=normal_rng(state)
      do j=1,p;z(j)=normal_rng(state);end do
      u=delta*abs(z0)+matmul(l,z)
      s=sqrt(chi_square_rng(state,nu)/nu)
      x(i,:)=mu+scale*u/s
    end do
    if(present(ok))ok=.true.
  end subroutine mvst_rng

  subroutine mvsnorm_rect_prob(lower,upper,mu,omega,alpha,prob,se,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:)
    real(dp),intent(out)::prob,se
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::draw(:,:)
    integer::n,i,count
    logical::good
    integer(i8)::seed_use
    n=100000;if(present(nsim))n=nsim
    seed_use=42424242_i8;if(present(seed))seed_use=seed
    call mvsnorm_rng(n,mu,omega,alpha,draw,seed_use,good)
    if(.not.good)then;prob=0.0_dp;se=huge(1.0_dp);if(present(ok))ok=.false.;return;end if
    count=0;do i=1,n;if(all(draw(i,:)>=lower).and.all(draw(i,:)<=upper))count=count+1;end do
    prob=real(count,dp)/real(n,dp);se=sqrt(max(prob*(1.0_dp-prob),0.0_dp)/real(n,dp))
    if(present(ok))ok=.true.
  end subroutine mvsnorm_rect_prob

  subroutine mvst_rect_prob(lower,upper,mu,omega,alpha,nu,prob,se,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:),nu
    real(dp),intent(out)::prob,se
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    real(dp),allocatable::draw(:,:)
    integer::n,i,count
    logical::good
    integer(i8)::seed_use
    n=100000;if(present(nsim))n=nsim
    seed_use=51515151_i8;if(present(seed))seed_use=seed
    call mvst_rng(n,mu,omega,alpha,nu,draw,seed_use,good)
    if(.not.good)then;prob=0.0_dp;se=huge(1.0_dp);if(present(ok))ok=.false.;return;end if
    count=0;do i=1,n;if(all(draw(i,:)>=lower).and.all(draw(i,:)<=upper))count=count+1;end do
    prob=real(count,dp)/real(n,dp);se=sqrt(max(prob*(1.0_dp-prob),0.0_dp)/real(n,dp))
    if(present(ok))ok=.true.
  end subroutine mvst_rect_prob

  function fit_multivariate_normal(x) result(fit)
    real(dp),intent(in)::x(:,:)
    type(skew_fit_result)::fit
    integer::p,i
    call sample_mean_cov(x,fit%location,fit%omega)
    p=size(x,2);allocate(fit%alpha(p));fit%alpha=0.0_dp;fit%nu=huge(1.0_dp)
    fit%loglik=0.0_dp
    do i=1,size(x,1);fit%loglik=fit%loglik+log(max(mvnorm_pdf(x(i,:),fit%location,fit%omega),tiny(1.0_dp)));end do
    fit%converged=.true.
  end function fit_multivariate_normal

  function fit_skew_normal(x,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    fit=fit_skew_core(x,.false.,huge(1.0_dp),.true.,max_iter,tol)
  end function fit_skew_normal

  function fit_skew_t(x,fixed_nu,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(in),optional::fixed_nu
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    if(present(fixed_nu))then
      fit=fit_skew_core(x,.true.,fixed_nu,.true.,max_iter,tol)
    else
      fit=fit_skew_core(x,.true.,8.0_dp,.false.,max_iter,tol)
    end if
  end function fit_skew_t

  function fit_skew_cauchy(x,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    fit=fit_skew_core(x,.true.,1.0_dp,.true.,max_iter,tol)
  end function fit_skew_cauchy

  function mv_fit(x,method,fixed_nu,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::method
    real(dp),intent(in),optional::fixed_nu
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    select case(trim(adjustl(method)))
    case('normal','norm');fit=fit_multivariate_normal(x)
    case('snorm','skew-normal');fit=fit_skew_normal(x,max_iter,tol)
    case('cauchy','skew-cauchy');fit=fit_skew_cauchy(x,max_iter,tol)
    case('st','skew-t')
      if(present(fixed_nu))then;fit=fit_skew_t(x,fixed_nu,max_iter,tol)
      else;fit=fit_skew_t(x,max_iter=max_iter,tol=tol);end if
    case default
      fit=fit_multivariate_normal(x);fit%converged=.false.
    end select
  end function mv_fit

  function fit_skew_core(x,is_t,nu_start,fix_nu,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:),nu_start
    logical,intent(in)::is_t,fix_nu
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    type(optimizer_result)::opt
    real(dp),allocatable::mean(:),cov(:,:),l(:,:),theta0(:),step(:),lower(:),upper(:)
    real(dp),allocatable::mu(:),om(:,:),al(:)
    integer::p,npar,i,j,k,ituse
    real(dp)::toluse
    logical::good
    call sample_mean_cov(x,mean,cov);p=size(x,2)
    do i=1,p;cov(i,i)=cov(i,i)+1.0e-6_dp;end do
    call cholesky_lower(cov,l,good)
    npar=p+p*(p+1)/2+p+merge(0,1,fix_nu .or. .not.is_t)
    allocate(theta0(npar),step(npar),lower(npar),upper(npar))
    k=0;theta0(1:p)=mean;k=p
    do i=1,p
      do j=1,i
        k=k+1
        if(i==j)then;theta0(k)=log(max(l(i,j),1.0e-6_dp));else;theta0(k)=l(i,j);end if
      end do
    end do
    theta0(k+1:k+p)=initial_shape(x);k=k+p
    if(is_t.and..not.fix_nu)then;k=k+1;theta0(k)=log(max(nu_start,0.2_dp));end if
    step=0.08_dp*max(1.0_dp,abs(theta0));lower=-huge(1.0_dp);upper=huge(1.0_dp)
    lower(1:p)=mean-20.0_dp*sqrt([(max(cov(i,i),1.0e-6_dp),i=1,p)])
    upper(1:p)=mean+20.0_dp*sqrt([(max(cov(i,i),1.0e-6_dp),i=1,p)])
    k=p
    do i=1,p
      do j=1,i
        k=k+1
        if(i==j)then;lower(k)=-12.0_dp;upper(k)=12.0_dp
        else;lower(k)=-50.0_dp;upper(k)=50.0_dp;end if
      end do
    end do
    lower(k+1:k+p)=-20.0_dp;upper(k+1:k+p)=20.0_dp;k=k+p
    if(is_t.and..not.fix_nu)then;k=k+1;lower(k)=log(0.2_dp);upper(k)=log(200.0_dp);end if
    ituse=2500;if(present(max_iter))ituse=max_iter
    toluse=1.0e-6_dp;if(present(tol))toluse=tol
    opt=nelder_mead(objective,theta0,step,toluse,ituse,lower,upper)
    call decode(opt%x,p,is_t,fix_nu,nu_start,mu,om,al,fit%nu)
    fit%location=mu;fit%omega=om;fit%alpha=al;fit%loglik=-opt%value
    fit%iterations=opt%iterations;fit%converged=opt%converged
    call numerical_hessian(objective,opt%x,fit%hessian)
    call inverse_general(fit%hessian,fit%covariance,good)
    if(.not.good)fit%covariance=0.0_dp
  contains
    function objective(theta) result(v)
      real(dp),intent(in)::theta(:)
      real(dp)::v,nu
      real(dp),allocatable::m(:),o(:,:),a(:)
      integer::ii
      logical::valid
      call decode(theta,p,is_t,fix_nu,nu_start,m,o,a,nu)
      valid=.true.;v=0.0_dp
      do ii=1,size(x,1)
        if(is_t)then
          v=v-log(max(mvst_pdf(x(ii,:),m,o,a,nu,valid),tiny(1.0_dp)))
        else
          v=v-log(max(mvsnorm_pdf(x(ii,:),m,o,a,valid),tiny(1.0_dp)))
        end if
        if(.not.valid)then;v=huge(1.0_dp)/100.0_dp;return;end if
      end do
    end function objective
  end function fit_skew_core


  function initial_shape(x) result(alpha0)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: alpha0(size(x,2))
    real(dp) :: m, sd, skew
    integer :: j, n
    n = size(x,1)
    do j = 1, size(x,2)
      m = sum(x(:,j))/real(n,dp)
      sd = sqrt(sum((x(:,j)-m)**2)/real(max(1,n-1),dp))
      if (sd <= sqrt(epsilon(1.0_dp))) then
        alpha0(j) = 0.0_dp
      else
        skew = sum(((x(:,j)-m)/sd)**3)/real(n,dp)
        alpha0(j) = max(-4.0_dp,min(4.0_dp,4.0_dp*skew))
      end if
    end do
  end function initial_shape

  subroutine decode(theta,p,is_t,fix_nu,nu_fixed,mu,omega,alpha,nu)
    real(dp),intent(in)::theta(:),nu_fixed
    integer,intent(in)::p
    logical,intent(in)::is_t,fix_nu
    real(dp),allocatable,intent(out)::mu(:),omega(:,:),alpha(:)
    real(dp),intent(out)::nu
    real(dp),allocatable::l(:,:)
    integer::i,j,k
    allocate(mu(p),omega(p,p),alpha(p),l(p,p));l=0.0_dp;mu=theta(1:p);k=p
    do i=1,p;do j=1,i;k=k+1;if(i==j)then;l(i,j)=exp(theta(k));else;l(i,j)=theta(k);end if;end do;end do
    omega=matmul(l,transpose(l));alpha=theta(k+1:k+p);k=k+p
    if(is_t)then
      if(fix_nu)then;nu=nu_fixed;else;nu=exp(theta(k+1));end if
    else;nu=huge(1.0_dp);end if
  end subroutine decode

  subroutine numerical_hessian(fun,x,hess)
    interface
      function fun(x) result(f)
        import dp
        real(dp),intent(in)::x(:)
        real(dp)::f
      end function fun
    end interface
    real(dp),intent(in)::x(:)
    real(dp),allocatable,intent(out)::hess(:,:)
    real(dp),allocatable::xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:),h(:)
    real(dp)::f0
    integer::n,i,j
    n=size(x);allocate(hess(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),h(n))
    h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x));h=max(h,1.0e-5_dp);f0=fun(x);hess=0.0_dp
    do i=1,n
      xp=x;xm=x;xp(i)=xp(i)+h(i);xm(i)=xm(i)-h(i)
      hess(i,i)=(fun(xp)-2.0_dp*f0+fun(xm))/(h(i)*h(i))
      do j=i+1,n
        xpp=x;xpm=x;xmp=x;xmm=x
        xpp(i)=xpp(i)+h(i);xpp(j)=xpp(j)+h(j)
        xpm(i)=xpm(i)+h(i);xpm(j)=xpm(j)-h(j)
        xmp(i)=xmp(i)-h(i);xmp(j)=xmp(j)+h(j)
        xmm(i)=xmm(i)-h(i);xmm(j)=xmm(j)-h(j)
        hess(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*h(i)*h(j))
        hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine numerical_hessian
end module fmultivar_skew
