! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_gh_fit
  use fbasics_kinds, only: dp
  use fbasics_stats, only: sample_mean, sample_sd, sample_quantile
  use fbasics_optimize, only: nelder_mead_bounded, numerical_hessian
  use fbasics_linalg, only: matrix_inverse
  use fbasics_distributions, only: distribution_fit, dgh, qgh, dhyp, qhyp, dnig, qnig, &
    dght, pght, qght, rght, dsgh, psgh, qsgh, rsgh
  implicit none
  private
  type, public :: robust_quantile_moments
    real(dp) :: median=0.0_dp,iqr=0.0_dp,skewness=0.0_dp,kurtosis=0.0_dp
  end type robust_quantile_moments
  public :: fit_gh, fit_hyp, fit_ght, fit_sgh, fit_snig, fit_sght
  public :: dsnig, psnig, qsnig, rsnig, dsght, psght, qsght, rsght
  public :: gh_robust_moments, hyp_robust_moments, ght_robust_moments
  public :: nig_robust_moments, sgh_robust_moments, snig_robust_moments
  real(dp), allocatable, save :: active_data(:)
  character(len=8), save :: active_family='gh'
contains
  subroutine fit_gh(x,fit,max_iter)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;integer,intent(in),optional::max_iter
    call fit_family(x,'gh',fit,max_iter)
  end subroutine
  subroutine fit_hyp(x,fit,max_iter)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;integer,intent(in),optional::max_iter
    call fit_family(x,'hyp',fit,max_iter)
  end subroutine
  subroutine fit_ght(x,fit,max_iter)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;integer,intent(in),optional::max_iter
    call fit_family(x,'ght',fit,max_iter)
  end subroutine
  subroutine fit_sgh(x,fit,max_iter)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;integer,intent(in),optional::max_iter
    call fit_family(x,'sgh',fit,max_iter)
  end subroutine
  subroutine fit_sght(x,fit,max_iter)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;integer,intent(in),optional::max_iter
    call fit_family(x,'ght',fit,max_iter)
  end subroutine
  subroutine fit_snig(x,fit,max_iter)
    real(dp),intent(in)::x(:);type(distribution_fit),intent(out)::fit;integer,intent(in),optional::max_iter
    call fit_family(x,'snig',fit,max_iter)
  end subroutine

  subroutine fit_family(x,family,fit,max_iter)
    real(dp),intent(in)::x(:)
    character(len=*),intent(in)::family
    type(distribution_fit),intent(out)::fit
    integer,intent(in),optional::max_iter
    real(dp),allocatable::start(:),lower(:),upper(:),best(:),h(:,:),cov(:,:)
    real(dp)::m,s,fbest
    logical::conv
    integer::info
    integer::k,imax
    active_family=adjustl(family)
    if(allocated(active_data))deallocate(active_data);allocate(active_data(size(x)));active_data=x
    m=sample_mean(x);s=max(sample_sd(x),1.0e-4_dp)
    select case(trim(active_family))
    case('gh')
      k=5;allocate(start(k),lower(k),upper(k));start=[1.5_dp/s,0.0_dp,s,m,-0.5_dp]
      lower=[1.0e-4_dp,-20.0_dp/s,1.0e-5_dp,minval(x)-5*s,-15.0_dp]
      upper=[20.0_dp/s,20.0_dp/s,20.0_dp*s,maxval(x)+5*s,15.0_dp]
    case('hyp')
      k=4;allocate(start(k),lower(k),upper(k));start=[1.5_dp/s,0.0_dp,s,m]
      lower=[1.0e-4_dp,-20.0_dp/s,1.0e-5_dp,minval(x)-5*s]
      upper=[20.0_dp/s,20.0_dp/s,20.0_dp*s,maxval(x)+5*s]
    case('ght')
      k=4;allocate(start(k),lower(k),upper(k));start=[0.0_dp,s,m,8.0_dp]
      lower=[-20.0_dp/s,1.0e-5_dp,minval(x)-5*s,0.5_dp]
      upper=[20.0_dp/s,20.0_dp*s,maxval(x)+5*s,50.0_dp]
    case('sgh')
      k=3;allocate(start(k),lower(k),upper(k));start=[1.0_dp,0.0_dp,-0.5_dp]
      lower=[0.05_dp,-0.98_dp,-15.0_dp];upper=[30.0_dp,0.98_dp,15.0_dp]
    case default
      k=2;allocate(start(k),lower(k),upper(k));start=[1.0_dp,0.0_dp]
      lower=[0.05_dp,-0.98_dp];upper=[30.0_dp,0.98_dp]
    end select
    imax=1200;if(present(max_iter))imax=max_iter
    call nelder_mead_bounded(gh_fit_objective,start,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-8_dp)
    fit%family=trim(active_family);fit%parameters=best;fit%loglik=-fbest;fit%converged=conv
    fit%aic=2.0_dp*k+2.0_dp*fbest;fit%bic=log(real(size(x),dp))*k+2.0_dp*fbest
    call numerical_hessian(gh_fit_objective,best,h);call matrix_inverse(h,cov,info);fit%hessian=h;if(info==0)fit%covariance=cov
  end subroutine fit_family

  real(dp) function gh_fit_objective(p) result(v)
    real(dp),intent(in)::p(:)
    integer::i
    real(dp)::d
    v=0.0_dp
    do i=1,size(active_data)
      select case(trim(active_family))
      case('gh')
        if(p(1)<=abs(p(2)).or.p(3)<=0.0_dp)then;v=huge(1.0_dp);return;end if
        d=dgh(active_data(i),p(1),p(2),p(3),p(4),p(5))
      case('hyp')
        if(p(1)<=abs(p(2)).or.p(3)<=0.0_dp)then;v=huge(1.0_dp);return;end if
        d=dhyp(active_data(i),p(1),p(2),p(3),p(4))
      case('ght')
        if(p(2)<=0.0_dp.or.p(4)<=0.0_dp)then;v=huge(1.0_dp);return;end if
        d=dght(active_data(i),p(1),p(2),p(3),p(4))
      case('sgh')
        if(p(1)<=0.0_dp.or.abs(p(2))>=1.0_dp)then;v=huge(1.0_dp);return;end if
        d=dsgh(active_data(i),p(1),p(2),p(3))
      case default
        if(p(1)<=0.0_dp.or.abs(p(2))>=1.0_dp)then;v=huge(1.0_dp);return;end if
        d=dsgh(active_data(i),p(1),p(2),-0.5_dp)
      end select
      if(d<=tiny(1.0_dp).or.d/=d)then;v=huge(1.0_dp);return;end if
      v=v-log(d)
    end do
  end function gh_fit_objective


  real(dp) function dsnig(x,zeta,rho) result(v)
    real(dp),intent(in)::x,zeta,rho;v=dsgh(x,zeta,rho,-0.5_dp)
  end function
  real(dp) function psnig(x,zeta,rho) result(v)
    real(dp),intent(in)::x,zeta,rho;v=psgh(x,zeta,rho,-0.5_dp)
  end function
  real(dp) function qsnig(p,zeta,rho) result(v)
    real(dp),intent(in)::p,zeta,rho;v=qsgh(p,zeta,rho,-0.5_dp)
  end function
  real(dp) function rsnig(zeta,rho) result(v)
    real(dp),intent(in)::zeta,rho;v=rsgh(zeta,rho,-0.5_dp)
  end function
  real(dp) function dsght(x,beta,delta,mu,nu) result(v)
    real(dp),intent(in)::x,beta,delta,mu,nu;v=dght(x,beta,delta,mu,nu)
  end function
  real(dp) function psght(x,beta,delta,mu,nu) result(v)
    real(dp),intent(in)::x,beta,delta,mu,nu;v=pght(x,beta,delta,mu,nu)
  end function
  real(dp) function qsght(p,beta,delta,mu,nu) result(v)
    real(dp),intent(in)::p,beta,delta,mu,nu;v=qght(p,beta,delta,mu,nu)
  end function
  real(dp) function rsght(beta,delta,mu,nu) result(v)
    real(dp),intent(in)::beta,delta,mu,nu;v=rght(beta,delta,mu,nu)
  end function

  function quantile_summary(q) result(r)
    real(dp),intent(in)::q(7)
    type(robust_quantile_moments)::r
    r%median=q(4);r%iqr=q(6)-q(2)
    if(q(6)>q(2))then
      r%skewness=(q(6)+q(2)-2.0_dp*q(4))/(q(6)-q(2))
      r%kurtosis=(q(7)-q(5)+q(3)-q(1))/(q(6)-q(2))
    end if
  end function quantile_summary

  function gh_robust_moments(alpha,beta,delta,mu,lambda) result(r)
    real(dp),intent(in)::alpha,beta,delta,mu,lambda;type(robust_quantile_moments)::r;real(dp)::q(7);integer::i
    do i=1,7;q(i)=qgh(real(i,dp)/8.0_dp,alpha,beta,delta,mu,lambda);end do;r=quantile_summary(q)
  end function
  function hyp_robust_moments(alpha,beta,delta,mu) result(r)
    real(dp),intent(in)::alpha,beta,delta,mu;type(robust_quantile_moments)::r;real(dp)::q(7);integer::i
    do i=1,7;q(i)=qhyp(real(i,dp)/8.0_dp,alpha,beta,delta,mu);end do;r=quantile_summary(q)
  end function
  function ght_robust_moments(beta,delta,mu,nu) result(r)
    real(dp),intent(in)::beta,delta,mu,nu;type(robust_quantile_moments)::r;real(dp)::q(7);integer::i
    do i=1,7;q(i)=qght(real(i,dp)/8.0_dp,beta,delta,mu,nu);end do;r=quantile_summary(q)
  end function
  function nig_robust_moments(alpha,beta,delta,mu) result(r)
    real(dp),intent(in)::alpha,beta,delta,mu;type(robust_quantile_moments)::r;real(dp)::q(7);integer::i
    do i=1,7;q(i)=qnig(real(i,dp)/8.0_dp,alpha,beta,delta,mu);end do;r=quantile_summary(q)
  end function
  function sgh_robust_moments(zeta,rho,lambda) result(r)
    real(dp),intent(in)::zeta,rho,lambda;type(robust_quantile_moments)::r;real(dp)::q(7);integer::i
    do i=1,7;q(i)=qsgh(real(i,dp)/8.0_dp,zeta,rho,lambda);end do;r=quantile_summary(q)
  end function
  function snig_robust_moments(zeta,rho) result(r)
    real(dp),intent(in)::zeta,rho;type(robust_quantile_moments)::r
    r=sgh_robust_moments(zeta,rho,-0.5_dp)
  end function
end module fbasics_gh_fit
