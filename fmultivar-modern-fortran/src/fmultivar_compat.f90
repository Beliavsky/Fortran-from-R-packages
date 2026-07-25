! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
module fmultivar_compat
  use fmultivar_kinds, only : dp, i8
  use fmultivar_distributions, only : mvnorm_pdf, mvnorm_rng, mvnorm_rect_prob, &
    mvnorm_equicoordinate_quantile, mvt_pdf, mvt_rng, mvt_rect_prob, &
    mvt_equicoordinate_quantile, elliptical2d_density
  use fmultivar_skew, only : skew_fit_result, mvsnorm_pdf, mvsnorm_rng, &
    mvsnorm_rect_prob, mvst_pdf, mvst_rng, mvst_rect_prob, fit_skew_normal, &
    fit_skew_t, fit_skew_cauchy, mv_fit
  use fmultivar_integration, only : integration_result, integrate2d_rule, adapt_integrate_nd
  use fmultivar_grid, only : grid_data, binning_result, make_grid_data, square_binning, hex_binning
  implicit none
  private
  public :: dmvnorm, pmvnorm, qmvnorm, rmvnorm
  public :: dmvt, pmvt, qmvt, rmvt
  public :: dmsn, pmsn, rmsn, dmst, pmst, rmst, dmsc, pmsc, rmsc
  public :: dmvsnorm, pmvsnorm, rmvsnorm, dmvst, pmvst, rmvst
  public :: msn_fit, mst_fit, msc_fit, mvfit
  public :: delliptical2d, integrate2d, adapt, griddata, squarebinning, hexbinning
  abstract interface
    function compat_fun2d(x,y) result(z)
      import dp
      real(dp), intent(in) :: x,y
      real(dp) :: z
    end function compat_fun2d
    function compat_fun_nd(x) result(z)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: z
    end function compat_fun_nd
  end interface
contains
  function dmvnorm(x,mu,omega,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:)
    logical,intent(out),optional::ok
    real(dp)::f
    logical::good
    f=mvnorm_pdf(x,mu,omega,good);if(present(ok))ok=good
  end function dmvnorm

  subroutine pmvnorm(lower,upper,mu,omega,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:)
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer::nuse
    integer(i8)::suse
    logical::good
    nuse=100000;if(present(nsim))nuse=nsim
    suse=712367821_i8;if(present(seed))suse=seed
    call mvnorm_rect_prob(lower,upper,mu,omega,prob,error,nuse,suse,good)
    if(present(ok))ok=good
  end subroutine pmvnorm

  function qmvnorm(prob,mu,omega,nsim) result(q)
    real(dp),intent(in)::prob,mu(:),omega(:,:)
    integer,intent(in),optional::nsim
    real(dp)::q
    if(present(nsim))then;q=mvnorm_equicoordinate_quantile(prob,mu,omega,nsim)
    else;q=mvnorm_equicoordinate_quantile(prob,mu,omega);end if
  end function qmvnorm

  subroutine rmvnorm(n,mu,omega,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer(i8)::suse
    logical::good
    suse=975318642_i8;if(present(seed))suse=seed
    call mvnorm_rng(n,mu,omega,x,suse,good);if(present(ok))ok=good
  end subroutine rmvnorm

  function dmvt(x,mu,omega,nu,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),nu
    logical,intent(out),optional::ok
    real(dp)::f
    logical::good
    f=mvt_pdf(x,mu,omega,nu,good);if(present(ok))ok=good
  end function dmvt

  subroutine pmvt(lower,upper,mu,omega,nu,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),nu
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer::nuse
    integer(i8)::suse
    logical::good
    nuse=100000;if(present(nsim))nuse=nsim
    suse=832761245_i8;if(present(seed))suse=seed
    call mvt_rect_prob(lower,upper,mu,omega,nu,prob,error,nuse,suse,good)
    if(present(ok))ok=good
  end subroutine pmvt

  function qmvt(prob,mu,omega,nu,nsim) result(q)
    real(dp),intent(in)::prob,mu(:),omega(:,:),nu
    integer,intent(in),optional::nsim
    real(dp)::q
    if(present(nsim))then;q=mvt_equicoordinate_quantile(prob,mu,omega,nu,nsim)
    else;q=mvt_equicoordinate_quantile(prob,mu,omega,nu);end if
  end function qmvt

  subroutine rmvt(n,mu,omega,nu,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),nu
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer(i8)::suse
    logical::good
    suse=864209753_i8;if(present(seed))suse=seed
    call mvt_rng(n,mu,omega,nu,x,suse,good);if(present(ok))ok=good
  end subroutine rmvt

  function dmsn(x,mu,omega,alpha,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:)
    logical,intent(out),optional::ok
    real(dp)::f
    logical::good
    f=mvsnorm_pdf(x,mu,omega,alpha,good);if(present(ok))ok=good
  end function dmsn

  subroutine pmsn(lower,upper,mu,omega,alpha,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:)
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer::nuse
    integer(i8)::suse
    logical::good
    nuse=100000;if(present(nsim))nuse=nsim
    suse=42424242_i8;if(present(seed))suse=seed
    call mvsnorm_rect_prob(lower,upper,mu,omega,alpha,prob,error,nuse,suse,good)
    if(present(ok))ok=good
  end subroutine pmsn

  subroutine rmsn(n,mu,omega,alpha,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer(i8)::suse
    logical::good
    suse=314159265_i8;if(present(seed))suse=seed
    call mvsnorm_rng(n,mu,omega,alpha,x,suse,good);if(present(ok))ok=good
  end subroutine rmsn

  function dmst(x,mu,omega,alpha,nu,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:),nu
    logical,intent(out),optional::ok
    real(dp)::f
    logical::good
    f=mvst_pdf(x,mu,omega,alpha,nu,good);if(present(ok))ok=good
  end function dmst

  subroutine pmst(lower,upper,mu,omega,alpha,nu,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:),nu
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer::nuse
    integer(i8)::suse
    logical::good
    nuse=100000;if(present(nsim))nuse=nsim
    suse=51515151_i8;if(present(seed))suse=seed
    call mvst_rect_prob(lower,upper,mu,omega,alpha,nu,prob,error,nuse,suse,good)
    if(present(ok))ok=good
  end subroutine pmst

  subroutine rmst(n,mu,omega,alpha,nu,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:),nu
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    integer(i8)::suse
    logical::good
    suse=271828182_i8;if(present(seed))suse=seed
    call mvst_rng(n,mu,omega,alpha,nu,x,suse,good);if(present(ok))ok=good
  end subroutine rmst

  function dmsc(x,mu,omega,alpha,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:)
    logical,intent(out),optional::ok
    real(dp)::f
    logical::good
    f=mvst_pdf(x,mu,omega,alpha,1.0_dp,good);if(present(ok))ok=good
  end function dmsc

  subroutine pmsc(lower,upper,mu,omega,alpha,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:)
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    call pmst(lower,upper,mu,omega,alpha,1.0_dp,prob,error,nsim,seed,ok)
  end subroutine pmsc

  subroutine rmsc(n,mu,omega,alpha,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    call rmst(n,mu,omega,alpha,1.0_dp,x,seed,ok)
  end subroutine rmsc

  function dmvsnorm(x,mu,omega,alpha,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:)
    logical,intent(out),optional::ok
    real(dp)::f
    f=dmsn(x,mu,omega,alpha,ok)
  end function dmvsnorm

  subroutine pmvsnorm(lower,upper,mu,omega,alpha,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:)
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    call pmsn(lower,upper,mu,omega,alpha,prob,error,nsim,seed,ok)
  end subroutine pmvsnorm

  subroutine rmvsnorm(n,mu,omega,alpha,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    call rmsn(n,mu,omega,alpha,x,seed,ok)
  end subroutine rmvsnorm

  function dmvst(x,mu,omega,alpha,nu,ok) result(f)
    real(dp),intent(in)::x(:),mu(:),omega(:,:),alpha(:),nu
    logical,intent(out),optional::ok
    real(dp)::f
    f=dmst(x,mu,omega,alpha,nu,ok)
  end function dmvst

  subroutine pmvst(lower,upper,mu,omega,alpha,nu,prob,error,nsim,seed,ok)
    real(dp),intent(in)::lower(:),upper(:),mu(:),omega(:,:),alpha(:),nu
    real(dp),intent(out)::prob,error
    integer,intent(in),optional::nsim
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    call pmst(lower,upper,mu,omega,alpha,nu,prob,error,nsim,seed,ok)
  end subroutine pmvst

  subroutine rmvst(n,mu,omega,alpha,nu,x,seed,ok)
    integer,intent(in)::n
    real(dp),intent(in)::mu(:),omega(:,:),alpha(:),nu
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    logical,intent(out),optional::ok
    call rmst(n,mu,omega,alpha,nu,x,seed,ok)
  end subroutine rmvst

  elemental function delliptical2d(x,y,rho,type_name,param1,param2) result(f)
    real(dp),intent(in)::x,y,rho
    character(len=*),intent(in)::type_name
    real(dp),intent(in),optional::param1,param2
    real(dp)::f
    f=elliptical2d_density(x,y,rho,type_name,param1,param2)
  end function delliptical2d

  function integrate2d(fun,error_target) result(res)
    procedure(compat_fun2d)::fun
    real(dp),intent(in),optional::error_target
    type(integration_result)::res
    if(present(error_target))then
      res=integrate2d_rule(fun,error_target)
    else
      res=integrate2d_rule(fun)
    end if
  end function integrate2d

  function adapt(lower,upper,fun,tol,max_points) result(res)
    real(dp),intent(in)::lower(:),upper(:)
    procedure(compat_fun_nd)::fun
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::max_points
    type(integration_result)::res
    real(dp)::toluse
    integer::nuse
    toluse=1.0e-4_dp;if(present(tol))toluse=tol
    nuse=262144;if(present(max_points))nuse=max_points
    res=adapt_integrate_nd(fun,lower,upper,toluse,nuse)
  end function adapt

  function griddata(x,y,z) result(data)
    real(dp),intent(in)::x(:),y(:),z(:,:)
    type(grid_data)::data
    data=make_grid_data(x,y,z)
  end function griddata

  function squarebinning(x,y,bins_x,bins_y) result(out)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::bins_x,bins_y
    type(binning_result)::out
    if(present(bins_x).and.present(bins_y))then
      out=square_binning(x,y,bins_x,bins_y)
    else if(present(bins_x))then
      out=square_binning(x,y,bins_x)
    else
      out=square_binning(x,y)
    end if
  end function squarebinning

  function hexbinning(x,y,bins) result(out)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in),optional::bins
    type(binning_result)::out
    if(present(bins))then;out=hex_binning(x,y,bins);else;out=hex_binning(x,y);end if
  end function hexbinning

  function msn_fit(x,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    fit=fit_skew_normal(x,max_iter,tol)
  end function msn_fit

  function mst_fit(x,fixed_nu,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(in),optional::fixed_nu
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    if(present(fixed_nu))then
      fit=fit_skew_t(x,fixed_nu,max_iter,tol)
    else
      fit=fit_skew_t(x,max_iter=max_iter,tol=tol)
    end if
  end function mst_fit

  function msc_fit(x,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    fit=fit_skew_cauchy(x,max_iter,tol)
  end function msc_fit

  function mvfit(x,method,fixed_nu,max_iter,tol) result(fit)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::method
    real(dp),intent(in),optional::fixed_nu
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(skew_fit_result)::fit
    fit=mv_fit(x,method,fixed_nu,max_iter,tol)
  end function mvfit
end module fmultivar_compat
