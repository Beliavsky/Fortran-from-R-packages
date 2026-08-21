! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_spline
  use flexsurv_kinds, only : dp
  use flexsurv_math, only : normal_pdf, normal_cdf, normal_quantile, &
    integrate_gauss_legendre, log1p_fs
  use flexsurv_splines2ns, only : splines2ns_basis, splines2ns_dbasis
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_finite
  implicit none
  private

  integer, parameter, public :: spline_scale_hazard = 1
  integer, parameter, public :: spline_scale_odds = 2
  integer, parameter, public :: spline_scale_normal = 3
  integer, parameter, public :: spline_time_log = 1
  integer, parameter, public :: spline_time_identity = 2
  integer, parameter, public :: spline_basis_rp = 1
  integer, parameter, public :: spline_basis_splines2ns = 2

  type, public :: survspline_model
    real(dp), allocatable :: gamma(:)
    real(dp), allocatable :: knots(:)
    integer :: scale = spline_scale_hazard
    integer :: timescale = spline_time_log
    integer :: basis = spline_basis_rp
  end type survspline_model

  public :: rp_basis, rp_dbasis, rp_basis_matrix, rp_dbasis_matrix
  public :: survspline_eta, survspline_deta_dt
  public :: survspline_pdf, survspline_logpdf, survspline_cdf, survspline_survival
  public :: survspline_hazard, survspline_cumhaz
  public :: survspline_quantile, survspline_random, survspline_rmst, survspline_mean
  public :: validate_survspline

contains

  pure subroutine rp_basis(knots, x, b)
    real(dp), intent(in) :: knots(:), x
    real(dp), intent(out) :: b(size(knots))
    integer :: i, nk
    real(dp) :: first, last, lam
    nk=size(knots)
    b=0.0_dp
    if(nk<2)return
    b(1)=1.0_dp
    b(2)=x
    first=knots(1);last=knots(nk)
    do i=1,nk-2
      lam=(last-knots(i+1))/(last-first)
      b(i+2)=cube_pos(x-knots(i+1))-lam*cube_pos(x-first) &
        -(1.0_dp-lam)*cube_pos(x-last)
    end do
  end subroutine rp_basis

  pure subroutine rp_dbasis(knots, x, b)
    real(dp), intent(in) :: knots(:), x
    real(dp), intent(out) :: b(size(knots))
    integer :: i, nk
    real(dp) :: first, last, lam
    nk=size(knots)
    b=0.0_dp
    if(nk<2)return
    b(2)=1.0_dp
    first=knots(1);last=knots(nk)
    do i=1,nk-2
      lam=(last-knots(i+1))/(last-first)
      b(i+2)=dcube_pos(x-knots(i+1))-lam*dcube_pos(x-first) &
        -(1.0_dp-lam)*dcube_pos(x-last)
    end do
  end subroutine rp_dbasis

  pure subroutine rp_basis_matrix(knots,x,b)
    real(dp),intent(in)::knots(:),x(:)
    real(dp),intent(out)::b(size(x),size(knots))
    integer::i
    do i=1,size(x);call rp_basis(knots,x(i),b(i,:));end do
  end subroutine rp_basis_matrix

  pure subroutine rp_dbasis_matrix(knots,x,b)
    real(dp),intent(in)::knots(:),x(:)
    real(dp),intent(out)::b(size(x),size(knots))
    integer::i
    do i=1,size(x);call rp_dbasis(knots,x(i),b(i,:));end do
  end subroutine rp_dbasis_matrix

  pure real(dp) function survspline_eta(model,t) result(eta)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp),allocatable::b(:)
    real(dp)::z
    if(t<=0.0_dp)then
      eta=-ieee_value(0.0_dp,ieee_positive_inf)
      return
    end if
    z=time_transform(t,model%timescale)
    allocate(b(size(model%knots)))
    if(model%basis==spline_basis_splines2ns)then
      call splines2ns_basis(model%knots,z,b)
    else
      call rp_basis(model%knots,z,b)
    end if
    eta=dot_product(model%gamma,b)
  end function survspline_eta

  pure real(dp) function survspline_deta_dt(model,t) result(deta)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp),allocatable::b(:)
    real(dp)::z
    if(t<=0.0_dp)then;deta=0.0_dp;return;end if
    z=time_transform(t,model%timescale)
    allocate(b(size(model%knots)))
    if(model%basis==spline_basis_splines2ns)then
      call splines2ns_dbasis(model%knots,z,b)
    else
      call rp_dbasis(model%knots,z,b)
    end if
    deta=dot_product(model%gamma,b)*dtime_transform(t,model%timescale)
  end function survspline_deta_dt

  pure real(dp) function survspline_survival(model,t) result(s)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp)::eta
    if(t<=0.0_dp)then;s=1.0_dp;return;end if
    eta=survspline_eta(model,t)
    select case(model%scale)
    case(spline_scale_hazard)
      if(eta>700.0_dp)then;s=0.0_dp;else;s=exp(-exp(eta));end if
    case(spline_scale_odds)
      if(eta>=0.0_dp)then;s=exp(-eta)/(1.0_dp+exp(-eta));else;s=1.0_dp/(1.0_dp+exp(eta));end if
    case(spline_scale_normal)
      s=normal_cdf(-eta)
    case default
      s=0.0_dp
    end select
  end function survspline_survival

  pure real(dp) function survspline_cdf(model,t) result(p)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    p=1.0_dp-survspline_survival(model,t)
  end function survspline_cdf

  pure real(dp) function survspline_hazard(model,t) result(h)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp)::eta,deta,s
    if(t<=0.0_dp)then;h=0.0_dp;return;end if
    eta=survspline_eta(model,t);deta=survspline_deta_dt(model,t)
    select case(model%scale)
    case(spline_scale_hazard)
      h=deta*exp(eta)
    case(spline_scale_odds)
      if(eta>=0.0_dp)then;h=deta/(1.0_dp+exp(-eta));else;h=deta*exp(eta)/(1.0_dp+exp(eta));end if
    case(spline_scale_normal)
      s=normal_cdf(-eta)
      if(s<=0.0_dp)then;h=ieee_value(0.0_dp,ieee_positive_inf);else;h=deta*normal_pdf(eta)/s;end if
    case default
      h=0.0_dp
    end select
    if(h<0.0_dp)h=0.0_dp
  end function survspline_hazard

  pure real(dp) function survspline_cumhaz(model,t) result(h)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp)::eta,s
    if(t<=0.0_dp)then;h=0.0_dp;return;end if
    eta=survspline_eta(model,t)
    select case(model%scale)
    case(spline_scale_hazard)
      h=exp(eta)
    case(spline_scale_odds)
      if(eta>35.0_dp)then;h=eta;else;h=log1p_fs(exp(eta));end if
    case(spline_scale_normal)
      s=normal_cdf(-eta)
      if(s<=0.0_dp)then;h=ieee_value(0.0_dp,ieee_positive_inf);else;h=-log(s);end if
    case default;h=0.0_dp
    end select
  end function survspline_cumhaz

  pure real(dp) function survspline_logpdf(model,t) result(logf)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp)::h,s
    if(t<=0.0_dp)then;logf=-ieee_value(0.0_dp,ieee_positive_inf);return;end if
    h=survspline_hazard(model,t);s=survspline_survival(model,t)
    if(h<=0.0_dp.or.s<=0.0_dp)then;logf=-ieee_value(0.0_dp,ieee_positive_inf);else;logf=log(h)+log(s);end if
  end function survspline_logpdf

  pure real(dp) function survspline_pdf(model,t) result(f)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    f=exp(survspline_logpdf(model,t))
  end function survspline_pdf

  real(dp) function survspline_quantile(model,p) result(q)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::p
    real(dp)::lo,hi,mid,pmid
    integer::it
    if(p<=0.0_dp)then;q=0.0_dp;return;end if
    if(p>=1.0_dp)then;q=ieee_value(0.0_dp,ieee_positive_inf);return;end if
    lo=max(tiny(1.0_dp),1.0e-12_dp)
    hi=1.0_dp
    do while(survspline_cdf(model,hi)<p .and. hi<1.0e100_dp)
      hi=2.0_dp*hi
    end do
    do it=1,220
      mid=0.5_dp*(lo+hi);pmid=survspline_cdf(model,mid)
      if(pmid<p)then;lo=mid;else;hi=mid;end if
      if(abs(hi-lo)<1.0e-11_dp*max(1.0_dp,mid))exit
    end do
    q=0.5_dp*(lo+hi)
  end function survspline_quantile

  real(dp) function survspline_random(model) result(x)
    type(survspline_model),intent(in)::model
    real(dp)::u
    call random_number(u);x=survspline_quantile(model,u)
  end function survspline_random

  real(dp) function survspline_rmst(model,t,start) result(v)
    type(survspline_model),intent(in)::model
    real(dp),intent(in)::t
    real(dp),intent(in),optional::start
    real(dp)::a,b,st,s0
    st=0.0_dp;if(present(start))st=max(0.0_dp,start)
    if(t<=st)then;v=0.0_dp;return;end if
    if(ieee_is_finite(t))then;b=t;else;b=survspline_quantile(model,1.0_dp-1.0e-10_dp);end if
    a=st;v=integrate_gauss_legendre(sf,a,b,96)
    if(st>0.0_dp)then;s0=survspline_survival(model,st);if(s0>0.0_dp)v=v/s0;end if
  contains
    real(dp) function sf(x) result(s)
      real(dp),intent(in)::x;s=survspline_survival(model,x)
    end function sf
  end function survspline_rmst

  real(dp) function survspline_mean(model) result(v)
    type(survspline_model),intent(in)::model
    v=survspline_rmst(model,ieee_value(0.0_dp,ieee_positive_inf))
  end function survspline_mean

  pure logical function validate_survspline(model,tmin,tmax,ncheck) result(ok)
    type(survspline_model),intent(in)::model
    real(dp),intent(in),optional::tmin,tmax
    integer,intent(in),optional::ncheck
    real(dp)::a,b,t
    integer::i,n
    ok=.false.
    if(.not.allocated(model%gamma).or..not.allocated(model%knots))return
    if(size(model%gamma)/=size(model%knots).or.size(model%knots)<2)return
    if(any(model%knots(2:)<=model%knots(:size(model%knots)-1)))return
    a=1.0e-5_dp;b=100.0_dp;if(present(tmin))a=max(tmin,tiny(1.0_dp));if(present(tmax))b=tmax
    n=200;if(present(ncheck))n=ncheck
    do i=0,n
      t=exp(log(a)+(log(b)-log(a))*real(i,dp)/real(n,dp))
      if(survspline_deta_dt(model,t)<-1.0e-10_dp)return
    end do
    ok=.true.
  end function validate_survspline

  pure real(dp) function time_transform(t,kind) result(z)
    real(dp),intent(in)::t;integer,intent(in)::kind
    if(kind==spline_time_log)then;z=log(t);else;z=t;end if
  end function time_transform
  pure real(dp) function dtime_transform(t,kind) result(z)
    real(dp),intent(in)::t;integer,intent(in)::kind
    if(kind==spline_time_log)then;z=1.0_dp/t;else;z=1.0_dp;end if
  end function dtime_transform
  pure real(dp) function cube_pos(x) result(y)
    real(dp),intent(in)::x;if(x<=0.0_dp)then;y=0.0_dp;else;y=x*x*x;end if
  end function cube_pos
  pure real(dp) function dcube_pos(x) result(y)
    real(dp),intent(in)::x;if(x<=0.0_dp)then;y=0.0_dp;else;y=3.0_dp*x*x;end if
  end function dcube_pos

end module flexsurv_spline
