! SPDX-License-Identifier: GPL-3.0-only
module tscopula_vtransforms
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use tscopula_kinds, only : dp
  use tscopula_math, only : beta_pdf, regularized_beta, beta_quantile, uniform_random, integrate_simpson
  implicit none
  private

  integer, parameter, public :: vt_symmetric = 1
  integer, parameter, public :: vt_degenerate = 2
  integer, parameter, public :: vt_linear = 3
  integer, parameter, public :: vt_2p = 4
  integer, parameter, public :: vt_2b = 5
  integer, parameter, public :: vt_3p = 6
  integer, parameter, public :: vt_3b = 7

  type, public :: vtransform_spec
    integer :: family = vt_linear
    real(dp) :: delta = 0.5_dp
    real(dp) :: kappa = 1.0_dp
    real(dp) :: xi = 1.0_dp
  end type vtransform_spec

  public :: vsymmetric, vdegenerate, vlinear, v2p, v2b, v3p, v3b
  public :: vtrans, vgradient, vinverse, vdownprob, stochinverse, pcoincide

  type :: coincidence_context
    type(vtransform_spec) :: transform
  end type coincidence_context
contains
  pure function vsymmetric() result(vt)
    type(vtransform_spec) :: vt
    vt%family = vt_symmetric
    vt%delta = 0.5_dp
  end function vsymmetric

  pure function vdegenerate() result(vt)
    type(vtransform_spec) :: vt
    vt%family = vt_degenerate
    vt%delta = 0.0_dp
  end function vdegenerate

  pure function vlinear(delta) result(vt)
    real(dp), intent(in), optional :: delta
    type(vtransform_spec) :: vt
    vt%family = vt_linear
    if (present(delta)) vt%delta = delta
  end function vlinear

  pure function v2p(delta,kappa) result(vt)
    real(dp),intent(in),optional::delta,kappa
    type(vtransform_spec)::vt
    vt%family=vt_2p;if(present(delta))vt%delta=delta;if(present(kappa))vt%kappa=kappa
  end function v2p

  pure function v2b(delta,kappa) result(vt)
    real(dp),intent(in),optional::delta,kappa
    type(vtransform_spec)::vt
    vt%family=vt_2b;if(present(delta))vt%delta=delta;if(present(kappa))vt%kappa=kappa
  end function v2b

  pure function v3p(delta,kappa,xi) result(vt)
    real(dp),intent(in),optional::delta,kappa,xi
    type(vtransform_spec)::vt
    vt%family=vt_3p;if(present(delta))vt%delta=delta;if(present(kappa))vt%kappa=kappa;if(present(xi))vt%xi=xi
  end function v3p

  pure function v3b(delta,kappa,xi) result(vt)
    real(dp),intent(in),optional::delta,kappa,xi
    type(vtransform_spec)::vt
    vt%family=vt_3b;if(present(delta))vt%delta=delta;if(present(kappa))vt%kappa=kappa;if(present(xi))vt%xi=xi
  end function v3b

  elemental real(dp) function vtrans(vt,u) result(value)
    type(vtransform_spec),intent(in)::vt
    real(dp),intent(in)::u
    real(dp)::d,k,x,a
    d=vt%delta;k=vt%kappa;x=vt%xi
    select case(vt%family)
    case(vt_symmetric);value=abs(2.0_dp*u-1.0_dp)
    case(vt_degenerate);value=u
    case(vt_linear)
      if(d<=0.0_dp)then;value=u
      else if(d>=1.0_dp)then;value=1.0_dp-u
      else if(u<=d)then;value=1.0_dp-u/d
      else;value=(u-d)/(1.0_dp-d);end if
    case(vt_2p)
      if(d<=0.0_dp)then;value=u
      else if(d>=1.0_dp)then;value=1.0_dp-u
      else if(u<=d)then
        if(u<=0.0_dp)then;value=1.0_dp;else;value=1.0_dp-u-(1.0_dp-d)*exp(-k*log(d/u));end if
      else
        value=u-d*exp(-(-log((1.0_dp-u)/(1.0_dp-d))/k))
      end if
    case(vt_2b)
      if(d<=0.0_dp)then;value=u
      else if(d>=1.0_dp)then;value=1.0_dp-u
      else if(u<=d)then;value=1.0_dp-u-(1.0_dp-d)*regularized_beta(u/d,k,1.0_dp/k)
      else;value=u-d*beta_quantile((1.0_dp-u)/(1.0_dp-d),k,1.0_dp/k);end if
    case(vt_3p)
      if(d<=0.0_dp)then;value=u
      else if(d>=1.0_dp)then;value=1.0_dp-u
      else if(u<=d)then
        if(u<=0.0_dp)then;value=1.0_dp;else;value=1.0_dp-u-(1.0_dp-d)*exp(-k*log(d/u)**x);end if
      else
        a=-log((1.0_dp-u)/(1.0_dp-d))/k
        value=u-d*exp(-(max(a,0.0_dp))**(1.0_dp/x))
      end if
    case(vt_3b)
      if(d<=0.0_dp)then;value=u
      else if(d>=1.0_dp)then;value=1.0_dp-u
      else if(u<=d)then;value=1.0_dp-u-(1.0_dp-d)*regularized_beta(u/d,k,x)
      else;value=u-d*beta_quantile((1.0_dp-u)/(1.0_dp-d),k,x);end if
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    value=min(max(value,0.0_dp),1.0_dp)
  end function vtrans

  elemental real(dp) function vgradient(vt,u) result(value)
    type(vtransform_spec),intent(in)::vt
    real(dp),intent(in)::u
    real(dp)::d,k,x,arg1,arg2,qb,db
    d=vt%delta;k=vt%kappa;x=vt%xi
    select case(vt%family)
    case(vt_symmetric);if(u<=0.5_dp)then;value=-2.0_dp;else;value=2.0_dp;end if
    case(vt_degenerate);value=1.0_dp
    case(vt_linear);if(u<=d)then;value=-1.0_dp/max(d,tiny(1.0_dp));else;value=1.0_dp/max(1.0_dp-d,tiny(1.0_dp));end if
    case(vt_2p)
      if(u<=d)then
        if(u<=0.0_dp)then
          if(k<1.0_dp)then;value=-huge(1.0_dp);else if(k>1.0_dp)then;value=-1.0_dp;else;value=-1.0_dp/max(d,tiny(1.0_dp));end if
        else;arg1=log(d/u);value=-1.0_dp-(1.0_dp-d)*exp(-k*arg1)*k/u;end if
      else
        if(u>=1.0_dp)then
          if(k>1.0_dp)then;value=huge(1.0_dp);else if(k<1.0_dp)then;value=1.0_dp;else;value=1.0_dp/max(1.0_dp-d,tiny(1.0_dp));end if
        else;arg2=log((1.0_dp-d)/(1.0_dp-u));value=1.0_dp+d*exp(-k*arg2)*k/(1.0_dp-u);end if
      end if
    case(vt_2b)
      if(u<=d)then;value=-1.0_dp-(1.0_dp-d)*beta_pdf(u/d,k,1.0_dp/k)/d
      else;qb=beta_quantile((1.0_dp-u)/(1.0_dp-d),k,1.0_dp/k);db=beta_pdf(qb,k,1.0_dp/k);value=1.0_dp+d/(db*(1.0_dp-d));end if
    case(vt_3p)
      if(u<=d)then
        if(u<=0.0_dp)then;value=-1.0_dp
        else;arg1=log(d/u);value=-1.0_dp-(1.0_dp-d)*exp(-k*arg1**x)*arg1**(x-1.0_dp)*x*k/u;end if
      else
        if(u>=1.0_dp)then;value=1.0_dp
        else;arg2=log((1.0_dp-d)/(1.0_dp-u));value=1.0_dp+d*exp(-(k*arg2)**(1.0_dp/x))*arg2**(1.0_dp/x-1.0_dp)*k**(1.0_dp/x)/(x*(1.0_dp-u));end if
      end if
    case(vt_3b)
      if(u<=d)then;value=-1.0_dp-(1.0_dp-d)*beta_pdf(u/d,k,x)/d
      else;qb=beta_quantile((1.0_dp-u)/(1.0_dp-d),k,x);db=beta_pdf(qb,k,x);value=1.0_dp+d/(db*(1.0_dp-d));end if
    case default;value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function vgradient

  elemental real(dp) function vinverse(vt,v) result(u)
    type(vtransform_spec),intent(in)::vt
    real(dp),intent(in)::v
    real(dp)::lo,hi,mid
    integer::iter
    select case(vt%family)
    case(vt_symmetric);u=0.5_dp*(1.0_dp-v)
    case(vt_degenerate);u=v
    case(vt_linear);u=vt%delta*(1.0_dp-v)
    case default
      lo=0.0_dp;hi=max(vt%delta,0.0_dp)
      do iter=1,100
        mid=0.5_dp*(lo+hi)
        if(vtrans(vt,mid)>v)then;lo=mid;else;hi=mid;end if
      end do
      u=0.5_dp*(lo+hi)
    end select
  end function vinverse

  elemental real(dp) function vdownprob(vt,v) result(value)
    type(vtransform_spec),intent(in)::vt;real(dp),intent(in)::v;real(dp)::u
    u=vinverse(vt,v);value=-1.0_dp/vgradient(vt,u);value=min(max(value,0.0_dp),1.0_dp)
  end function vdownprob

  function stochinverse(vt,v,w) result(u)
    type(vtransform_spec),intent(in)::vt;real(dp),intent(in)::v(:);real(dp),intent(in),optional::w(:)
    real(dp),allocatable::u(:);real(dp)::lower,prob,draw;integer::i
    allocate(u(size(v)))
    do i=1,size(v)
      lower=vinverse(vt,v(i));prob=vdownprob(vt,v(i));if(present(w))then;draw=w(i);else;draw=uniform_random();end if
      if(draw<=prob)then;u(i)=lower;else;u(i)=lower+v(i);end if
    end do
  end function stochinverse

  real(dp) function pcoincide(vt) result(value)
    type(vtransform_spec),intent(in)::vt
    type(coincidence_context)::context;real(dp)::variance
    if(vt%family==vt_symmetric)then;value=0.5_dp;return;end if
    context%transform=vt
    variance=integrate_simpson(coincidence_integrand,0.0_dp,1.0_dp,context,1.0e-8_dp)
    value=vt%delta**2+(1.0_dp-vt%delta)**2+2.0_dp*variance
  end function pcoincide

  real(dp) function coincidence_integrand(v,context_any) result(value)
    real(dp),intent(in)::v;class(*),intent(inout)::context_any
    select type(context=>context_any)
    type is(coincidence_context);value=(vdownprob(context%transform,v)-context%transform%delta)**2
    class default;value=0.0_dp
    end select
  end function coincidence_integrand
end module tscopula_vtransforms
