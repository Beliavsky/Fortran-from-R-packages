! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_curve_ops
  use yc_kinds, only : dp
  use yc_types, only : curve_t, series_t
  use yc_utils, only : valid_positive_vector, linear_interpolate, lower_string
  use yc_splines, only : evaluate_cubic_spline, fit_cubic_spline_coefficients
  use yc_models, only : ns_rate_scalar, sv_rate_scalar, ns_forward_scalar, sv_forward_scalar
  implicit none
  private
  public :: yc_predict, yc_interpolate, yc_discount, yc_forward

contains

  function yc_predict(curve, maturities) result(out)
    type(curve_t), intent(in) :: curve
    real(dp), intent(in) :: maturities(:)
    type(series_t) :: out
    integer :: i
    if (.not. curve%ok .or. .not. valid_positive_vector(maturities)) then
      out%ok=.false.; out%message='Invalid curve or maturities.'; return
    end if
    allocate(out%x,source=maturities); allocate(out%y(size(maturities)))
    do i=1,size(maturities)
      select case(trim(curve%method))
      case('nelson_siegel')
        out%y(i)=ns_rate_scalar(maturities(i),curve%beta0,curve%beta1,curve%beta2,curve%tau)
      case('svensson')
        out%y(i)=sv_rate_scalar(maturities(i),curve%beta0,curve%beta1,curve%beta2,curve%beta3,curve%tau1,curve%tau2)
      case('cubic_spline')
        out%y(i)=evaluate_cubic_spline(curve%maturities,curve%rates,curve%spline_b,curve%spline_c,curve%spline_d,maturities(i))
      case('observed')
        out%y(i)=linear_interpolate(curve%maturities,curve%rates,maturities(i))
      case default
        out%ok=.false.;out%message='Unknown curve method.';return
      end select
    end do
  end function yc_predict

  function yc_interpolate(curve,maturities,method) result(out)
    type(curve_t),intent(in)::curve
    real(dp),intent(in)::maturities(:)
    character(len=*),intent(in),optional::method
    type(series_t)::out
    character(len=16)::meth
    real(dp),allocatable::logdf(:),b(:),c(:),d(:)
    logical::ok
    integer::i
    if(trim(curve%method)/='observed')then
      out=yc_predict(curve,maturities);return
    end if
    if(.not.valid_positive_vector(maturities))then;out%ok=.false.;out%message='Invalid maturities.';return;end if
    meth='linear';if(present(method))meth=trim(lower_string(method))
    allocate(out%x,source=maturities);allocate(out%y(size(maturities)))
    select case(meth)
    case('linear')
      do i=1,size(maturities);out%y(i)=linear_interpolate(curve%maturities,curve%rates,maturities(i));end do
    case('log_linear')
      allocate(logdf(size(curve%maturities)));logdf=-curve%rates*curve%maturities
      do i=1,size(maturities)
        out%y(i)=-linear_interpolate(curve%maturities,logdf,maturities(i))/maturities(i)
      end do
    case('cubic')
      call fit_cubic_spline_coefficients(curve%maturities,curve%rates,'natural',b,c,d,ok)
      if(.not.ok)then;out%ok=.false.;out%message='Cubic interpolation failed.';return;end if
      do i=1,size(maturities)
        out%y(i)=evaluate_cubic_spline(curve%maturities,curve%rates,b,c,d,maturities(i))
      end do
    case default
      out%ok=.false.;out%message='Unknown interpolation method.'
    end select
  end function yc_interpolate

  function yc_discount(curve,maturities,compounding) result(out)
    type(curve_t),intent(in)::curve
    real(dp),intent(in),optional::maturities(:)
    character(len=*),intent(in),optional::compounding
    type(series_t)::out
    type(series_t)::pred
    character(len=20)::comp
    comp='continuous';if(present(compounding))comp=trim(lower_string(compounding))
    if(present(maturities))then
      pred=yc_predict(curve,maturities)
    else
      pred=yc_predict(curve,curve%maturities)
    end if
    if(.not.pred%ok)then;out=pred;return;end if
    out%x=pred%x;allocate(out%y(size(pred%y)))
    select case(comp)
    case('continuous')
      out%y=exp(-pred%y*pred%x)
    case('annual')
      if(any(1.0_dp+pred%y<=0.0_dp))then;out%ok=.false.;out%message='Invalid annual rates.';return;end if
      out%y=(1.0_dp+pred%y)**(-pred%x)
    case('semi_annual')
      if(any(1.0_dp+pred%y/2.0_dp<=0.0_dp))then;out%ok=.false.;out%message='Invalid semi-annual rates.';return;end if
      out%y=(1.0_dp+pred%y/2.0_dp)**(-2.0_dp*pred%x)
    case default
      out%ok=.false.;out%message='Unknown compounding convention.'
    end select
  end function yc_discount

  function yc_forward(curve,maturities,horizon) result(out)
    type(curve_t),intent(in)::curve
    real(dp),intent(in),optional::maturities(:)
    real(dp),intent(in),optional::horizon
    type(series_t)::out
    type(series_t)::p0,p1
    real(dp),allocatable::m(:)
    real(dp)::dm,r,rp
    integer::i
    if(present(maturities))then;allocate(m,source=maturities);else;allocate(m,source=curve%maturities);end if
    if(.not.valid_positive_vector(m))then;out%ok=.false.;out%message='Invalid maturities.';return;end if
    allocate(out%x,source=m);allocate(out%y(size(m)))
    if(present(horizon))then
      if(horizon<=0.0_dp)then;out%ok=.false.;out%message='Horizon must be positive.';return;end if
      p0=yc_predict(curve,m);p1=yc_predict(curve,m+horizon)
      if(.not.p0%ok.or..not.p1%ok)then;out%ok=.false.;out%message='Prediction failed.';return;end if
      out%y=(p1%y*(m+horizon)-p0%y*m)/horizon
      return
    end if
    do i=1,size(m)
      select case(trim(curve%method))
      case('nelson_siegel')
        out%y(i)=ns_forward_scalar(m(i),curve%beta0,curve%beta1,curve%beta2,curve%tau)
      case('svensson')
        out%y(i)=sv_forward_scalar(m(i),curve%beta0,curve%beta1,curve%beta2,curve%beta3,curve%tau1,curve%tau2)
      case default
        p0=yc_predict(curve,[m(i)])
        dm=max(m(i)*1.0e-6_dp,1.0e-8_dp)
        p1=yc_predict(curve,[m(i)+dm])
        r=p0%y(1);rp=p1%y(1)
        out%y(i)=r+m(i)*(rp-r)/dm
      end select
    end do
  end function yc_forward

end module yc_curve_ops
