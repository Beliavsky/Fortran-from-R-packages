! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_standardize
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_spec, predict_survival, predict_hazard
  use flexsurv_math, only : integrate_gauss_legendre
  implicit none
  private
  public :: standsurv_survival, standsurv_hazard, standsurv_rmst
  public :: standsurv_quantile, standardized_contrast
contains

  real(dp) function standsurv_survival(spec,theta,t,weights) result(s)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    real(dp),intent(in),optional::weights(:)
    integer::i,n
    real(dp)::w,sw
    n=prediction_rows(spec);s=0.0_dp;sw=0.0_dp
    do i=1,n
      w=1.0_dp;if(present(weights))w=weights(i)
      s=s+w*predict_survival(spec,theta,i,t);sw=sw+w
    end do
    if(sw>0.0_dp)s=s/sw
  end function standsurv_survival

  real(dp) function standsurv_hazard(spec,theta,t,weights) result(h)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    real(dp),intent(in),optional::weights(:)
    integer::i,n
    real(dp)::w,sw,surv,hs
    n=prediction_rows(spec);h=0.0_dp;sw=0.0_dp
    ! Standardized hazard is the hazard of the standardized survival curve:
    ! sum w_i S_i h_i / sum w_i S_i.
    do i=1,n
      w=1.0_dp;if(present(weights))w=weights(i)
      surv=predict_survival(spec,theta,i,t);hs=predict_hazard(spec,theta,i,t)
      h=h+w*surv*hs;sw=sw+w*surv
    end do
    if(sw>0.0_dp)h=h/sw
  end function standsurv_hazard

  real(dp) function standsurv_rmst(spec,theta,t,weights) result(r)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),t
    real(dp),intent(in),optional::weights(:)
    if(t<=0.0_dp)then;r=0.0_dp;return;end if
    r=integrate_gauss_legendre(sfun,0.0_dp,t,64)
  contains
    real(dp) function sfun(x) result(v)
      real(dp),intent(in)::x
      v=standsurv_survival(spec,theta,x,weights)
    end function sfun
  end function standsurv_rmst

  real(dp) function standsurv_quantile(spec,theta,p,weights,upper) result(q)
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:),p
    real(dp),intent(in),optional::weights(:),upper
    real(dp)::lo,hi,mid
    integer::it
    if(p<=0.0_dp)then;q=0.0_dp;return;end if
    hi=1.0_dp;if(present(upper))hi=max(upper,1.0_dp)
    do while(1.0_dp-standsurv_survival(spec,theta,hi,weights)<p.and.hi<1.0e12_dp)
      hi=2.0_dp*hi
    end do
    lo=0.0_dp
    do it=1,100
      mid=0.5_dp*(lo+hi)
      if(1.0_dp-standsurv_survival(spec,theta,mid,weights)<p)then;lo=mid;else;hi=mid;end if
    end do
    q=0.5_dp*(lo+hi)
  end function standsurv_quantile

  pure real(dp) function standardized_contrast(a,b,kind) result(v)
    real(dp),intent(in)::a,b
    integer,intent(in)::kind
    select case(kind)
    case(1);v=a-b
    case(2);if(b/=0.0_dp)then;v=a/b;else;v=huge(1.0_dp);end if
    case default;v=a-b
    end select
  end function standardized_contrast

  integer function prediction_rows(spec) result(n)
    type(flexsurv_spec),intent(in)::spec
    integer::i
    n=1
    do i=1,size(spec%reg)
      if(allocated(spec%reg(i)%x))then
        if(size(spec%reg(i)%x,1)>0)n=max(n,size(spec%reg(i)%x,1))
      end if
    end do
  end function prediction_rows

end module flexsurv_standardize
