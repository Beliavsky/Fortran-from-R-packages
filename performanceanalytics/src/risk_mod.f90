! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module risk_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, skewness_value, kurtosis_value, quantile_type7
  use probability_mod, only: normal_pdf, normal_quantile, chi_square_cdf_1, chi_square_cdf_2
  implicit none
  private
  public :: value_at_risk, expected_shortfall, historical_var, historical_es
  public :: gaussian_var, gaussian_es, modified_var, modified_es
  public :: var_backtest_result, var_backtest, component_var, marginal_var
  public :: sample_var, sample_cvar, rachev_tail_ratio
  type :: var_backtest_result
    integer :: observations=0
    integer :: violations=0
    real(dp) :: expected_violations=0.0_dp
    real(dp) :: violation_ratio=0.0_dp
    real(dp) :: kupiec_lr=0.0_dp
    real(dp) :: kupiec_pvalue=1.0_dp
    real(dp) :: independence_lr=0.0_dp
    real(dp) :: independence_pvalue=1.0_dp
    real(dp) :: conditional_lr=0.0_dp
    real(dp) :: conditional_pvalue=1.0_dp
  end type var_backtest_result
contains
  real(dp) function historical_var(r,p) result(v)
    real(dp),intent(in)::r(:),p
    v=-quantile_type7(r,1.0_dp-p)
  end function historical_var

  real(dp) function historical_es(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp)::q
    integer::n
    q=quantile_type7(r,1.0_dp-p); n=count(r<=q)
    if(n==0) then; v=-q; else; v=-sum(r,mask=r<=q)/real(n,dp); end if
  end function historical_es

  pure real(dp) function gaussian_var(r,p) result(v)
    real(dp),intent(in)::r(:),p
    v=-mean_value(r)-normal_quantile(1.0_dp-p)*sd_value(r)
  end function gaussian_var

  pure real(dp) function gaussian_es(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp)::z
    z=normal_quantile(1.0_dp-p)
    v=-mean_value(r)+sd_value(r)*normal_pdf(z)/(1.0_dp-p)
  end function gaussian_es

  pure real(dp) function modified_var(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp)::z,s,k,zcf
    z=normal_quantile(1.0_dp-p); s=skewness_value(r,1); k=kurtosis_value(r,1,.true.)
    zcf=z+(z*z-1.0_dp)*s/6.0_dp+(z**3-3.0_dp*z)*k/24.0_dp- &
      (2.0_dp*z**3-5.0_dp*z)*s*s/36.0_dp
    v=-mean_value(r)-zcf*sd_value(r)
  end function modified_var

  pure real(dp) function modified_es(r,p) result(v)
    real(dp),intent(in)::r(:),p
    real(dp)::z,h,s,k,e
    z=normal_quantile(1.0_dp-p); s=skewness_value(r,1); k=kurtosis_value(r,1,.true.)
    h=z+(z*z-1.0_dp)*s/6.0_dp+(z**3-3.0_dp*z)*k/24.0_dp- &
      (2.0_dp*z**3-5.0_dp*z)*s*s/36.0_dp
    e=normal_pdf(h)*(1.0_dp+h**3*s/6.0_dp+ &
      (h**6-9.0_dp*h**4+9.0_dp*h*h+3.0_dp)*s*s/72.0_dp+ &
      (h**4-2.0_dp*h*h-1.0_dp)*k/24.0_dp)/(1.0_dp-p)
    v=-mean_value(r)-sd_value(r)*min(-e,h)
  end function modified_es

  real(dp) function value_at_risk(r,p,method) result(v)
    real(dp),intent(in)::r(:),p
    character(len=*),intent(in),optional::method
    character(len=24)::m
    m='modified'; if(present(method)) m=trim(adjustl(method))
    select case(trim(m))
    case('historical','historic'); v=historical_var(r,p)
    case('gaussian','normal'); v=gaussian_var(r,p)
    case default; v=modified_var(r,p)
    end select
  end function value_at_risk

  real(dp) function expected_shortfall(r,p,method) result(v)
    real(dp),intent(in)::r(:),p
    character(len=*),intent(in),optional::method
    character(len=24)::m
    m='modified'; if(present(method)) m=trim(adjustl(method))
    select case(trim(m))
    case('historical','historic'); v=historical_es(r,p)
    case('gaussian','normal'); v=gaussian_es(r,p)
    case default; v=modified_es(r,p)
    end select
  end function expected_shortfall

  real(dp) function sample_var(r,p,lower_tail) result(v)
    real(dp),intent(in)::r(:),p
    logical,intent(in),optional::lower_tail
    logical::lo
    lo=.true.; if(present(lower_tail)) lo=lower_tail
    if(lo) then; v=quantile_type7(r,p); else; v=quantile_type7(r,1.0_dp-p); end if
  end function sample_var

  real(dp) function sample_cvar(r,p,lower_tail) result(v)
    real(dp),intent(in)::r(:),p
    logical,intent(in),optional::lower_tail
    logical::lo
    real(dp)::q
    integer::n
    lo=.true.; if(present(lower_tail)) lo=lower_tail
    if(lo) then
      q=quantile_type7(r,p); n=count(r<=q)
      if(n==0) then; v=q; else; v=sum(r,mask=r<=q)/real(n,dp); end if
    else
      q=quantile_type7(r,1.0_dp-p); n=count(r>=q)
      if(n==0) then; v=q; else; v=sum(r,mask=r>=q)/real(n,dp); end if
    end if
  end function sample_cvar

  real(dp) function rachev_tail_ratio(r,p_lower,p_upper) result(v)
    real(dp),intent(in)::r(:),p_lower,p_upper
    real(dp)::lower,upper
    lower=-sample_cvar(r,p_lower,.true.); upper=sample_cvar(r,p_upper,.false.)
    if(lower<=tiny(1.0_dp)) then; v=huge(1.0_dp); else; v=upper/lower; end if
  end function rachev_tail_ratio

  subroutine var_backtest(r,var_series,p,result)
    real(dp),intent(in)::r(:),var_series(:),p
    type(var_backtest_result),intent(out)::result
    logical,allocatable::hit(:)
    integer::n,i,n00,n01,n10,n11,x
    real(dp)::phat,lruc,p01,p11,p1,ll0,ll1,eps
    n=min(size(r),size(var_series)); allocate(hit(n)); hit=r(:n)<-var_series(:n)
    x=count(hit); result%observations=n; result%violations=x
    result%expected_violations=real(n,dp)*(1.0_dp-p)
    if(result%expected_violations>0.0_dp) result%violation_ratio=real(x,dp)/result%expected_violations
    eps=1.0e-14_dp; phat=max(eps,min(1.0_dp-eps,real(x,dp)/max(real(n,dp),1.0_dp)))
    lruc=-2.0_dp*((real(n-x,dp)*log(max(p,eps))+real(x,dp)*log(max(1.0_dp-p,eps)))- &
      (real(n-x,dp)*log(1.0_dp-phat)+real(x,dp)*log(phat)))
    result%kupiec_lr=max(0.0_dp,lruc); result%kupiec_pvalue=1.0_dp-chi_square_cdf_1(result%kupiec_lr)
    n00=0;n01=0;n10=0;n11=0
    do i=2,n
      if(.not.hit(i-1).and..not.hit(i)) n00=n00+1
      if(.not.hit(i-1).and.hit(i)) n01=n01+1
      if(hit(i-1).and..not.hit(i)) n10=n10+1
      if(hit(i-1).and.hit(i)) n11=n11+1
    end do
    p01=real(n01,dp)/max(real(n00+n01,dp),1.0_dp)
    p11=real(n11,dp)/max(real(n10+n11,dp),1.0_dp)
    p1=real(n01+n11,dp)/max(real(n00+n01+n10+n11,dp),1.0_dp)
    p01=max(eps,min(1.0_dp-eps,p01)); p11=max(eps,min(1.0_dp-eps,p11)); p1=max(eps,min(1.0_dp-eps,p1))
    ll0=real(n00+n10,dp)*log(1.0_dp-p1)+real(n01+n11,dp)*log(p1)
    ll1=real(n00,dp)*log(1.0_dp-p01)+real(n01,dp)*log(p01)+ &
      real(n10,dp)*log(1.0_dp-p11)+real(n11,dp)*log(p11)
    result%independence_lr=max(0.0_dp,-2.0_dp*(ll0-ll1))
    result%independence_pvalue=1.0_dp-chi_square_cdf_1(result%independence_lr)
    result%conditional_lr=result%kupiec_lr+result%independence_lr
    result%conditional_pvalue=1.0_dp-chi_square_cdf_2(result%conditional_lr)
  end subroutine var_backtest

  real(dp) function marginal_var(asset_returns,weights,p,method,asset_index) result(v)
    real(dp),intent(in)::asset_returns(:,:),weights(:),p
    character(len=*),intent(in)::method
    integer,intent(in)::asset_index
    real(dp),allocatable::wp(:),wm(:),rp(:),rm(:)
    real(dp)::h
    h=1.0e-5_dp; allocate(wp(size(weights)),wm(size(weights)))
    allocate(rp(size(asset_returns,1)),rm(size(asset_returns,1)))
    wp=weights;wm=weights;wp(asset_index)=wp(asset_index)+h;wm(asset_index)=wm(asset_index)-h
    rp=matmul(asset_returns,wp);rm=matmul(asset_returns,wm)
    v=(value_at_risk(rp,p,method)-value_at_risk(rm,p,method))/(2.0_dp*h)
  end function marginal_var

  subroutine component_var(asset_returns,weights,p,method,components,total)
    real(dp),intent(in)::asset_returns(:,:),weights(:),p
    character(len=*),intent(in)::method
    real(dp),intent(out)::components(:),total
    real(dp),allocatable::pr(:)
    integer::i
    allocate(pr(size(asset_returns,1)));pr=matmul(asset_returns,weights)
    total=value_at_risk(pr,p,method)
    do i=1,min(size(weights),size(components))
      components(i)=weights(i)*marginal_var(asset_returns,weights,p,method,i)
    end do
  end subroutine component_var
end module risk_mod
