! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module performance_ratios_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, skewness_value, kurtosis_value
  use statistics_mod, only: lower_partial_moment, upper_partial_moment, mean_absolute_deviation
  use probability_mod, only: normal_cdf, normal_quantile
  use returns_mod, only: annualized_return, cumulative_return
  use drawdown_mod, only: max_drawdown, average_drawdown, pain_index, ulcer_index
  use drawdown_mod, only: conditional_drawdown_at_risk, drawdown_episode, find_drawdowns
  use risk_mod, only: rachev_tail_ratio, modified_var, modified_es
  implicit none
  private
  public :: sharpe_ratio, annualized_sharpe_ratio, adjusted_sharpe_ratio
  public :: probabilistic_sharpe_ratio, minimum_track_record_length
  public :: sortino_ratio, downside_deviation, upside_risk, upside_potential_ratio
  public :: downside_frequency, upside_frequency, omega_ratio, omega_sharpe_ratio
  public :: kappa_ratio, downside_sharpe_ratio, bernardo_ledoit_ratio, prospect_ratio
  public :: calmar_ratio, sterling_ratio, burke_ratio, martin_ratio, pain_ratio
  public :: information_ratio, tracking_error, active_premium, modigliani_m2
  public :: skewness_kurtosis_ratio, volatility_skewness, dratio, smoothing_index
  public :: hurst_index, kelly_ratio, modified_sharpe_ratio
contains
  pure real(dp) function sharpe_ratio(r,rf) result(v)
    real(dp),intent(in)::r(:),rf(:)
    real(dp),allocatable::e(:)
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n)
    if(sd_value(e)<=tiny(1.0_dp)) then;v=0.0_dp;else;v=mean_value(e)/sd_value(e);end if
  end function sharpe_ratio

  pure real(dp) function annualized_sharpe_ratio(r,rf,scale,geometric) result(v)
    real(dp),intent(in)::r(:),rf(:),scale
    logical,intent(in),optional::geometric
    real(dp),allocatable::e(:)
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n)
    if(sd_value(e)<=tiny(1.0_dp)) then;v=0.0_dp;return;end if
    if(present(geometric)) then
      if(geometric) then;v=annualized_return(e,scale,.true.)/(sd_value(e)*sqrt(scale));else;v=mean_value(e)*sqrt(scale)/sd_value(e);end if
    else
      v=annualized_return(e,scale,.true.)/(sd_value(e)*sqrt(scale))
    end if
  end function annualized_sharpe_ratio

  pure real(dp) function adjusted_sharpe_ratio(r,rf) result(v)
    real(dp),intent(in)::r(:),rf(:)
    real(dp)::sr,s,k
    real(dp),allocatable::e(:)
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n)
    sr=sharpe_ratio(r,rf);s=skewness_value(e,1);k=kurtosis_value(e,1,.true.)
    v=sr*(1.0_dp+s*sr/6.0_dp-k*sr*sr/24.0_dp)
  end function adjusted_sharpe_ratio

  pure real(dp) function probabilistic_sharpe_ratio(r,rf,benchmark_sr) result(v)
    real(dp),intent(in)::r(:),rf(:),benchmark_sr
    real(dp),allocatable::e(:)
    real(dp)::sr,s,k,den
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n)
    sr=sharpe_ratio(r,rf);s=skewness_value(e,1);k=kurtosis_value(e,1,.false.)
    den=sqrt(max(1.0e-15_dp,1.0_dp-s*sr+(k-1.0_dp)*sr*sr/4.0_dp))
    v=normal_cdf((sr-benchmark_sr)*sqrt(real(max(1,n-1),dp))/den)
  end function probabilistic_sharpe_ratio

  pure real(dp) function minimum_track_record_length(r,rf,benchmark_sr,confidence) result(v)
    real(dp),intent(in)::r(:),rf(:),benchmark_sr,confidence
    real(dp),allocatable::e(:)
    real(dp)::sr,s,k,z,den
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n)
    sr=sharpe_ratio(r,rf);s=skewness_value(e,1);k=kurtosis_value(e,1,.false.);z=normal_quantile(confidence)
    den=(sr-benchmark_sr)**2
    if(den<=tiny(1.0_dp)) then;v=huge(1.0_dp);else;v=1.0_dp+(1.0_dp-s*sr+(k-1.0_dp)*sr*sr/4.0_dp)*(z*z)/den;end if
  end function minimum_track_record_length

  pure real(dp) function downside_deviation(r,mar,subset) result(v)
    real(dp),intent(in)::r(:),mar
    logical,intent(in),optional::subset
    logical::sub
    integer::n
    sub=.false.;if(present(subset))sub=subset
    if(sub) then;n=count(r<mar);else;n=size(r);end if
    if(n==0) then;v=0.0_dp;else;v=sqrt(sum((min(r-mar,0.0_dp))**2)/real(n,dp));end if
  end function downside_deviation

  pure real(dp) function upside_risk(r,mar,subset) result(v)
    real(dp),intent(in)::r(:),mar
    logical,intent(in),optional::subset
    logical::sub
    integer::n
    sub=.false.;if(present(subset))sub=subset
    if(sub) then;n=count(r>mar);else;n=size(r);end if
    if(n==0) then;v=0.0_dp;else;v=sqrt(sum((max(r-mar,0.0_dp))**2)/real(n,dp));end if
  end function upside_risk

  pure real(dp) function sortino_ratio(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    real(dp)::d
    d=downside_deviation(r,mar)
    if(d<=tiny(1.0_dp))then;v=0.0_dp;else;v=(mean_value(r)-mar)/d;end if
  end function sortino_ratio

  pure real(dp) function upside_potential_ratio(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    real(dp)::d
    d=downside_deviation(r,mar)
    if(d<=tiny(1.0_dp))then;v=0.0_dp;else;v=upper_partial_moment(r,mar,1)/d;end if
  end function upside_potential_ratio

  pure real(dp) function downside_frequency(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    if(size(r)==0)then;v=0.0_dp;else;v=real(count(r<mar),dp)/real(size(r),dp);end if
  end function downside_frequency

  pure real(dp) function upside_frequency(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    if(size(r)==0)then;v=0.0_dp;else;v=real(count(r>mar),dp)/real(size(r),dp);end if
  end function upside_frequency

  pure real(dp) function omega_ratio(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    real(dp)::l
    l=lower_partial_moment(r,mar,1)
    if(l<=tiny(1.0_dp))then;v=huge(1.0_dp);else;v=upper_partial_moment(r,mar,1)/l;end if
  end function omega_ratio

  pure real(dp) function omega_sharpe_ratio(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    real(dp)::o
    o=omega_ratio(r,mar);v=o-1.0_dp
  end function omega_sharpe_ratio

  pure real(dp) function kappa_ratio(r,mar,order) result(v)
    real(dp),intent(in)::r(:),mar
    integer,intent(in)::order
    real(dp)::l
    l=lower_partial_moment(r,mar,order)
    if(l<=tiny(1.0_dp))then;v=0.0_dp;else;v=(mean_value(r)-mar)/l**(1.0_dp/real(order,dp));end if
  end function kappa_ratio

  pure real(dp) function downside_sharpe_ratio(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    real(dp)::d
    d=downside_deviation(r,mar)
    if(d<=tiny(1.0_dp))then;v=0.0_dp;else;v=mean_value(r-mar)/d;end if
  end function downside_sharpe_ratio

  pure real(dp) function bernardo_ledoit_ratio(r) result(v)
    real(dp),intent(in)::r(:)
    real(dp)::loss
    loss=-sum(min(r,0.0_dp))
    if(loss<=tiny(1.0_dp))then;v=huge(1.0_dp);else;v=sum(max(r,0.0_dp))/loss;end if
  end function bernardo_ledoit_ratio

  pure real(dp) function prospect_ratio(r,mar,loss_aversion) result(v)
    real(dp),intent(in)::r(:),mar,loss_aversion
    real(dp),allocatable::d(:)
    real(dp)::risk
    allocate(d(size(r)));d=r-mar;risk=downside_deviation(r,mar)
    if(risk<=tiny(1.0_dp))then;v=0.0_dp;else;v=(sum(max(d,0.0_dp))+loss_aversion*sum(min(d,0.0_dp)))/real(size(r),dp)/risk;end if
  end function prospect_ratio

  real(dp) function calmar_ratio(r,scale) result(v)
    real(dp),intent(in)::r(:),scale
    real(dp)::d
    d=max_drawdown(r)
    if(d<=tiny(1.0_dp))then;v=0.0_dp;else;v=annualized_return(r,scale)/d;end if
  end function calmar_ratio

  real(dp) function sterling_ratio(r,scale,excess) result(v)
    real(dp),intent(in)::r(:),scale,excess
    real(dp)::d
    d=average_drawdown(r)+excess
    if(d<=tiny(1.0_dp))then;v=0.0_dp;else;v=annualized_return(r,scale)/d;end if
  end function sterling_ratio

  real(dp) function burke_ratio(r,scale,modified) result(v)
    real(dp),intent(in)::r(:),scale
    logical,intent(in),optional::modified
    real(dp),allocatable::dd(:)
    real(dp)::den
    integer::i,n
    type(drawdown_episode),allocatable::ep(:)
    logical::m
    call find_drawdowns(r,ep,n);m=.false.;if(present(modified))m=modified
    if(n==0)then;v=0.0_dp;return;end if
    allocate(dd(n));dd=[(-ep(i)%depth,i=1,n)];den=sqrt(sum(dd*dd))
    if(m)den=den/sqrt(real(n,dp))
    if(den<=tiny(1.0_dp))then;v=0.0_dp;else;v=annualized_return(r,scale)/den;end if
  end function burke_ratio

  real(dp) function martin_ratio(r,scale,rf_annual) result(v)
    real(dp),intent(in)::r(:),scale,rf_annual
    real(dp)::u
    u=ulcer_index(r)
    if(u<=tiny(1.0_dp))then;v=0.0_dp;else;v=(annualized_return(r,scale)-rf_annual)/u;end if
  end function martin_ratio

  real(dp) function pain_ratio(r,scale,rf_annual) result(v)
    real(dp),intent(in)::r(:),scale,rf_annual
    real(dp)::p
    p=pain_index(r)
    if(p<=tiny(1.0_dp))then;v=0.0_dp;else;v=(annualized_return(r,scale)-rf_annual)/p;end if
  end function pain_ratio

  pure real(dp) function tracking_error(r,b,scale) result(v)
    real(dp),intent(in)::r(:),b(:),scale
    integer::n
    n=min(size(r),size(b));v=sd_value(r(:n)-b(:n))*sqrt(scale)
  end function tracking_error

  pure real(dp) function active_premium(r,b,scale) result(v)
    real(dp),intent(in)::r(:),b(:),scale
    integer::n
    n=min(size(r),size(b));v=mean_value(r(:n)-b(:n))*scale
  end function active_premium

  pure real(dp) function information_ratio(r,b,scale) result(v)
    real(dp),intent(in)::r(:),b(:),scale
    real(dp)::te
    te=tracking_error(r,b,scale)
    if(te<=tiny(1.0_dp))then;v=0.0_dp;else;v=active_premium(r,b,scale)/te;end if
  end function information_ratio

  pure real(dp) function modigliani_m2(r,rf,b,scale) result(v)
    real(dp),intent(in)::r(:),rf(:),b(:),scale
    real(dp)::sr
    sr=annualized_sharpe_ratio(r,rf,scale,.false.);v=mean_value(rf)*scale+sr*sd_value(b)*sqrt(scale)
  end function modigliani_m2

  pure real(dp) function skewness_kurtosis_ratio(r) result(v)
    real(dp),intent(in)::r(:)
    real(dp)::k
    k=kurtosis_value(r,1,.false.)
    if(abs(k)<=tiny(1.0_dp))then;v=0.0_dp;else;v=skewness_value(r,1)/k;end if
  end function skewness_kurtosis_ratio

  pure real(dp) function volatility_skewness(r,mar) result(v)
    real(dp),intent(in)::r(:),mar
    real(dp)::d,u
    d=downside_deviation(r,mar);u=upside_risk(r,mar)
    if(u+d<=tiny(1.0_dp))then;v=0.0_dp;else;v=(u-d)/(u+d);end if
  end function volatility_skewness

  pure real(dp) function dratio(r) result(v)
    real(dp),intent(in)::r(:)
    real(dp)::loss
    loss=-sum(min(r,0.0_dp))
    if(loss<=tiny(1.0_dp))then;v=huge(1.0_dp);else;v=sum(max(r,0.0_dp))/loss;end if
  end function dratio

  pure real(dp) function smoothing_index(r) result(v)
    real(dp),intent(in)::r(:)
    real(dp)::rho
    integer::n
    n=size(r)
    if(n<3)then;v=0.0_dp;return;end if
    rho=sum((r(2:n)-mean_value(r(2:n)))*(r(1:n-1)-mean_value(r(1:n-1))))/ &
      max(sum((r(1:n-1)-mean_value(r(1:n-1)))**2),tiny(1.0_dp))
    v=rho
  end function smoothing_index

  real(dp) function hurst_index(r,min_block) result(v)
    real(dp),intent(in)::r(:)
    integer,intent(in),optional::min_block
    integer::n,b,k,nb,i,j
    real(dp),allocatable::lx(:),ly(:),seg(:),cum(:)
    real(dp)::rs,slope,intercept,den
    n=size(r);b=8;if(present(min_block))b=max(4,min_block)
    allocate(lx(32),ly(32));k=0
    do while(b<=n/2 .and. k<32)
      nb=n/b;rs=0.0_dp
      do i=1,nb
        allocate(seg(b),cum(b));seg=r((i-1)*b+1:i*b);seg=seg-mean_value(seg)
        cum(1)=seg(1);do j=2,b;cum(j)=cum(j-1)+seg(j);end do
        if(sd_value(seg)>tiny(1.0_dp))rs=rs+(maxval(cum)-minval(cum))/sd_value(seg)
        deallocate(seg,cum)
      end do
      rs=rs/real(nb,dp)
      if(rs>0.0_dp)then;k=k+1;lx(k)=log(real(b,dp));ly(k)=log(rs);end if
      b=b*2
    end do
    if(k<2)then;v=0.5_dp;return;end if
    den=sum((lx(:k)-mean_value(lx(:k)))**2)
    if(den<=tiny(1.0_dp))then;v=0.5_dp;else;slope=sum((lx(:k)-mean_value(lx(:k)))*(ly(:k)-mean_value(ly(:k))))/den;intercept=mean_value(ly(:k))-slope*mean_value(lx(:k));v=slope+0.0_dp*intercept;end if
  end function hurst_index

  pure real(dp) function kelly_ratio(r,rf) result(v)
    real(dp),intent(in)::r(:),rf(:)
    real(dp),allocatable::e(:)
    real(dp)::var
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n);var=sd_value(e)**2
    if(var<=tiny(1.0_dp))then;v=0.0_dp;else;v=mean_value(e)/var;end if
  end function kelly_ratio

  real(dp) function modified_sharpe_ratio(r,rf,p,use_es) result(v)
    real(dp),intent(in)::r(:),rf(:),p
    logical,intent(in),optional::use_es
    real(dp),allocatable::e(:)
    real(dp)::risk
    logical::es
    integer::n
    n=min(size(r),size(rf));allocate(e(n));e=r(:n)-rf(:n);es=.false.;if(present(use_es))es=use_es
    if(es)then;risk=modified_es(e,p);else;risk=modified_var(e,p);end if
    if(risk<=tiny(1.0_dp))then;v=0.0_dp;else;v=mean_value(e)/risk;end if
  end function modified_sharpe_ratio
end module performance_ratios_mod
