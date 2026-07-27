! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module capm_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, covariance_value, variance_value, ols_fit
  implicit none
  private
  public :: sfm_result, sfm_fit, capm_alpha, capm_beta, capm_beta_bull, capm_beta_bear
  public :: jensen_alpha, systematic_risk, specific_risk, total_risk
  public :: selectivity, net_selectivity, treynor_ratio, appraisal_ratio
  public :: market_timing_result, treynor_mazuy_fit, henriksson_merton_fit
  type :: sfm_result
    real(dp)::alpha=0.0_dp
    real(dp)::beta=0.0_dp
    real(dp)::r_squared=0.0_dp
    real(dp)::residual_sd=0.0_dp
    real(dp),allocatable::fitted(:)
    real(dp),allocatable::residuals(:)
  end type sfm_result
  type :: market_timing_result
    real(dp)::alpha=0.0_dp
    real(dp)::beta=0.0_dp
    real(dp)::gamma=0.0_dp
    real(dp)::r_squared=0.0_dp
  end type market_timing_result
contains
  subroutine sfm_fit(ra,rb,rf,result)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    type(sfm_result),intent(out)::result
    real(dp),allocatable::x(:,:),y(:),b(:),res(:),bm(:)
    real(dp)::r2
    logical::ok
    integer::n
    n=min(size(ra),size(rb));allocate(y(n),bm(n),x(n,2),b(2),res(n))
    y=ra(:n);bm=rb(:n)
    if(present(rf)) then;y=y-rf(:n);bm=bm-rf(:n);end if
    x(:,1)=1.0_dp;x(:,2)=bm
    call ols_fit(x,y,b,res,r2,ok)
    if(.not.ok) return
    result%alpha=b(1);result%beta=b(2);result%r_squared=r2;result%residual_sd=sd_value(res)
    allocate(result%fitted(n),result%residuals(n));result%fitted=matmul(x,b);result%residuals=res
  end subroutine sfm_fit

  real(dp) function capm_alpha(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    type(sfm_result)::fit
    if(present(rf)) then;call sfm_fit(ra,rb,rf,fit);else;call sfm_fit(ra,rb,result=fit);end if
    v=fit%alpha
  end function capm_alpha

  real(dp) function capm_beta(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    type(sfm_result)::fit
    if(present(rf)) then;call sfm_fit(ra,rb,rf,fit);else;call sfm_fit(ra,rb,result=fit);end if
    v=fit%beta
  end function capm_beta

  real(dp) function capm_beta_bull(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    real(dp),allocatable::a(:),b(:),f(:)
    logical,allocatable::mask(:)
    integer::n
    n=min(size(ra),size(rb));allocate(mask(n));mask=rb(:n)>0.0_dp
    allocate(a(count(mask)),b(count(mask)));a=pack(ra(:n),mask);b=pack(rb(:n),mask)
    if(present(rf)) then;allocate(f(count(mask)));f=pack(rf(:n),mask);v=capm_beta(a,b,f);else;v=capm_beta(a,b);end if
  end function capm_beta_bull

  real(dp) function capm_beta_bear(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    real(dp),allocatable::a(:),b(:),f(:)
    logical,allocatable::mask(:)
    integer::n
    n=min(size(ra),size(rb));allocate(mask(n));mask=rb(:n)<0.0_dp
    allocate(a(count(mask)),b(count(mask)));a=pack(ra(:n),mask);b=pack(rb(:n),mask)
    if(present(rf)) then;allocate(f(count(mask)));f=pack(rf(:n),mask);v=capm_beta(a,b,f);else;v=capm_beta(a,b);end if
  end function capm_beta_bear

  real(dp) function jensen_alpha(ra,rb,rf,scale) result(v)
    real(dp),intent(in)::ra(:),rb(:),rf(:),scale
    v=capm_alpha(ra,rb,rf)*scale
  end function jensen_alpha

  real(dp) function systematic_risk(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    real(dp)::beta
    if(present(rf)) then;beta=capm_beta(ra,rb,rf);else;beta=capm_beta(ra,rb);end if
    v=abs(beta)*sd_value(rb)
  end function systematic_risk

  real(dp) function specific_risk(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:)
    real(dp),intent(in),optional::rf(:)
    type(sfm_result)::fit
    if(present(rf)) then;call sfm_fit(ra,rb,rf,fit);else;call sfm_fit(ra,rb,result=fit);end if
    v=fit%residual_sd
  end function specific_risk

  pure real(dp) function total_risk(ra) result(v)
    real(dp),intent(in)::ra(:)
    v=sd_value(ra)
  end function total_risk

  real(dp) function selectivity(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:),rf(:)
    v=mean_value(ra-rf)-capm_beta(ra,rb,rf)*mean_value(rb-rf)
  end function selectivity

  real(dp) function net_selectivity(ra,rb,rf) result(v)
    real(dp),intent(in)::ra(:),rb(:),rf(:)
    real(dp)::beta,diversification
    beta=capm_beta(ra,rb,rf)
    diversification=sd_value(ra)-abs(beta)*sd_value(rb)
    v=selectivity(ra,rb,rf)-diversification*mean_value(rb-rf)/max(sd_value(rb),tiny(1.0_dp))
  end function net_selectivity

  real(dp) function treynor_ratio(ra,rb,rf,scale) result(v)
    real(dp),intent(in)::ra(:),rb(:),rf(:),scale
    real(dp)::beta
    beta=capm_beta(ra,rb,rf)
    if(abs(beta)<=tiny(1.0_dp)) then;v=0.0_dp;else;v=mean_value(ra-rf)*scale/beta;end if
  end function treynor_ratio

  real(dp) function appraisal_ratio(ra,rb,rf,scale) result(v)
    real(dp),intent(in)::ra(:),rb(:),rf(:),scale
    type(sfm_result)::fit
    call sfm_fit(ra,rb,rf,fit)
    if(fit%residual_sd<=tiny(1.0_dp)) then;v=0.0_dp;else;v=fit%alpha*sqrt(scale)/fit%residual_sd;end if
  end function appraisal_ratio

  subroutine treynor_mazuy_fit(ra,rb,rf,result)
    real(dp),intent(in)::ra(:),rb(:),rf(:)
    type(market_timing_result),intent(out)::result
    real(dp),allocatable::x(:,:),y(:),b(:),res(:),m(:)
    real(dp)::r2
    logical::ok
    integer::n
    n=min(size(ra),min(size(rb),size(rf)));allocate(x(n,3),y(n),b(3),res(n),m(n))
    y=ra(:n)-rf(:n);m=rb(:n)-rf(:n);x(:,1)=1.0_dp;x(:,2)=m;x(:,3)=m*m
    call ols_fit(x,y,b,res,r2,ok);if(.not.ok)return
    result%alpha=b(1);result%beta=b(2);result%gamma=b(3);result%r_squared=r2
  end subroutine treynor_mazuy_fit

  subroutine henriksson_merton_fit(ra,rb,rf,result)
    real(dp),intent(in)::ra(:),rb(:),rf(:)
    type(market_timing_result),intent(out)::result
    real(dp),allocatable::x(:,:),y(:),b(:),res(:),m(:)
    real(dp)::r2
    logical::ok
    integer::n
    n=min(size(ra),min(size(rb),size(rf)));allocate(x(n,3),y(n),b(3),res(n),m(n))
    y=ra(:n)-rf(:n);m=rb(:n)-rf(:n);x(:,1)=1.0_dp;x(:,2)=m;x(:,3)=max(0.0_dp,-m)
    call ols_fit(x,y,b,res,r2,ok);if(.not.ok)return
    result%alpha=b(1);result%beta=b(2);result%gamma=b(3);result%r_squared=r2
  end subroutine henriksson_merton_fit
end module capm_mod
