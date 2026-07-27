! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module analytics_extensions_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, skewness_value, kurtosis_value, ols_fit
  use returns_mod, only: annualized_return
  use drawdown_mod, only: max_drawdown
  use risk_mod, only: value_at_risk, expected_shortfall
  use performance_ratios_mod, only: sortino_ratio, information_ratio
  use capm_mod, only: sfm_result, sfm_fit
  implicit none
  private
  public :: performance_summary, compute_performance_summary
  public :: rolling_sfm, expanding_sfm, capture_ratios
  public :: conditional_capm_result, conditional_capm_fit, outperformance_probabilities

  type :: conditional_capm_result
    integer :: lags = 0
    integer :: observations = 0
    real(dp) :: r_squared = 0.0_dp
    logical :: success = .false.
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: residuals(:)
  end type conditional_capm_result

  type :: performance_summary
    integer :: observations = 0
    real(dp) :: mean_return = 0.0_dp
    real(dp) :: annualized_return = 0.0_dp
    real(dp) :: annualized_volatility = 0.0_dp
    real(dp) :: skewness = 0.0_dp
    real(dp) :: excess_kurtosis = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: sortino = 0.0_dp
    real(dp) :: maximum_drawdown = 0.0_dp
    real(dp) :: var_value = 0.0_dp
    real(dp) :: es_value = 0.0_dp
    real(dp) :: benchmark_information_ratio = 0.0_dp
  end type performance_summary
contains
  subroutine compute_performance_summary(r,scale,p,result,benchmark,rf,method)
    real(dp),intent(in)::r(:),scale,p
    type(performance_summary),intent(out)::result
    real(dp),intent(in),optional::benchmark(:),rf
    character(len=*),intent(in),optional::method
    real(dp)::risk_free
    character(len=20)::m
    result%observations=size(r);result%mean_return=mean_value(r)
    result%annualized_return=annualized_return(r,scale,.true.)
    result%annualized_volatility=sd_value(r)*sqrt(scale)
    result%skewness=skewness_value(r,1);result%excess_kurtosis=kurtosis_value(r,1,.true.)
    risk_free=0.0_dp;if(present(rf))risk_free=rf
    result%sharpe=(mean_value(r)-risk_free)/max(sd_value(r),tiny(1.0_dp))*sqrt(scale)
    result%sortino=sortino_ratio(r,risk_free)*sqrt(scale)
    result%maximum_drawdown=max_drawdown(r,.true.)
    m='modified';if(present(method))m=method
    result%var_value=value_at_risk(r,p,trim(m));result%es_value=expected_shortfall(r,p,trim(m))
    if(present(benchmark))result%benchmark_information_ratio=information_ratio(r,benchmark,scale)
  end subroutine compute_performance_summary


  subroutine conditional_capm_fit(asset,market,rf,z,lags,result)
    real(dp),intent(in)::asset(:),market(:),rf(:),z(:,:)
    integer,intent(in)::lags
    type(conditional_capm_result),intent(out)::result
    real(dp),allocatable::x(:,:),y(:),beta(:),resid(:),zc(:,:)
    integer::n,q,nobs,t,j,l,col,row
    logical::ok
    n=min(size(asset),min(size(market),min(size(rf),size(z,1))))
    q=size(z,2);nobs=max(0,n-lags)
    allocate(x(nobs,2+2*q*lags),y(nobs),beta(2+2*q*lags),resid(nobs),zc(n,q))
    do j=1,q;zc(:,j)=z(:n,j)-mean_value(z(:n,j));end do
    do row=1,nobs
      t=lags+row;y(row)=asset(t)-rf(t);x(row,1)=1.0_dp;x(row,2)=market(t)-rf(t);col=2
      do l=1,lags;do j=1,q;col=col+1;x(row,col)=zc(t-l,j);end do;end do
      do l=1,lags;do j=1,q;col=col+1;x(row,col)=zc(t-l,j)*(market(t)-rf(t));end do;end do
    end do
    if(nobs>size(beta))then
      call ols_fit(x,y,beta,resid,result%r_squared,ok)
    else
      beta=0.0_dp;resid=0.0_dp;result%r_squared=0.0_dp;ok=.false.
    end if
    result%lags=lags;result%observations=nobs;result%success=ok
    allocate(result%coefficients(size(beta)),result%residuals(size(resid)));result%coefficients=beta;result%residuals=resid
  end subroutine conditional_capm_fit

  subroutine outperformance_probabilities(asset,benchmark,period_lengths,wins,losses,total,probability)
    real(dp),intent(in)::asset(:),benchmark(:)
    integer,intent(in)::period_lengths(:)
    integer,intent(out)::wins(:),losses(:),total(:)
    real(dp),intent(out)::probability(:)
    real(dp)::ra,rb
    integer::n,j,t,h
    n=min(size(asset),size(benchmark));wins=0;losses=0;total=0;probability=0.0_dp
    do j=1,min(size(period_lengths),size(wins),size(losses),size(total),size(probability))
      h=period_lengths(j)
      if(h<=0 .or. h>n)cycle
      do t=h,n
        ra=product(1.0_dp+asset(t-h+1:t))-1.0_dp
        rb=product(1.0_dp+benchmark(t-h+1:t))-1.0_dp
        if(ra>rb)then;wins(j)=wins(j)+1
        else if(ra<rb)then;losses(j)=losses(j)+1
        end if
      end do
      total(j)=wins(j)+losses(j)
      if(total(j)>0)probability(j)=real(wins(j),dp)/real(total(j),dp)
    end do
  end subroutine outperformance_probabilities

  subroutine rolling_sfm(ra,rb,rf,width,alpha,beta,r2,nout)
    real(dp),intent(in)::ra(:),rb(:),rf(:)
    integer,intent(in)::width
    real(dp),allocatable,intent(out)::alpha(:),beta(:),r2(:)
    integer,intent(out)::nout
    type(sfm_result)::fit
    integer::n,i
    n=min(size(ra),min(size(rb),size(rf)));nout=max(0,n-width+1)
    allocate(alpha(nout),beta(nout),r2(nout))
    do i=1,nout
      call sfm_fit(ra(i:i+width-1),rb(i:i+width-1),rf(i:i+width-1),fit)
      alpha(i)=fit%alpha;beta(i)=fit%beta;r2(i)=fit%r_squared
    end do
  end subroutine rolling_sfm

  subroutine expanding_sfm(ra,rb,rf,min_obs,alpha,beta,r2,nout)
    real(dp),intent(in)::ra(:),rb(:),rf(:)
    integer,intent(in)::min_obs
    real(dp),allocatable,intent(out)::alpha(:),beta(:),r2(:)
    integer,intent(out)::nout
    type(sfm_result)::fit
    integer::n,i,endp
    n=min(size(ra),min(size(rb),size(rf)));nout=max(0,n-min_obs+1)
    allocate(alpha(nout),beta(nout),r2(nout))
    do i=1,nout
      endp=min_obs+i-1;call sfm_fit(ra(:endp),rb(:endp),rf(:endp),fit)
      alpha(i)=fit%alpha;beta(i)=fit%beta;r2(i)=fit%r_squared
    end do
  end subroutine expanding_sfm

  subroutine capture_ratios(asset,benchmark,up_capture,down_capture,capture_ratio)
    real(dp),intent(in)::asset(:),benchmark(:)
    real(dp),intent(out)::up_capture,down_capture,capture_ratio
    integer::n,nup,ndown
    real(dp)::aup,bup,adown,bdown
    n=min(size(asset),size(benchmark));nup=count(benchmark(:n)>0.0_dp);ndown=count(benchmark(:n)<0.0_dp)
    if(nup>0)then
      aup=product(1.0_dp+pack(asset(:n),benchmark(:n)>0.0_dp))**(1.0_dp/real(nup,dp))-1.0_dp
      bup=product(1.0_dp+pack(benchmark(:n),benchmark(:n)>0.0_dp))**(1.0_dp/real(nup,dp))-1.0_dp
      if(abs(bup)>tiny(1.0_dp))then;up_capture=aup/bup;else;up_capture=0.0_dp;end if
    else;up_capture=0.0_dp;end if
    if(ndown>0)then
      adown=product(1.0_dp+pack(asset(:n),benchmark(:n)<0.0_dp))**(1.0_dp/real(ndown,dp))-1.0_dp
      bdown=product(1.0_dp+pack(benchmark(:n),benchmark(:n)<0.0_dp))**(1.0_dp/real(ndown,dp))-1.0_dp
      if(abs(bdown)>tiny(1.0_dp))then;down_capture=adown/bdown;else;down_capture=0.0_dp;end if
    else;down_capture=0.0_dp;end if
    if(abs(down_capture)>tiny(1.0_dp))then;capture_ratio=up_capture/down_capture;else;capture_ratio=0.0_dp;end if
  end subroutine capture_ratios
end module analytics_extensions_mod
