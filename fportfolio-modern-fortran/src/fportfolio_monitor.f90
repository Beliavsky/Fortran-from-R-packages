! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fPortfolio contributors and modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it under GPL version 2 or later.
module fportfolio_monitor
  use fportfolio_kinds, only: dp
  use fportfolio_risk, only: drawdown_series
  use fportfolio_statistics, only: column_means, sample_covariance
  implicit none
  private
  public :: ema_indicator, macd_indicator, drawdown_indicator, rebalancing_statistics, &
            turning_points, rolling_stability
contains
  function ema_indicator(series,lambda,initial_count) result(out)
    real(dp),intent(in)::series(:),lambda
    integer,intent(in),optional::initial_count
    real(dp)::out(size(series))
    integer::i,k
    k=min(size(series),10);if(present(initial_count))k=max(1,min(size(series),initial_count))
    out(1)=sum(series(1:k))/real(k,dp)
    do i=2,size(series);out(i)=(1.0_dp-lambda)*series(i)+lambda*out(i-1);end do
  end function ema_indicator

  subroutine macd_indicator(index,lambda_fast,lambda_slow,lambda_signal,macd,signal,histogram,position)
    real(dp),intent(in)::index(:),lambda_fast,lambda_slow,lambda_signal
    real(dp),allocatable,intent(out)::macd(:),signal(:),histogram(:),position(:)
    real(dp),allocatable::ef(:),es(:),log_index(:)
    integer::n
    n=size(index);allocate(macd(n),signal(n),histogram(n),position(n),ef(n),es(n),log_index(n))
    log_index=log(max(index,tiny(1.0_dp)));ef=ema_indicator(log_index,lambda_fast);es=ema_indicator(log_index,lambda_slow)
    macd=ef-es;signal=ema_indicator(macd,lambda_signal);histogram=macd-signal;position=sign(1.0_dp,histogram)
  end subroutine macd_indicator

  subroutine drawdown_indicator(returns,lambda_fast,lambda_slow,lambda_signal,drawdown,macd,signal,histogram,position)
    real(dp),intent(in)::returns(:),lambda_fast,lambda_slow,lambda_signal
    real(dp),allocatable,intent(out)::drawdown(:),macd(:),signal(:),histogram(:),position(:)
    real(dp),allocatable::ef(:),es(:)
    integer::n
    n=size(returns);allocate(drawdown(n),macd(n),signal(n),histogram(n),position(n),ef(n),es(n))
    drawdown=drawdown_series(returns);ef=ema_indicator(drawdown,lambda_fast);es=ema_indicator(drawdown,lambda_slow)
    macd=ef-es;signal=ema_indicator(macd,lambda_signal);histogram=macd-signal
    position=1.0_dp
    where(macd<0.0_dp .and. histogram<0.0_dp)position=0.0_dp
  end subroutine drawdown_indicator

  pure function rebalancing_statistics(returns,positions) result(stats)
    real(dp),intent(in)::returns(:),positions(:)
    real(dp)::stats(4),forecast(size(positions))
    forecast=0.0_dp
    if(size(positions)>1)forecast(2:)=positions(:size(positions)-1)
    stats(1)=sum(abs(returns));stats(2)=sum(returns*positions);stats(3)=sum(returns*forecast);stats(4)=sum(returns)
  end function rebalancing_statistics

  subroutine turning_points(series,smoothing,indices,directions)
    real(dp),intent(in)::series(:),smoothing
    integer,allocatable,intent(out)::indices(:),directions(:)
    real(dp),allocatable::smooth(:)
    integer,allocatable::tmpi(:),tmpd(:)
    integer::i,n,k
    n=size(series);allocate(smooth(n),tmpi(n),tmpd(n));smooth=ema_indicator(series,max(0.0_dp,min(0.999_dp,smoothing)))
    k=0
    do i=2,n-1
      if(smooth(i)>smooth(i-1) .and. smooth(i)>=smooth(i+1))then;k=k+1;tmpi(k)=i;tmpd(k)=-1
      else if(smooth(i)<smooth(i-1) .and. smooth(i)<=smooth(i+1))then;k=k+1;tmpi(k)=i;tmpd(k)=1
      end if
    end do
    allocate(indices(k),directions(k));if(k>0)then;indices=tmpi(:k);directions=tmpd(:k);end if
  end subroutine turning_points

  subroutine rolling_stability(data,width,mean_distance,cov_distance)
    real(dp),intent(in)::data(:,:)
    integer,intent(in)::width
    real(dp),allocatable,intent(out)::mean_distance(:),cov_distance(:)
    real(dp)::base_mu(size(data,2)),base_cov(size(data,2),size(data,2))
    real(dp)::mu(size(data,2)),cov(size(data,2),size(data,2))
    integer::t,n
    n=size(data,1);allocate(mean_distance(n),cov_distance(n));mean_distance=0.0_dp;cov_distance=0.0_dp
    base_mu=column_means(data(:width,:));base_cov=sample_covariance(data(:width,:))
    do t=width,n
      mu=column_means(data(t-width+1:t,:));cov=sample_covariance(data(t-width+1:t,:))
      mean_distance(t)=sqrt(sum((mu-base_mu)**2));cov_distance(t)=sqrt(sum((cov-base_cov)**2))
    end do
  end subroutine rolling_stability
end module fportfolio_monitor
