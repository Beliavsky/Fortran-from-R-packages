! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency
  use highfrequency_kinds, only: dp, pi
  use highfrequency_types
  use highfrequency_stats, only: normal_cdf, normal_quantile, normal_pdf
  use highfrequency_linalg, only: make_psd, covariance_to_correlation
  use highfrequency_data, only: make_returns, aggregate_last, aggregate_sum
  use highfrequency_data, only: previous_tick, refresh_time_pair, trade_direction
  use highfrequency_data, only: liquidity_measures, merge_same_timestamp, make_ohlcv
  use highfrequency_realized
  use highfrequency_models
  use highfrequency_jumps
  use highfrequency_leadlag
  use highfrequency_spot
  use highfrequency_remedi
  use highfrequency_cleaning
  implicit none
  private

  public :: dp, pi
  public :: jump_test_result, iv_inference_result, har_model, heavy_model
  public :: liquidity_result, lead_lag_result
  public :: make_returns, aggregate_last, aggregate_sum, previous_tick
  public :: refresh_time_pair, trade_direction, liquidity_measures
  public :: merge_same_timestamp, make_ohlcv, make_psd
  public :: realized_variance, realized_covariance, realized_skewness
  public :: realized_kurtosis, realized_semivariance, realized_semicovariance
  public :: bipower_variation, bipower_covariance, realized_quarticity
  public :: tripower_quarticity, quadpower_variation, minimum_realized_variance
  public :: median_realized_variance, minimum_realized_quarticity
  public :: median_realized_quarticity, multipower_variation
  public :: threshold_covariance, realized_kernel_covariance
  public :: two_scale_variance, two_scale_covariance
  public :: hayashi_yoshida_covariance, realized_beta, noise_variance
  public :: modulated_realized_covariance
  public :: preaveraged_covariance, average_realized_covariance
  public :: fit_har, har_forecast, fit_heavy, heavy_forecast, heavy_recursion
  public :: bns_jump_test, aj_jump_test, abd_jump_test, iv_inference
  public :: lead_lag, spot_volatility, spot_drift, drift_burst_statistic
  public :: preaverage_returns, remedi, choose_remedi_kn
  public :: no_zero_prices_mask, no_zero_quotes_mask
  public :: nonnegative_spread_mask, maximum_spread_mask, price_outlier_mask
  public :: match_trades_quotes, spread_prices, business_time_groups
  public :: normal_cdf, normal_quantile, normal_pdf
  public :: rrvar, rcov, rskew, rkurt, rsvar, rsemicov
  public :: rbpvar, rbpcov, rquar, rtpquar, rqpvar
  public :: rminrvar, rmedrvar, rminrquar, rmedrquar, rmpvar
  public :: rthresholdcov, rkernelcov, rtsvar, rtscov, rhy_cov, rbeta

  interface rrvar
    module procedure rrvar_vector
    module procedure rrvar_matrix
  end interface rrvar

contains

  pure real(dp) function rrvar_vector(r) result(value)
    real(dp),intent(in)::r(:)
    value=realized_variance(r)
  end function rrvar_vector

  pure function rrvar_matrix(r) result(value)
    real(dp),intent(in)::r(:,:)
    real(dp)::value(size(r,2))
    value=realized_variance(r)
  end function rrvar_matrix

  pure function rcov(r,correlation) result(value)
    real(dp),intent(in)::r(:,:)
    logical,intent(in),optional::correlation
    real(dp)::value(size(r,2),size(r,2))
    value=realized_covariance(r,correlation)
  end function rcov

  pure real(dp) function rskew(r) result(value)
    real(dp),intent(in)::r(:)
    value=realized_skewness(r)
  end function rskew

  pure real(dp) function rkurt(r) result(value)
    real(dp),intent(in)::r(:)
    value=realized_kurtosis(r)
  end function rkurt

  pure function rsvar(r) result(value)
    real(dp),intent(in)::r(:)
    real(dp)::value(2)
    value=realized_semivariance(r)
  end function rsvar

  pure function rsemicov(r,side) result(value)
    real(dp),intent(in)::r(:,:)
    integer,intent(in),optional::side
    real(dp)::value(size(r,2),size(r,2))
    value=realized_semicovariance(r,side)
  end function rsemicov

  pure real(dp) function rbpvar(r) result(value)
    real(dp),intent(in)::r(:)
    value=bipower_variation(r)
  end function rbpvar

  pure function rbpcov(r,correlation) result(value)
    real(dp),intent(in)::r(:,:)
    logical,intent(in),optional::correlation
    real(dp)::value(size(r,2),size(r,2))
    value=bipower_covariance(r,correlation)
  end function rbpcov

  pure real(dp) function rquar(r) result(value)
    real(dp),intent(in)::r(:)
    value=realized_quarticity(r)
  end function rquar

  pure real(dp) function rtpquar(r) result(value)
    real(dp),intent(in)::r(:)
    value=tripower_quarticity(r)
  end function rtpquar

  pure real(dp) function rqpvar(r) result(value)
    real(dp),intent(in)::r(:)
    value=quadpower_variation(r)
  end function rqpvar

  pure real(dp) function rminrvar(r) result(value)
    real(dp),intent(in)::r(:)
    value=minimum_realized_variance(r)
  end function rminrvar

  pure real(dp) function rmedrvar(r) result(value)
    real(dp),intent(in)::r(:)
    value=median_realized_variance(r)
  end function rmedrvar

  pure real(dp) function rminrquar(r) result(value)
    real(dp),intent(in)::r(:)
    value=minimum_realized_quarticity(r)
  end function rminrquar

  pure real(dp) function rmedrquar(r) result(value)
    real(dp),intent(in)::r(:)
    value=median_realized_quarticity(r)
  end function rmedrquar

  pure real(dp) function rmpvar(r,m,p) result(value)
    real(dp),intent(in)::r(:),p
    integer,intent(in)::m
    value=multipower_variation(r,m,p)
  end function rmpvar

  function rthresholdcov(r,correlation) result(value)
    real(dp),intent(in)::r(:,:)
    logical,intent(in),optional::correlation
    real(dp)::value(size(r,2),size(r,2))
    value=threshold_covariance(r,correlation)
  end function rthresholdcov

  function rkernelcov(r,bandwidth,kernel,correlation,force_psd) result(value)
    real(dp),intent(in)::r(:,:)
    integer,intent(in)::bandwidth
    character(len=*),intent(in),optional::kernel
    logical,intent(in),optional::correlation,force_psd
    real(dp)::value(size(r,2),size(r,2))
    value=realized_kernel_covariance(r,bandwidth,kernel,correlation,force_psd)
  end function rkernelcov

  pure real(dp) function rtsvar(prices,k,j) result(value)
    real(dp),intent(in)::prices(:)
    integer,intent(in)::k,j
    value=two_scale_variance(prices,k,j)
  end function rtsvar

  pure real(dp) function rtscov(px,py,k,j) result(value)
    real(dp),intent(in)::px(:),py(:)
    integer,intent(in)::k,j
    value=two_scale_covariance(px,py,k,j)
  end function rtscov

  pure real(dp) function rhy_cov(t1,p1,t2,p2) result(value)
    integer,intent(in)::t1(:),t2(:)
    real(dp),intent(in)::p1(:),p2(:)
    value=hayashi_yoshida_covariance(t1,p1,t2,p2)
  end function rhy_cov

  pure real(dp) function rbeta(asset,index_returns) result(value)
    real(dp),intent(in)::asset(:),index_returns(:)
    value=realized_beta(asset,index_returns)
  end function rbeta

end module highfrequency
