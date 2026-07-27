! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo_metrics
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_quiet_nan
  use pbo_kinds, only : dp
  use pbo_stats, only : mean_value, sample_sd
  implicit none
  private
  public :: column_mean, column_sum, sharpe_ratio, sharpe_ratio_rf
  public :: omega_ratio, omega_ratio_threshold
contains
  subroutine column_mean(data, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: values(:)
    integer :: j
    do j = 1, size(data,2)
      values(j) = mean_value(data(:,j))
    end do
  end subroutine column_mean

  subroutine column_sum(data, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: values(:)
    integer :: j
    do j = 1, size(data,2)
      values(j) = sum(data(:,j))
    end do
  end subroutine column_sum

  subroutine sharpe_ratio(data, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: values(:)
    call sharpe_ratio_rf(data, 0.0_dp, values)
  end subroutine sharpe_ratio

  subroutine sharpe_ratio_rf(data, risk_free, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(in) :: risk_free
    real(dp), intent(out) :: values(:)
    real(dp) :: sd_value
    integer :: j
    do j = 1, size(data,2)
      sd_value = sample_sd(data(:,j) - risk_free)
      if (sd_value > 0.0_dp) then
        values(j) = mean_value(data(:,j) - risk_free) / sd_value
      else
        values(j) = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
    end do
  end subroutine sharpe_ratio_rf

  subroutine omega_ratio(data, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(out) :: values(:)
    call omega_ratio_threshold(data, 0.0_dp, values)
  end subroutine omega_ratio

  subroutine omega_ratio_threshold(data, threshold, values)
    real(dp), intent(in) :: data(:,:)
    real(dp), intent(in) :: threshold
    real(dp), intent(out) :: values(:)
    real(dp) :: gains, losses
    integer :: j
    do j = 1, size(data,2)
      gains = sum(max(data(:,j) - threshold, 0.0_dp))
      losses = sum(max(threshold - data(:,j), 0.0_dp))
      if (losses > 0.0_dp) then
        values(j) = gains / losses
      else if (gains > 0.0_dp) then
        values(j) = ieee_value(0.0_dp, ieee_positive_inf)
      else
        values(j) = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
    end do
  end subroutine omega_ratio_threshold
end module pbo_metrics
