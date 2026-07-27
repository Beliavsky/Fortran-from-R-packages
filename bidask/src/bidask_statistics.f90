! SPDX-License-Identifier: MIT
! Based on bidask 2.1.5, Copyright (c) 2024 Emanuele Guidotti.
module bidask_statistics
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use bidask_kinds, only: dp
  implicit none
  private
  public :: nan_dp, mean_values, sum_values, signed_root, safe_log

contains

  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure real(dp) function mean_values(x, na_rm) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: na_rm
    integer :: i, n
    real(dp) :: total

    total = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        total = total + x(i)
        n = n + 1
      else if (.not. na_rm) then
        value = nan_dp()
        return
      end if
    end do
    if (n == 0) then
      value = nan_dp()
    else
      value = total / real(n, dp)
    end if
  end function mean_values

  pure real(dp) function sum_values(x, na_rm) result(value)
    real(dp), intent(in) :: x(:)
    logical, intent(in) :: na_rm
    integer :: i, n

    value = 0.0_dp
    n = 0
    do i = 1, size(x)
      if (ieee_is_finite(x(i))) then
        value = value + x(i)
        n = n + 1
      else if (.not. na_rm) then
        value = nan_dp()
        return
      end if
    end do
    if (n == 0) value = nan_dp()
  end function sum_values

  pure real(dp) function signed_root(x, keep_sign) result(value)
    real(dp), intent(in) :: x
    logical, intent(in) :: keep_sign
    if (.not. ieee_is_finite(x)) then
      value = nan_dp()
    else if (keep_sign .and. x < 0.0_dp) then
      value = -sqrt(abs(x))
    else
      value = sqrt(abs(x))
    end if
  end function signed_root

  pure real(dp) function safe_log(x) result(value)
    real(dp), intent(in) :: x
    if (ieee_is_finite(x) .and. x > 0.0_dp) then
      value = log(x)
    else
      value = nan_dp()
    end if
  end function safe_log

end module bidask_statistics
