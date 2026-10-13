! SPDX-License-Identifier: GPL-2.0-only
module propcis_common
  !! Shared result type and domain validation for PropCIs.
  use r_kinds, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
  implicit none
  private
  public :: dp, ci_result, valid_counts, valid_level, valid_paired, nan_value, infinity
  integer, parameter, public :: ci_ok = 0, ci_invalid = 1, ci_numerical = 2
  type, public :: ci_result
    !! Bounds and optional point estimate; status 0 success, 1 invalid input, 2 numerical failure.
    real(dp) :: lower = 0.0_dp
    real(dp) :: upper = 0.0_dp
    real(dp) :: estimate = 0.0_dp
    real(dp) :: conf_level = 0.95_dp
    integer :: status = ci_ok
  end type ci_result
contains
  pure elemental logical function valid_counts(x, n) result(ok)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Strictly positive number of trials.
    ok = n > 0 .and. x >= 0 .and. x <= n
  end function valid_counts
  pure elemental logical function valid_paired(b, c, n) result(ok)
    integer, intent(in) :: b !! Discordant pairs of the first kind, nonnegative.
    integer, intent(in) :: c !! Discordant pairs of the second kind, nonnegative.
    integer, intent(in) :: n !! Total pairs, positive and at least b plus c.
    ok = n > 0 .and. b >= 0 .and. c >= 0
    if (ok) ok = b <= n .and. c <= n - b
  end function valid_paired
  pure elemental logical function valid_level(level) result(ok)
    real(dp), intent(in) :: level !! Confidence or posterior probability, strictly between zero and one.
    ok = ieee_is_finite(level)
    if (ok) ok = level > 0.0_dp .and. level < 1.0_dp
  end function valid_level
  pure real(dp) function nan_value() result(value)
    value = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_value
  pure real(dp) function infinity() result(value)
    value = ieee_value(0.0_dp, ieee_positive_inf)
  end function infinity
end module propcis_common
