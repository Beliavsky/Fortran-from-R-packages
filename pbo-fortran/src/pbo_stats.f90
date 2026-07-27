! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
module pbo_stats
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use pbo_kinds, only : dp
  implicit none
  private
  public :: mean_value, sample_sd, first_argmax, average_rank_of
  public :: fit_line, empirical_cdf, round_significant, all_finite
contains
  pure function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    if (size(x) == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = sum(x) / real(size(x), dp)
    end if
  end function mean_value

  pure function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, xbar
    if (size(x) < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    xbar = mean_value(x)
    value = sqrt(sum((x - xbar)**2) / real(size(x) - 1, dp))
  end function sample_sd

  pure function all_finite(x) result(ok)
    real(dp), intent(in) :: x(:)
    logical :: ok
    ok = all(ieee_is_finite(x))
  end function all_finite

  pure function first_argmax(x) result(index_max)
    real(dp), intent(in) :: x(:)
    integer :: index_max, i
    index_max = 1
    do i = 2, size(x)
      if (x(i) > x(index_max)) index_max = i
    end do
  end function first_argmax

  pure function average_rank_of(x, index_value) result(rank_value)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: index_value
    real(dp) :: rank_value, target
    integer :: n_less, n_equal
    target = x(index_value)
    n_less = count(x < target)
    n_equal = count((.not. (x < target)) .and. (.not. (x > target)))
    rank_value = 1.0_dp + real(n_less, dp) + 0.5_dp * real(n_equal - 1, dp)
  end function average_rank_of

  subroutine fit_line(x, y, intercept, slope, r2, adjusted_r2, success)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: intercept, slope, r2, adjusted_r2
    logical, intent(out) :: success
    real(dp) :: xbar, ybar, ssx, sst, sse
    real(dp), allocatable :: fitted(:)
    integer :: n

    success = .false.
    intercept = ieee_value(0.0_dp, ieee_quiet_nan)
    slope = intercept
    r2 = intercept
    adjusted_r2 = intercept
    n = size(x)
    if (n /= size(y) .or. n < 2) return
    if (.not. all_finite(x) .or. .not. all_finite(y)) return
    xbar = mean_value(x)
    ybar = mean_value(y)
    ssx = sum((x - xbar)**2)
    if (ssx <= epsilon(1.0_dp) * max(1.0_dp, sum(x*x))) return
    slope = sum((x - xbar) * (y - ybar)) / ssx
    intercept = ybar - slope * xbar
    allocate(fitted(n))
    fitted = intercept + slope * x
    sse = sum((y - fitted)**2)
    sst = sum((y - ybar)**2)
    if (sst <= epsilon(1.0_dp) * max(1.0_dp, sum(y*y))) then
      r2 = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      r2 = 1.0_dp - sse / sst
    end if
    if (n > 2 .and. ieee_is_finite(r2)) then
      adjusted_r2 = 1.0_dp - (1.0_dp - r2) * real(n - 1, dp) / real(n - 2, dp)
    end if
    success = .true.
  end subroutine fit_line

  pure function empirical_cdf(sample, x) result(value)
    real(dp), intent(in) :: sample(:), x
    real(dp) :: value
    if (size(sample) == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = real(count(sample <= x), dp) / real(size(sample), dp)
    end if
  end function empirical_cdf

  pure function round_significant(x, digits) result(value)
    real(dp), intent(in) :: x
    integer, intent(in) :: digits
    real(dp) :: value, scale
    integer :: exponent10
    if (abs(x) < tiny(1.0_dp) .or. .not. ieee_is_finite(x)) then
      value = x
      return
    end if
    exponent10 = floor(log10(abs(x)))
    scale = 10.0_dp**real(digits - 1 - exponent10, dp)
    value = anint(x * scale) / scale
  end function round_significant
end module pbo_stats
