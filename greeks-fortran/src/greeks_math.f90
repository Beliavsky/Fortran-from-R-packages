! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_math
  use greeks_kinds, only: dp, pi
  implicit none
  private
  public :: normal_pdf, normal_cdf, sample_mean, sample_se
  public :: control_variate_intercept, control_variate_stats, safe_divide
contains
  pure elemental function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure elemental function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    real(dp) :: value
    value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure function sample_mean(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    if (size(x) == 0) then
      value = 0.0_dp
    else
      value = sum(x)/real(size(x), dp)
    end if
  end function sample_mean

  pure function sample_se(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value, avg
    integer :: n
    n = size(x)
    if (n <= 1) then
      value = 0.0_dp
    else
      avg = sum(x)/real(n, dp)
      value = sqrt(sum((x-avg)**2)/real(n-1, dp)/real(n, dp))
    end if
  end function sample_se

  pure subroutine control_variate_stats(y, x, estimate, se)
    real(dp), intent(in) :: y(:), x(:)
    real(dp), intent(out) :: estimate, se
    real(dp) :: mx, my, denom, beta
    real(dp), allocatable :: adjusted(:)
    integer :: n
    n = min(size(y), size(x))
    if (n <= 0) then
      estimate = 0.0_dp
      se = 0.0_dp
      return
    end if
    mx = sum(x(1:n))/real(n, dp)
    my = sum(y(1:n))/real(n, dp)
    denom = sum((x(1:n)-mx)**2)
    if (denom <= epsilon(1.0_dp)*max(1.0_dp, sum(x(1:n)**2))) then
      beta = 0.0_dp
    else
      beta = sum((x(1:n)-mx)*(y(1:n)-my))/denom
    end if
    allocate(adjusted(n))
    adjusted = y(1:n)-beta*x(1:n)
    estimate = sum(adjusted)/real(n,dp)
    se = sample_se(adjusted)
  end subroutine control_variate_stats

  pure function control_variate_intercept(y, x) result(value)
    real(dp), intent(in) :: y(:), x(:)
    real(dp) :: value, mx, my, denom, beta
    integer :: n
    n = min(size(y), size(x))
    if (n <= 1) then
      value = 0.0_dp
      if (n == 1) value = y(1)
      return
    end if
    mx = sum(x(1:n))/real(n, dp)
    my = sum(y(1:n))/real(n, dp)
    denom = sum((x(1:n)-mx)**2)
    if (denom <= epsilon(1.0_dp)*max(1.0_dp, sum(x(1:n)**2))) then
      value = my
    else
      beta = sum((x(1:n)-mx)*(y(1:n)-my))/denom
      value = my - beta*mx
    end if
  end function control_variate_intercept

  pure elemental function safe_divide(numerator, denominator, fallback) result(value)
    real(dp), intent(in) :: numerator, denominator, fallback
    real(dp) :: value
    if (abs(denominator) <= tiny(1.0_dp)) then
      value = fallback
    else
      value = numerator/denominator
    end if
  end function safe_divide
end module greeks_math
