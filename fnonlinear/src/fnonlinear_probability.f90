! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2000 Adrian Trapletti
! Copyright (C) 2005-2026 Rmetrics contributors
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is part of a modern Fortran translation of fNonlinear and is
! distributed under the GNU General Public License version 2 or later.
module fnonlinear_probability
  use chaos_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: normal_cdf, normal_survival, chi_square_survival, f_survival
contains
  pure real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_survival(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp * erfc(x / sqrt(2.0_dp))
  end function normal_survival

  real(dp) function chi_square_survival(x, degrees_freedom) result(value)
    real(dp), intent(in) :: x, degrees_freedom
    if (x < 0.0_dp .or. degrees_freedom <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    value = regularized_gamma_q(0.5_dp * degrees_freedom, 0.5_dp * x)
  end function chi_square_survival

  real(dp) function f_survival(x, df1, df2) result(value)
    real(dp), intent(in) :: x, df1, df2
    real(dp) :: z
    if (x < 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    z = df2 / (df2 + df1 * x)
    value = regularized_beta(z, 0.5_dp * df2, 0.5_dp * df1)
  end function f_survival

  real(dp) function regularized_gamma_q(a, x) result(q)
    real(dp), intent(in) :: a, x
    real(dp), parameter :: eps = 3.0e-14_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    integer, parameter :: max_iter = 10000
    real(dp) :: ap, sum_value, delta, b, c, d, h, an, p
    integer :: n

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      q = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      q = 1.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      ap = a
      sum_value = 1.0_dp / a
      delta = sum_value
      do n = 1, max_iter
        ap = ap + 1.0_dp
        delta = delta * x / ap
        sum_value = sum_value + delta
        if (abs(delta) <= abs(sum_value) * eps) exit
      end do
      p = sum_value * exp(-x + a * log(x) - log_gamma(a))
      q = max(0.0_dp, min(1.0_dp, 1.0_dp - p))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / max(abs(b), fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do n = 1, max_iter
        an = -real(n, dp) * (real(n, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < fpmin) d = sign(fpmin, d)
        c = b + an / c
        if (abs(c) < fpmin) c = sign(fpmin, c)
        d = 1.0_dp / d
        delta = d * c
        h = h * delta
        if (abs(delta - 1.0_dp) <= eps) exit
      end do
      q = exp(-x + a * log(x) - log_gamma(a)) * h
      q = max(0.0_dp, min(1.0_dp, q))
    end if
  end function regularized_gamma_q

  real(dp) function regularized_beta(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    real(dp) :: front
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      value = 1.0_dp
      return
    end if
    front = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
      a * log(x) + b * log(1.0_dp - x))
    if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      value = front * beta_continued_fraction(a, b, x) / a
    else
      value = 1.0_dp - front * beta_continued_fraction(b, a, 1.0_dp - x) / b
    end if
    value = max(0.0_dp, min(1.0_dp, value))
  end function regularized_beta

  real(dp) function beta_continued_fraction(a, b, x) result(h)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: max_iter = 10000
    real(dp), parameter :: eps = 3.0e-14_dp
    real(dp), parameter :: fpmin = 1.0e-300_dp
    real(dp) :: qab, qap, qam, c, d, aa, delta
    integer :: m, m2

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = sign(fpmin, d)
    d = 1.0_dp / d
    h = d
    do m = 1, max_iter
      m2 = 2 * m
      aa = real(m, dp) * (b - real(m, dp)) * x / &
        ((qam + real(m2, dp)) * (a + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = sign(fpmin, d)
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = sign(fpmin, c)
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
        ((a + real(m2, dp)) * (qap + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = sign(fpmin, d)
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = sign(fpmin, c)
      d = 1.0_dp / d
      delta = d * c
      h = h * delta
      if (abs(delta - 1.0_dp) <= eps) exit
    end do
  end function beta_continued_fraction
end module fnonlinear_probability
