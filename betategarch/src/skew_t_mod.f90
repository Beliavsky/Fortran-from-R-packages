! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2011-2025 Genaro Sucarrat
! Copyright (C) 2026 contributors to the Modern Fortran translation
!
! This file is part of betategarch-modern-fortran.
! It is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License version 2 only.

module skew_t_mod
  use betategarch_kinds, only : dp
  use betategarch_math, only : beta_fn, log_beta_fn, pi, signum
  use betategarch_rng, only : random_student_t
  implicit none
  private

  public :: skew_t_random, skew_t_pdf, skew_t_logpdf
  public :: skew_t_mean, skew_t_variance, skew_t_skewness, skew_t_kurtosis
  public :: skew_t_raw_moment

contains

  subroutine skew_t_random(x, df, skew)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in) :: df, skew

    real(dp) :: zstar, weight, z, u
    integer :: i, s

    if (df <= 0.0_dp .or. skew <= 0.0_dp) error stop "skew_t_random: invalid df or skew"
    weight = skew/(skew + 1.0_dp/skew)
    do i = 1, size(x)
      zstar = random_student_t(df)
      call random_number(u)
      z = -weight + u
      s = signum(z)
      x(i) = -abs(zstar)/(skew**s)*real(s, dp)
    end do
  end subroutine skew_t_random

  pure elemental function skew_t_logpdf(y, df, scale, skew) result(value)
    real(dp), intent(in) :: y, df, scale, skew
    real(dp) :: value

    integer :: s
    real(dp) :: log_numerator, log_denom1, log_denom2

    if (df <= 0.0_dp .or. scale <= 0.0_dp .or. skew <= 0.0_dp) then
      value = -huge(1.0_dp)
      return
    end if
    s = signum(y)
    log_numerator = log(2.0_dp) - log(skew + 1.0_dp/skew)
    log_denom1 = log_beta_fn(0.5_dp, 0.5_dp*df) + log(scale) + 0.5_dp*log(df)
    log_denom2 = 0.5_dp*(df + 1.0_dp)*log(1.0_dp + y*y/(skew**(2*s)*df*scale*scale))
    value = log_numerator - log_denom1 - log_denom2
  end function skew_t_logpdf

  pure elemental function skew_t_pdf(y, df, scale, skew) result(value)
    real(dp), intent(in) :: y, df, scale, skew
    real(dp) :: value

    value = exp(skew_t_logpdf(y, df, scale, skew))
  end function skew_t_pdf

  pure function skew_t_raw_moment(order, df, skew) result(value)
    integer, intent(in) :: order
    real(dp), intent(in) :: df, skew
    real(dp) :: value

    real(dp) :: abs_moment, sign_factor

    if (order < 0 .or. df <= real(order, dp) .or. skew <= 0.0_dp) then
      value = huge(1.0_dp)
      return
    end if
    if (order == 0) then
      value = 1.0_dp
      return
    end if

    abs_moment = df**(0.5_dp*order) * exp(log_gamma(0.5_dp*(order + 1)) + &
      log_gamma(0.5_dp*(df - order)) - log_gamma(0.5_dp) - log_gamma(0.5_dp*df))
    if (mod(order, 2) == 0) then
      sign_factor = skew**(order + 1) + skew**(-(order + 1))
    else
      sign_factor = skew**(order + 1) - skew**(-(order + 1))
    end if
    value = abs_moment*sign_factor/(skew + 1.0_dp/skew)
  end function skew_t_raw_moment

  pure function skew_t_mean(df, skew) result(value)
    real(dp), intent(in) :: df, skew
    real(dp) :: value

    value = (skew - 1.0_dp/skew)*sqrt(df)*beta_fn(0.5_dp*(df - 1.0_dp), 0.5_dp)/pi
  end function skew_t_mean

  pure function skew_t_variance(df, skew) result(value)
    real(dp), intent(in) :: df, skew
    real(dp) :: value

    real(dp) :: mu, raw2

    mu = skew_t_mean(df, skew)
    raw2 = skew_t_raw_moment(2, df, skew)
    value = raw2 - mu*mu
  end function skew_t_variance

  pure function skew_t_skewness(df, skew) result(value)
    real(dp), intent(in) :: df, skew
    real(dp) :: value

    real(dp) :: mu, variance, raw2, raw3

    mu = skew_t_mean(df, skew)
    variance = skew_t_variance(df, skew)
    raw2 = skew_t_raw_moment(2, df, skew)
    raw3 = skew_t_raw_moment(3, df, skew)
    value = (raw3 - 3.0_dp*mu*raw2 + 2.0_dp*mu**3)/variance**1.5_dp
  end function skew_t_skewness

  pure function skew_t_kurtosis(df, skew) result(value)
    real(dp), intent(in) :: df, skew
    real(dp) :: value

    real(dp) :: mu, variance, raw2, raw3, raw4

    mu = skew_t_mean(df, skew)
    variance = skew_t_variance(df, skew)
    raw2 = skew_t_raw_moment(2, df, skew)
    raw3 = skew_t_raw_moment(3, df, skew)
    raw4 = skew_t_raw_moment(4, df, skew)
    value = (raw4 - 4.0_dp*mu*raw3 + 6.0_dp*mu*mu*raw2 - 3.0_dp*mu**4)/(variance*variance)
  end function skew_t_kurtosis

end module skew_t_mod
