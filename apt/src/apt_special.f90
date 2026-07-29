! SPDX-License-Identifier: GPL-2.0-or-later
module apt_special
  use apt_kinds, only : dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: normal_cdf, student_t_cdf, f_cdf, chi_square_cdf

contains

  pure real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function beta_cf(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: max_iter = 400
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    integer :: m, m2
    real(dp) :: aa, c, d, del, h, qab, qam, qap

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, max_iter
      m2 = 2 * m
      aa = real(m, dp) * (b - real(m, dp)) * x / &
        ((qam + real(m2, dp)) * (a + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
        ((a + real(m2, dp)) * (qap + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_cf

  pure real(dp) function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt, xx

    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    xx = min(1.0_dp, max(0.0_dp, x))
    if (xx <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (xx >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
      a * log(xx) + b * log(1.0_dp - xx))
    if (xx < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      p = bt * beta_cf(a, b, xx) / a
    else
      p = 1.0_dp - bt * beta_cf(b, a, 1.0_dp - xx) / b
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function regularized_beta

  pure real(dp) function student_t_cdf(x, nu) result(p)
    real(dp), intent(in) :: x, nu
    real(dp) :: z, ib
    if (nu <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      p = 0.5_dp
      return
    end if
    z = nu / (nu + x * x)
    ib = regularized_beta(z, 0.5_dp * nu, 0.5_dp)
    if (x > 0.0_dp) then
      p = 1.0_dp - 0.5_dp * ib
    else
      p = 0.5_dp * ib
    end if
  end function student_t_cdf

  pure real(dp) function f_cdf(x, d1, d2) result(p)
    real(dp), intent(in) :: x, d1, d2
    real(dp) :: z
    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (d1 <= 0.0_dp .or. d2 <= 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      z = d1 * x / (d1 * x + d2)
      p = regularized_beta(z, 0.5_dp * d1, 0.5_dp * d2)
    end if
  end function f_cdf

  pure real(dp) function regularized_gamma_p(a, x) result(p)
    real(dp), intent(in) :: a, x
    integer, parameter :: max_iter = 1000
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    integer :: n, i
    real(dp) :: ap, del, sumv, b, c, d, h, an

    if (a <= 0.0_dp .or. x < 0.0_dp) then
      p = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    else if (abs(x) <= tiny(1.0_dp)) then
      p = 0.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp / a
      del = sumv
      do n = 1, max_iter
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) < abs(sumv) * eps) exit
      end do
      p = sumv * exp(-x + a * log(x) - log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp / fpmin
      d = 1.0_dp / b
      h = d
      do i = 1, max_iter
        an = -real(i, dp) * (real(i, dp) - a)
        b = b + 2.0_dp
        d = an * d + b
        if (abs(d) < fpmin) d = fpmin
        c = b + an / c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp / d
        del = d * c
        h = h * del
        if (abs(del - 1.0_dp) < eps) exit
      end do
      p = 1.0_dp - exp(-x + a * log(x) - log_gamma(a)) * h
    end if
    p = min(1.0_dp, max(0.0_dp, p))
  end function regularized_gamma_p

  pure real(dp) function chi_square_cdf(x, df) result(p)
    real(dp), intent(in) :: x, df
    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = regularized_gamma_p(0.5_dp * df, 0.5_dp * x)
    end if
  end function chi_square_cdf

end module apt_special
