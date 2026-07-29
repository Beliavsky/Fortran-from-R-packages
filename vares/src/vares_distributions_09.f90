! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_09
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: ddweibull, pdweibull, vardweibull, esdweibull, dexpweibull, pexpweibull
  public :: varexpweibull, esexpweibull, dgenpowerweibull, pgenpowerweibull, vargenpowerweibull, esgenpowerweibull
  public :: dchen, pchen, varchen, eschen, dxie, pxie
  public :: varxie, esxie, vartl, estl, varrs, esrs
  public :: varfr, esfr, varhl, eshl, dloglog, ploglog
  public :: varloglog, esloglog, dexplog, pexplog, varexplog, esexplog
  public :: dexpgeo, pexpgeo, varexpgeo, esexpgeo
contains
  pure elemental function ddweibull(x, c, mu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: c, mu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: c_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    c_v = 1.0_dp
    if (present(c)) c_v = c
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (0.5_dp * weibull_pdf((x - mu_v), shape=c_v, scale=sigma_v))
    if (((log_pdf_v))) res = (weibull_pdf((x - mu_v), shape=c_v, scale=sigma_v, log_pdf=.true.) - log(2.0_dp))
  end function ddweibull

  pure elemental function pdweibull(x, c, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: c, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: c_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    c_v = 1.0_dp
    if (present(c)) c_v = c
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= mu_v)))) res = (0.5_dp * exp(((-((mu_v - x) ** c_v)) &
      & / (sigma_v ** c_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > mu_v)))) res = (1.0_dp - (0.5_dp * exp(((-((x - mu_v) &
      & ** c_v)) / (sigma_v ** c_v)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= mu_v)))) res = (1.0_dp - (0.5_dp * &
      & exp(((-((mu_v - x) ** c_v)) / (sigma_v ** c_v)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > mu_v)))) res = (0.5_dp * exp(((-((x - mu_v) ** &
      & c_v)) / (sigma_v ** c_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= mu_v)))) res = (((-((mu_v - x) ** c_v)) / (sigma_v ** &
      & c_v)) - log(2.0_dp))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > mu_v)))) res = log((1.0_dp - (0.5_dp * exp(((-((x - mu_v) &
      & ** c_v)) / (sigma_v ** c_v))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= mu_v)))) res = log((1.0_dp - (0.5_dp * exp(((-((mu_v &
      & - x) ** c_v)) / (sigma_v ** c_v))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > mu_v)))) res = (((-((x - mu_v) ** c_v)) / (sigma_v ** &
      & c_v)) - log(2.0_dp))
  end function pdweibull

  pure elemental function vardweibull(p, c, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: c, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: c_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    c_v = 1.0_dp
    if (present(c)) c_v = c
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= 0.5_dp))) res = (mu_v - (sigma_v * ((-log((2.0_dp * pp))) ** (1.0_dp / c_v))))
    if (((pp > 0.5_dp))) res = (mu_v + (sigma_v * ((-log((2.0_dp * (1.0_dp - pp)))) ** (1.0_dp / c_v))))
  end function vardweibull

  pure elemental function esdweibull(p, c, mu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: c, mu, sigma
    real(dp) :: res
    real(dp) :: c_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    c_v = 1.0_dp
    if (present(c)) c_v = c
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vardweibull(p * s_quad**3, c_v, mu_v, sigma_v)
    end do
  end function esdweibull

  pure elemental function dexpweibull(x, a, alpha, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, alpha, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((a_v * weibull_pdf(x, shape=alpha_v, scale=sigma_v)) * (weibull_cdf(x, &
      & shape=alpha_v, scale=sigma_v) ** (a_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((log(a_v) + weibull_pdf(x, shape=alpha_v, scale=sigma_v, log_pdf=.true.)) + ((a_v - &
      & 1.0_dp) * log(weibull_cdf(x, shape=alpha_v, scale=sigma_v))))
  end function dexpweibull

  pure elemental function pexpweibull(x, a, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (weibull_cdf(x, shape=alpha_v, scale=sigma_v) ** a_v)
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (a_v * log(weibull_cdf(x, shape=alpha_v, scale=sigma_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v)))
  end function pexpweibull

  pure elemental function varexpweibull(p, a, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (sigma_v * ((-log((1.0_dp - (pp ** (1.0_dp / a_v))))) ** (1.0_dp / alpha_v)))
  end function varexpweibull

  pure elemental function esexpweibull(p, a, alpha, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, alpha, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexpweibull(p * s_quad**3, a_v, alpha_v, sigma_v)
    end do
  end function esexpweibull

  pure elemental function dgenpowerweibull(x, a, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * theta_v) * (x ** (a_v - 1.0_dp))) * ((1.0_dp + (x ** a_v)) ** (theta_v &
      & - 1.0_dp))) * exp((1.0_dp - ((1.0_dp + (x ** a_v)) ** theta_v))))
    if (((log_pdf_v))) res = ((((log((a_v * theta_v)) + ((a_v - 1.0_dp) * log(x))) + ((theta_v - 1.0_dp) * &
      & log((1.0_dp + (x ** a_v))))) + 1.0_dp) - ((1.0_dp + (x ** a_v)) ** theta_v))
  end function dgenpowerweibull

  pure elemental function pgenpowerweibull(x, a, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp((1.0_dp - ((1.0_dp + (x ** a_v)) ** &
      & theta_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp((1.0_dp - ((1.0_dp + (x ** a_v)) ** theta_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp((1.0_dp - ((1.0_dp + (x ** a_v)) ** theta_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp + (x ** a_v)) ** theta_v))
  end function pgenpowerweibull

  pure elemental function vargenpowerweibull(p, a, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((((1.0_dp - log((1.0_dp - pp))) ** (1.0_dp / theta_v)) - 1.0_dp) ** (1.0_dp / a_v))
  end function vargenpowerweibull

  pure elemental function esgenpowerweibull(p, a, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenpowerweibull(p * s_quad**3, a_v, theta_v)
    end do
  end function esgenpowerweibull

  pure elemental function dchen(x, b, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b, lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: lambda_v
    logical :: log_pdf_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((b_v * lambda_v) * (x ** (b_v - 1.0_dp))) * exp((x ** b_v))) * exp((lambda_v &
      & - (lambda_v * exp((x ** b_v))))))
    if (((log_pdf_v))) res = ((((log((b_v * lambda_v)) + ((b_v - 1.0_dp) * log(x))) + (x ** b_v)) + lambda_v) - &
      & (lambda_v * exp((x ** b_v))))
  end function dchen

  pure elemental function pchen(x, b, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp((lambda_v - (lambda_v * exp((x ** b_v))))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp((lambda_v - (lambda_v * exp((x ** b_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp((lambda_v - (lambda_v * exp((x ** b_v)))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (lambda_v - (lambda_v * exp((x ** b_v))))
  end function pchen

  pure elemental function varchen(p, b, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (log((1.0_dp - ((1.0_dp / lambda_v) * log((1.0_dp - pp))))) ** (1.0_dp / b_v))
  end function varchen

  pure elemental function eschen(p, b, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, lambda
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varchen(p * s_quad**3, b_v, lambda_v)
    end do
  end function eschen

  pure elemental function dxie(x, a, b, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: lambda_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((b_v * lambda_v) * ((x / a_v) ** (b_v - 1.0_dp))) * exp((((a_v * lambda_v) + &
      & ((x / a_v) ** b_v)) - ((a_v * lambda_v) * exp(((x / a_v) ** b_v))))))
    if (((log_pdf_v))) res = ((((log((b_v * lambda_v)) + ((b_v - 1.0_dp) * log((x / a_v)))) + (a_v * lambda_v)) + &
      & ((x / a_v) ** b_v)) - ((a_v * lambda_v) * exp(((x / a_v) ** b_v))))
  end function dxie

  pure elemental function pxie(x, a, b, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp(((a_v * lambda_v) - ((a_v * lambda_v) * &
      & exp(((x / a_v) ** b_v))))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp(((a_v * lambda_v) - ((a_v * lambda_v) * exp(((x &
      & / a_v) ** b_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp(((a_v * lambda_v) - ((a_v * lambda_v) * &
      & exp(((x / a_v) ** b_v)))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((a_v * lambda_v) - ((a_v * lambda_v) * exp(((x / a_v) ** &
      & b_v))))
  end function pxie

  pure elemental function varxie(p, a, b, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (a_v * (log((1.0_dp - (log((1.0_dp - pp)) / (lambda_v * a_v)))) ** (1.0_dp / b_v)))
  end function varxie

  pure elemental function esxie(p, a, b, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, lambda
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varxie(p * s_quad**3, a_v, b_v, lambda_v)
    end do
  end function esxie

  pure elemental function vartl(p, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (((pp ** lambda_v) - ((1.0_dp - pp) ** lambda_v)) / lambda_v)
  end function vartl

  pure elemental function estl(p, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda
    real(dp) :: res
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vartl(p * s_quad**3, lambda_v)
    end do
  end function estl

  pure elemental function varrs(p, b, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (((pp ** b_v) - ((1.0_dp - pp) ** c_v)) / d_v)
  end function varrs

  pure elemental function esrs(p, b, c, d) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, c, d
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    integer :: i
    real(dp) :: s_quad
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varrs(p * s_quad**3, b_v, c_v, d_v)
    end do
  end function esrs

  pure elemental function varfr(p, a, b, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / a_v) * ((((pp ** b_v) - 1.0_dp) / b_v) - ((((1.0_dp - pp) ** c_v) - 1.0_dp) / c_v)))
  end function varfr

  pure elemental function esfr(p, a, b, c) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varfr(p * s_quad**3, a_v, b_v, c_v)
    end do
  end function esfr

  pure elemental function varhl(p, a, b, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((c_v * (pp ** a_v)) / ((1.0_dp - pp) ** b_v))
  end function varhl

  pure elemental function eshl(p, a, b, c) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhl(p * s_quad**3, a_v, b_v, c_v)
    end do
  end function eshl

  pure elemental function dloglog(x, a, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    lambda_v = 2.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * log(lambda_v)) * (x ** (a_v - 1.0_dp))) * (lambda_v ** (x ** a_v))) * &
      & exp((1.0_dp - (lambda_v ** (x ** a_v)))))
    if (((log_pdf_v))) res = (((((log(a_v) + log(log(lambda_v))) + ((a_v - 1.0_dp) * log(x))) + ((x ** a_v) * &
      & log(lambda_v))) + 1.0_dp) - (lambda_v ** (x ** a_v)))
  end function dloglog

  pure elemental function ploglog(x, a, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    lambda_v = 2.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp((1.0_dp - (lambda_v ** (x ** a_v)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp((1.0_dp - (lambda_v ** (x ** a_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp((1.0_dp - (lambda_v ** (x ** a_v))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - (lambda_v ** (x ** a_v)))
  end function ploglog

  pure elemental function varloglog(p, a, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    lambda_v = 2.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((log((1.0_dp - log((1.0_dp - pp)))) / log(lambda_v)) ** (1.0_dp / a_v))
  end function varloglog

  pure elemental function esloglog(p, a, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, lambda
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    lambda_v = 2.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varloglog(p * s_quad**3, a_v, lambda_v)
    end do
  end function esloglog

  pure elemental function dexplog(x, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((-b_v) * (1.0_dp - a_v)) * exp(((-b_v) * x))) / (log(a_v) * (1.0_dp - &
      & ((1.0_dp - a_v) * exp(((-b_v) * x))))))
    if (((log_pdf_v))) res = ((((log(b_v) + log((1.0_dp - a_v))) - (b_v * x)) - log((-log(a_v)))) - log((1.0_dp - &
      & ((1.0_dp - a_v) * exp(((-b_v) * x))))))
  end function dexplog

  pure elemental function pexplog(x, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - (log((1.0_dp - ((1.0_dp - a_v) * exp(((-b_v) * &
      & x))))) / log(a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (log((1.0_dp - ((1.0_dp - a_v) * exp(((-b_v) * &
      & x))))) / log(a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - (log((1.0_dp - ((1.0_dp - a_v) * exp(((-b_v) * &
      & x))))) / log(a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (log((-log((1.0_dp - ((1.0_dp - a_v) * exp(((-b_v) * &
      & x))))))) - log((-log(a_v))))
  end function pexplog

  pure elemental function varexplog(p, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 0.5_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((-(1.0_dp / b_v)) * log(((1.0_dp - (a_v ** (1.0_dp - pp))) / (1.0_dp - a_v))))
  end function varexplog

  pure elemental function esexplog(p, a, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.5_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexplog(p * s_quad**3, a_v, b_v)
    end do
  end function esexplog

  pure elemental function dexpgeo(x, theta, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: theta, lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: lambda_v
    logical :: log_pdf_v
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((lambda_v * theta_v) * exp(((-lambda_v) * x))) * ((1.0_dp - ((1.0_dp - &
      & theta_v) * exp(((-lambda_v) * x)))) ** (-2.0_dp)))
    if (((log_pdf_v))) res = ((log((lambda_v * theta_v)) - (lambda_v * x)) - (2.0_dp * log((1.0_dp - ((1.0_dp - &
      & theta_v) * exp(((-lambda_v) * x)))))))
  end function dexpgeo

  pure elemental function pexpgeo(x, theta, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: theta, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((theta_v * exp(((-lambda_v) * x))) * ((1.0_dp - ((1.0_dp &
      & - theta_v) * exp(((-lambda_v) * x)))) ** (-1.0_dp)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - exp(((-lambda_v) * x))) * ((1.0_dp - &
      & ((1.0_dp - theta_v) * exp(((-lambda_v) * x)))) ** (-1.0_dp)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((log(theta_v) - (lambda_v * x)) - log((1.0_dp - ((1.0_dp - &
      & theta_v) * exp(((-lambda_v) * x))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (log((1.0_dp - exp(((-lambda_v) * x)))) - log((1.0_dp - &
      & ((1.0_dp - theta_v) * exp(((-lambda_v) * x))))))
  end function pexpgeo

  pure elemental function varexpgeo(p, theta, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: theta, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((-(1.0_dp / lambda_v)) * log((pp / (theta_v + ((1.0_dp - theta_v) * pp)))))
  end function varexpgeo

  pure elemental function esexpgeo(p, theta, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: theta, lambda
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexpgeo(p * s_quad**3, theta_v, lambda_v)
    end do
  end function esexpgeo

end module vares_distributions_09
