! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_10
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dexppois, pexppois, varexppois, esexppois, dtl2, ptl2
  public :: vartl2, estl2, dquad, pquad, varquad, esquad
  public :: dschabe, pschabe, varschabe, esschabe, dbs, pbs
  public :: varbs, esbs, dgev, pgev, vargev, esgev
contains
  pure elemental function dexppois(x, b, lambda, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((b_v * lambda_v) * exp(((((-b_v) * x) - lambda_v) + (lambda_v * exp(((-b_v) * &
      & x)))))) / (1.0_dp - exp((-lambda_v))))
    if (((log_pdf_v))) res = ((((log((b_v * lambda_v)) - (b_v * x)) - lambda_v) + (lambda_v * exp(((-b_v) * x)))) - &
      & log((1.0_dp - exp((-lambda_v)))))
  end function dexppois

  pure elemental function pexppois(x, b, lambda, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp - exp(((-lambda_v) + (lambda_v * exp(((-b_v) * &
      & x)))))) / (1.0_dp - exp((-lambda_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((exp(((-lambda_v) + (lambda_v * exp(((-b_v) * &
      & x))))) - exp((-lambda_v))) / (1.0_dp - exp((-lambda_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (log((1.0_dp - exp(((-lambda_v) + (lambda_v * exp(((-b_v) * &
      & x))))))) - log((1.0_dp - exp((-lambda_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (log((exp(((-lambda_v) + (lambda_v * exp(((-b_v) * x))))) &
      & - exp((-lambda_v)))) - log((1.0_dp - exp((-lambda_v)))))
  end function pexppois

  pure elemental function varexppois(p, b, lambda, log_p, lower_tail) result(res)
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
    res = ((-(1.0_dp / b_v)) * log((((1.0_dp / lambda_v) * log(((1.0_dp - pp) + (pp * exp((-lambda_v)))))) + 1.0_dp)))
  end function varexppois

  pure elemental function esexppois(p, b, lambda) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexppois(p * s_quad**3, b_v, lambda_v)
    end do
  end function esexppois

  pure elemental function dtl2(x, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: b_v
    logical :: log_pdf_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((2.0_dp * b_v) * ((x * (2.0_dp - x)) ** (b_v - 1.0_dp))) * (1.0_dp - x))
    if (((log_pdf_v))) res = ((log((2.0_dp * b_v)) + ((b_v - 1.0_dp) * log((x * (2.0_dp - x))))) + log((1.0_dp - x)))
  end function dtl2

  pure elemental function ptl2(x, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((x * (2.0_dp - x)) ** b_v)
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((x * (2.0_dp - x)) ** b_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (b_v * log((x * (2.0_dp - x))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((x * (2.0_dp - x)) ** b_v)))
  end function ptl2

  pure elemental function vartl2(p, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (1.0_dp - sqrt((1.0_dp - (pp ** (1.0_dp / b_v)))))
  end function vartl2

  pure elemental function estl2(p, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b
    real(dp) :: res
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
    b_v = 1.0_dp
    if (present(b)) b_v = b
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vartl2(p * s_quad**3, b_v)
    end do
  end function estl2

  pure elemental function dquad(x, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    real(dp) :: alpha_local
    real(dp) :: beta_local
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    alpha_local = (12.0_dp / ((b_v - a_v) ** 3.0_dp))
    beta_local = ((a_v + b_v) / 2.0_dp)
    res = x
    if (((.not. log_pdf_v))) res = (alpha_local * ((x - beta_local) ** 2.0_dp))
    if (((log_pdf_v))) res = (log(alpha_local) + (2.0_dp * log((x - beta_local))))
  end function dquad

  pure elemental function pquad(x, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: alpha_local
    real(dp) :: beta_local
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    alpha_local = (12.0_dp / ((b_v - a_v) ** 3.0_dp))
    beta_local = ((a_v + b_v) / 2.0_dp)
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((alpha_local / 3.0_dp) * (((x - beta_local) ** 3.0_dp) + &
      & ((beta_local - a_v) ** 3.0_dp)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((alpha_local / 3.0_dp) * (((x - &
      & beta_local) ** 3.0_dp) + ((beta_local - a_v) ** 3.0_dp))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (log((alpha_local / 3.0_dp)) + log((((x - beta_local) ** &
      & 3.0_dp) + ((beta_local - a_v) ** 3.0_dp))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((alpha_local / 3.0_dp) * (((x - &
      & beta_local) ** 3.0_dp) + ((beta_local - a_v) ** 3.0_dp)))))
  end function pquad

  pure elemental function varquad(p, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    real(dp) :: alpha_local
    real(dp) :: beta_local
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    alpha_local = (12.0_dp / ((b_v - a_v) ** 3.0_dp))
    beta_local = ((a_v + b_v) / 2.0_dp)
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (beta_local + (r_sign((((3.0_dp * pp) / alpha_local) - ((beta_local - a_v) ** 3.0_dp))) * (abs((((3.0_dp &
      & * pp) / alpha_local) - ((beta_local - a_v) ** 3.0_dp))) ** (1.0_dp / 3.0_dp))))
  end function varquad

  pure elemental function esquad(p, a, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.0_dp
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varquad(p * s_quad**3, a_v, b_v)
    end do
  end function esquad

  pure elemental function dschabe(x, gamma, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: gamma, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    gamma_v = 0.5_dp
    if (present(gamma)) gamma_v = gamma
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((1.0_dp / theta_v) * ((2.0_dp * gamma_v) + (((1.0_dp - gamma_v) * x) / &
      & theta_v))) * ((gamma_v + (x / theta_v)) ** (-2.0_dp)))
    if (((log_pdf_v))) res = (((-log(theta_v)) + log(((2.0_dp * gamma_v) + (((1.0_dp - gamma_v) * x) / theta_v)))) &
      & - (2.0_dp * log((gamma_v + (x / theta_v)))))
  end function dschabe

  pure elemental function pschabe(x, gamma, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: gamma, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    gamma_v = 0.5_dp
    if (present(gamma)) gamma_v = gamma
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (((1.0_dp + gamma_v) * x) / (x + (gamma_v * theta_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((gamma_v * (theta_v - x)) / (x + (gamma_v * &
      & theta_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((log((1.0_dp + gamma_v)) + log(x)) - log((x + (gamma_v * &
      & theta_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((log(gamma_v) + log((theta_v - x))) - log((x + (gamma_v &
      & * theta_v))))
  end function pschabe

  pure elemental function varschabe(p, gamma, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: gamma, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    gamma_v = 0.5_dp
    if (present(gamma)) gamma_v = gamma
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (((pp * gamma_v) * theta_v) / ((1.0_dp + gamma_v) - pp))
  end function varschabe

  pure elemental function esschabe(p, gamma, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: gamma, theta
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    gamma_v = 0.5_dp
    if (present(gamma)) gamma_v = gamma
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varschabe(p * s_quad**3, gamma_v, theta_v)
    end do
  end function esschabe

  pure elemental function dbs(x, gamma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: gamma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: gamma_v
    logical :: log_pdf_v
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((x ** (1.0_dp / 2.0_dp)) + (x ** ((-1.0_dp) / 2.0_dp))) * normal_pdf((((x ** &
      & (1.0_dp / 2.0_dp)) - (x ** ((-1.0_dp) / 2.0_dp))) / gamma_v))) / ((2.0_dp * gamma_v) * x))
    if (((log_pdf_v))) res = (((log(((x ** (1.0_dp / 2.0_dp)) + (x ** ((-1.0_dp) / 2.0_dp)))) + normal_pdf((((x ** &
      & (1.0_dp / 2.0_dp)) - (x ** ((-1.0_dp) / 2.0_dp))) / gamma_v), log_pdf=.true.)) - log((2.0_dp * gamma_v))) - &
      & log(x))
  end function dbs

  pure elemental function pbs(x, gamma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: gamma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: gamma_v
    logical :: log_p_v
    logical :: lower_tail_v
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = normal_cdf((((x ** (1.0_dp / 2.0_dp)) - (x ** ((-1.0_dp) / 2.0_dp))) / gamma_v), log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbs

  pure elemental function varbs(p, gamma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: gamma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: gamma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (0.25_dp * (((gamma_v * normal_quantile(pp)) + sqrt((4.0_dp + ((gamma_v ** 2.0_dp) * (normal_quantile(pp) &
      & ** 2.0_dp))))) ** 2.0_dp))
  end function varbs

  pure elemental function esbs(p, gamma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: gamma
    real(dp) :: res
    real(dp) :: gamma_v
    integer :: i
    real(dp) :: s_quad
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbs(p * s_quad**3, gamma_v)
    end do
  end function esbs

  pure elemental function dgev(x, mu, sigma, xi, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma, xi
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: xi_v
    logical :: log_pdf_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    xi_v = 1.0_dp
    if (present(xi)) xi_v = xi
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((1.0_dp / sigma_v) * ((1.0_dp + ((xi_v * (x - mu_v)) / sigma_v)) ** &
      & (((-1.0_dp) / xi_v) - 1.0_dp))) * exp((-((1.0_dp + ((xi_v * (x - mu_v)) / sigma_v)) ** ((-1.0_dp) / xi_v)))))
    if (((log_pdf_v))) res = (((-log(sigma_v)) - (((1.0_dp / xi_v) + 1.0_dp) * log((1.0_dp + ((xi_v * (x - mu_v)) / &
      & sigma_v))))) - ((1.0_dp + ((xi_v * (x - mu_v)) / sigma_v)) ** ((-1.0_dp) / xi_v)))
  end function dgev

  pure elemental function pgev(x, mu, sigma, xi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma, xi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: xi_v
    logical :: log_p_v
    logical :: lower_tail_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    xi_v = 1.0_dp
    if (present(xi)) xi_v = xi
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = exp((-((1.0_dp + ((xi_v * (x - mu_v)) / sigma_v)) ** &
      & ((-1.0_dp) / xi_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - exp((-((1.0_dp + ((xi_v * (x - mu_v)) / &
      & sigma_v)) ** ((-1.0_dp) / xi_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (-((1.0_dp + ((xi_v * (x - mu_v)) / sigma_v)) ** ((-1.0_dp) / &
      & xi_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - exp((-((1.0_dp + ((xi_v * (x - mu_v)) / &
      & sigma_v)) ** ((-1.0_dp) / xi_v))))))
  end function pgev

  pure elemental function vargev(p, mu, sigma, xi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma, xi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: xi_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    xi_v = 1.0_dp
    if (present(xi)) xi_v = xi
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((mu_v - (sigma_v / xi_v)) + ((sigma_v / xi_v) * ((-log(pp)) ** (-xi_v))))
  end function vargev

  pure elemental function esgev(p, mu, sigma, xi) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma, xi
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: xi_v
    integer :: i
    real(dp) :: s_quad
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    xi_v = 1.0_dp
    if (present(xi)) xi_v = xi
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargev(p * s_quad**3, mu_v, sigma_v, xi_v)
    end do
  end function esgev

end module vares_distributions_10
