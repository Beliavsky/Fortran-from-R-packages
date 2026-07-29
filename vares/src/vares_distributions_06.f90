! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_06
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dmlaplace, pmlaplace, varmlaplace, esmlaplace, dloglaplace, ploglaplace
  public :: varloglaplace, esloglaplace, dasylaplace, pasylaplace, varasylaplace, esasylaplace
  public :: dasypower, pasypower, varasypower, esasypower, dlogistic, plogistic
  public :: varlogistic, eslogistic, dsecant, psecant, varsecant, essecant
  public :: dgenlogis, pgenlogis, vargenlogis, esgenlogis, dgenlogis3, pgenlogis3
  public :: vargenlogis3, esgenlogis3, dgenlogis4, pgenlogis4, vargenlogis4, esgenlogis4
  public :: dhalflogis, phalflogis, varhalflogis, eshalflogis, dloglogis, ploglogis
  public :: varloglogis, esloglogis, dkumloglogis, pkumloglogis, varkumloglogis, eskumloglogis
contains
  pure elemental function dmlaplace(x, theta, phi, psi, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: theta, phi, psi
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: phi_v
    real(dp) :: psi_v
    logical :: log_pdf_v
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    psi_v = 1.0_dp
    if (present(psi)) psi_v = psi
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= theta_v)))) res = ((1.0_dp / (2.0_dp * psi_v)) * exp(((x - theta_v) / &
      & psi_v)))
    if ((((.not. log_pdf_v)) .and. ((x > theta_v)))) res = ((1.0_dp / (2.0_dp * phi_v)) * exp(((theta_v - x) / &
      & phi_v)))
    if ((((log_pdf_v)) .and. ((x <= theta_v)))) res = (((x - theta_v) / psi_v) - log((2.0_dp * psi_v)))
    if ((((log_pdf_v)) .and. ((x > theta_v)))) res = (((theta_v - x) / phi_v) - log((2.0_dp * phi_v)))
  end function dmlaplace

  pure elemental function pmlaplace(x, theta, phi, psi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: theta, phi, psi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: phi_v
    real(dp) :: psi_v
    logical :: log_p_v
    logical :: lower_tail_v
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    psi_v = 1.0_dp
    if (present(psi)) psi_v = psi
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (0.5_dp * exp(((x - theta_v) / &
      & psi_v)))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = (1.0_dp - (0.5_dp * exp(((theta_v - &
      & x) / phi_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = (1.0_dp - (0.5_dp * exp(((x &
      & - theta_v) / psi_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = (0.5_dp * exp(((theta_v - x) &
      & / phi_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (((x - theta_v) / psi_v) - log(2.0_dp))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = log((1.0_dp - (0.5_dp * exp(((theta_v - &
      & x) / phi_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = log((1.0_dp - (0.5_dp * exp(((x - &
      & theta_v) / psi_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = (((theta_v - x) / phi_v) - &
      & log(2.0_dp))
  end function pmlaplace

  pure elemental function varmlaplace(p, theta, phi, psi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: theta, phi, psi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: phi_v
    real(dp) :: psi_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    psi_v = 1.0_dp
    if (present(psi)) psi_v = psi
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= 0.5_dp))) res = (theta_v + (psi_v * log((2.0_dp * pp))))
    if (((pp > 0.5_dp))) res = (theta_v - (phi_v * log((2.0_dp * (1.0_dp - pp)))))
  end function varmlaplace

  pure elemental function esmlaplace(p, theta, phi, psi) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: theta, phi, psi
    real(dp) :: res
    real(dp) :: theta_v
    real(dp) :: phi_v
    real(dp) :: psi_v
    integer :: i
    real(dp) :: s_quad
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    psi_v = 1.0_dp
    if (present(psi)) psi_v = psi
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varmlaplace(p * s_quad**3, theta_v, phi_v, psi_v)
    end do
  end function esmlaplace

  pure elemental function dloglaplace(x, a, b, delta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, delta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: delta_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    delta_v = 0.0_dp
    if (present(delta)) delta_v = delta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= delta_v)))) res = (((a_v * b_v) * (x ** (b_v - 1.0_dp))) / ((delta_v ** &
      & b_v) * (a_v + b_v)))
    if ((((.not. log_pdf_v)) .and. ((x > delta_v)))) res = (((a_v * b_v) * (delta_v ** a_v)) / ((x ** (a_v + &
      & 1.0_dp)) * (a_v + b_v)))
    if ((((log_pdf_v)) .and. ((x <= delta_v)))) res = (((log((a_v * b_v)) + ((b_v - 1.0_dp) * log(x))) - (b_v * &
      & log(delta_v))) - log((a_v + b_v)))
    if ((((log_pdf_v)) .and. ((x > delta_v)))) res = (((log((a_v * b_v)) - ((a_v + 1.0_dp) * log(x))) + (a_v * &
      & log(delta_v))) - log((a_v + b_v)))
  end function dloglaplace

  pure elemental function ploglaplace(x, a, b, delta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, delta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: delta_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    delta_v = 0.0_dp
    if (present(delta)) delta_v = delta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= delta_v)))) res = ((a_v / (a_v + b_v)) * ((x / &
      & delta_v) ** b_v))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > delta_v)))) res = (1.0_dp - ((b_v / (a_v + b_v)) * &
      & ((delta_v / x) ** a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= delta_v)))) res = (1.0_dp - ((a_v / (a_v + &
      & b_v)) * ((x / delta_v) ** b_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > delta_v)))) res = ((b_v / (a_v + b_v)) * &
      & ((delta_v / x) ** a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= delta_v)))) res = (log((a_v / (a_v + b_v))) + (b_v * &
      & log((x / delta_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > delta_v)))) res = log((1.0_dp - ((b_v / (a_v + b_v)) * &
      & ((delta_v / x) ** a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= delta_v)))) res = log((1.0_dp - ((a_v / (a_v + b_v)) &
      & * ((x / delta_v) ** b_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > delta_v)))) res = (log((b_v / (a_v + b_v))) + (a_v * &
      & log((delta_v / x))))
  end function ploglaplace

  pure elemental function varloglaplace(p, a, b, delta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, delta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: delta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    delta_v = 0.0_dp
    if (present(delta)) delta_v = delta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= (a_v / (a_v + b_v))))) res = (delta_v * (((pp * (a_v + b_v)) / a_v) ** (1.0_dp / b_v)))
    if (((pp > (a_v / (a_v + b_v))))) res = (delta_v * ((((1.0_dp - pp) * (a_v + b_v)) / b_v) ** ((-1.0_dp) / a_v)))
  end function varloglaplace

  pure elemental function esloglaplace(p, a, b, delta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, delta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: delta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    delta_v = 0.0_dp
    if (present(delta)) delta_v = delta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varloglaplace(p * s_quad**3, a_v, b_v, delta_v)
    end do
  end function esloglaplace

  pure elemental function dasylaplace(x, tau, kappa, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: tau, kappa, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: tau_v
    real(dp) :: kappa_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    tau_v = 1.0_dp
    if (present(tau)) tau_v = tau
    kappa_v = 1.0_dp
    if (present(kappa)) kappa_v = kappa
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x >= theta_v)))) res = (((sqrt(2.0_dp) * kappa_v) * exp(((((-sqrt(2.0_dp)) * &
      & kappa_v) * abs((x - theta_v))) / tau_v))) / (tau_v * (1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((.not. log_pdf_v)) .and. ((x < theta_v)))) res = (((sqrt(2.0_dp) * kappa_v) * exp((((-sqrt(2.0_dp)) * &
      & abs((x - theta_v))) / (kappa_v * tau_v)))) / (tau_v * (1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((log_pdf_v)) .and. ((x >= theta_v)))) res = (((((0.5_dp * log(2.0_dp)) + log(kappa_v)) - (((sqrt(2.0_dp) &
      & * kappa_v) * abs((x - theta_v))) / tau_v)) - log(tau_v)) - log((1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((log_pdf_v)) .and. ((x < theta_v)))) res = (((((0.5_dp * log(2.0_dp)) + log(kappa_v)) - ((sqrt(2.0_dp) * &
      & abs((x - theta_v))) / (kappa_v * tau_v))) - log(tau_v)) - log((1.0_dp + (kappa_v ** 2.0_dp))))
  end function dasylaplace

  pure elemental function pasylaplace(x, tau, kappa, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: tau, kappa, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: tau_v
    real(dp) :: kappa_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    tau_v = 1.0_dp
    if (present(tau)) tau_v = tau
    kappa_v = 1.0_dp
    if (present(kappa)) kappa_v = kappa
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x >= theta_v)))) res = (1.0_dp - (exp((((sqrt(2.0_dp) * &
      & kappa_v) * (theta_v - x)) / tau_v)) / (1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x < theta_v)))) res = (((kappa_v ** 2.0_dp) * &
      & exp(((sqrt(2.0_dp) * (x - theta_v)) / (kappa_v * tau_v)))) / (1.0_dp + (kappa_v ** 2.0_dp)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x >= theta_v)))) res = (exp((((sqrt(2.0_dp) * &
      & kappa_v) * (theta_v - x)) / tau_v)) / (1.0_dp + (kappa_v ** 2.0_dp)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x < theta_v)))) res = (1.0_dp - (((kappa_v ** &
      & 2.0_dp) * exp(((sqrt(2.0_dp) * (x - theta_v)) / (kappa_v * tau_v)))) / (1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x >= theta_v)))) res = log((1.0_dp - (exp((((sqrt(2.0_dp) * &
      & kappa_v) * (theta_v - x)) / tau_v)) / (1.0_dp + (kappa_v ** 2.0_dp)))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x < theta_v)))) res = (((2.0_dp * log(kappa_v)) + &
      & ((sqrt(2.0_dp) * (x - theta_v)) / (kappa_v * tau_v))) - log((1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x >= theta_v)))) res = ((((sqrt(2.0_dp) * kappa_v) * &
      & (theta_v - x)) / tau_v) - log((1.0_dp + (kappa_v ** 2.0_dp))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x < theta_v)))) res = log((1.0_dp - (((kappa_v ** 2.0_dp) &
      & * exp(((sqrt(2.0_dp) * (x - theta_v)) / (kappa_v * tau_v)))) / (1.0_dp + (kappa_v ** 2.0_dp)))))
  end function pasylaplace

  pure elemental function varasylaplace(p, tau, kappa, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: tau, kappa, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: tau_v
    real(dp) :: kappa_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    tau_v = 1.0_dp
    if (present(tau)) tau_v = tau
    kappa_v = 1.0_dp
    if (present(kappa)) kappa_v = kappa
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp >= ((kappa_v ** 2.0_dp) / (1.0_dp + (kappa_v ** 2.0_dp)))))) res = (theta_v - ((tau_v * log(((1.0_dp - &
      & pp) * (1.0_dp + (kappa_v ** 2.0_dp))))) / (sqrt(2.0_dp) * kappa_v)))
    if (((pp < ((kappa_v ** 2.0_dp) / (1.0_dp + (kappa_v ** 2.0_dp)))))) res = (theta_v + (((tau_v * kappa_v) * &
      & log((pp * (1.0_dp + (kappa_v ** (-2.0_dp)))))) / sqrt(2.0_dp)))
  end function varasylaplace

  pure elemental function esasylaplace(p, tau, kappa, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: tau, kappa, theta
    real(dp) :: res
    real(dp) :: tau_v
    real(dp) :: kappa_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    tau_v = 1.0_dp
    if (present(tau)) tau_v = tau
    kappa_v = 1.0_dp
    if (present(kappa)) kappa_v = kappa
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varasylaplace(p * s_quad**3, tau_v, kappa_v, theta_v)
    end do
  end function esasylaplace

  pure elemental function dasypower(x, a, lambda, delta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, lambda, delta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    real(dp) :: delta_v
    logical :: log_pdf_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    delta_v = 1.0_dp
    if (present(delta)) delta_v = delta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= 0.0_dp)))) res = (((delta_v ** (1.0_dp / lambda_v)) * exp(((-delta_v) * &
      & (((-x) / a_v) ** lambda_v)))) / gamma_fn((1.0_dp + (1.0_dp / lambda_v))))
    if ((((.not. log_pdf_v)) .and. ((x > 0.0_dp)))) res = (((delta_v ** (1.0_dp / lambda_v)) * exp(((-delta_v) * &
      & ((x / (1.0_dp - a_v)) ** lambda_v)))) / gamma_fn((1.0_dp + (1.0_dp / lambda_v))))
    if ((((log_pdf_v)) .and. ((x <= 0.0_dp)))) res = ((((1.0_dp / lambda_v) * log(delta_v)) - (delta_v * (((-x) / &
      & a_v) ** lambda_v))) - log_gamma_fn((1.0_dp + (1.0_dp / lambda_v))))
    if ((((log_pdf_v)) .and. ((x > 0.0_dp)))) res = ((((1.0_dp / lambda_v) * log(delta_v)) - (delta_v * ((x / &
      & (1.0_dp - a_v)) ** lambda_v))) - log_gamma_fn((1.0_dp + (1.0_dp / lambda_v))))
  end function dasypower

  pure elemental function pasypower(x, a, lambda, delta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, lambda, delta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    real(dp) :: delta_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    delta_v = 1.0_dp
    if (present(delta)) delta_v = delta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= 0.0_dp)))) res = (a_v - (a_v * gamma_cdf((delta_v * &
      & (((-x) / a_v) ** lambda_v)), shape=(1.0_dp / lambda_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > 0.0_dp)))) res = (a_v + ((1.0_dp - a_v) * &
      & gamma_cdf((delta_v * ((x / (1.0_dp - a_v)) ** lambda_v)), shape=(1.0_dp / lambda_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= 0.0_dp)))) res = ((1.0_dp - a_v) + (a_v * &
      & gamma_cdf((delta_v * (((-x) / a_v) ** lambda_v)), shape=(1.0_dp / lambda_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > 0.0_dp)))) res = ((1.0_dp - a_v) * (1.0_dp - &
      & gamma_cdf((delta_v * ((x / (1.0_dp - a_v)) ** lambda_v)), shape=(1.0_dp / lambda_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= 0.0_dp)))) res = log((a_v - (a_v * gamma_cdf((delta_v * &
      & (((-x) / a_v) ** lambda_v)), shape=(1.0_dp / lambda_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > 0.0_dp)))) res = log((a_v + ((1.0_dp - a_v) * &
      & gamma_cdf((delta_v * ((x / (1.0_dp - a_v)) ** lambda_v)), shape=(1.0_dp / lambda_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= 0.0_dp)))) res = log(((1.0_dp - a_v) + (a_v * &
      & gamma_cdf((delta_v * (((-x) / a_v) ** lambda_v)), shape=(1.0_dp / lambda_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > 0.0_dp)))) res = (log((1.0_dp - a_v)) + log((1.0_dp - &
      & gamma_cdf((delta_v * ((x / (1.0_dp - a_v)) ** lambda_v)), shape=(1.0_dp / lambda_v)))))
  end function pasypower

  pure elemental function varasypower(p, a, lambda, delta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, lambda, delta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    real(dp) :: delta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 0.5_dp
    if (present(a)) a_v = a
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    delta_v = 1.0_dp
    if (present(delta)) delta_v = delta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= a_v))) res = (((-a_v) * (gamma_quantile((1.0_dp - (pp / a_v)), shape=(1.0_dp / lambda_v)) ** &
      & (1.0_dp / lambda_v))) / (delta_v ** (1.0_dp / lambda_v)))
    if (((pp > a_v))) res = (((1.0_dp - a_v) * (gamma_quantile((1.0_dp - ((1.0_dp - pp) / (1.0_dp - a_v))), &
      & shape=(1.0_dp / lambda_v)) ** (1.0_dp / lambda_v))) / (delta_v ** (1.0_dp / lambda_v)))
  end function varasypower

  pure elemental function esasypower(p, a, lambda, delta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, lambda, delta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    real(dp) :: delta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.5_dp
    if (present(a)) a_v = a
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    delta_v = 1.0_dp
    if (present(delta)) delta_v = delta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varasypower(p * s_quad**3, a_v, lambda_v, delta_v)
    end do
  end function esasypower

  pure elemental function dlogistic(x, mu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = logistic_pdf(x, location=mu_v, scale=sigma_v, log_pdf=log_pdf_v)
  end function dlogistic

  pure elemental function plogistic(x, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = logistic_cdf(x, location=mu_v, scale=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function plogistic

  pure elemental function varlogistic(p, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = logistic_quantile(pp, location=mu_v, scale=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varlogistic

  pure elemental function eslogistic(p, mu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlogistic(p * s_quad**3, mu_v, sigma_v)
    end do
  end function eslogistic

  pure elemental function dsecant(x, log_pdf) result(res)
    real(dp), intent(in) :: x
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    logical :: log_pdf_v
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (0.5_dp / cosh(((pi * x) / 2.0_dp)))
    if (((log_pdf_v))) res = ((-log(2.0_dp)) - log(cosh(((pi * x) / 2.0_dp))))
  end function dsecant

  pure elemental function psecant(x, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    logical :: log_p_v
    logical :: lower_tail_v
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((2.0_dp / pi) * atan(exp(((pi * x) / 2.0_dp))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((2.0_dp / pi) * atan(exp(((pi * x) / &
      & 2.0_dp)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((log(2.0_dp) - log(pi)) + log(atan(exp(((pi * x) / 2.0_dp)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((2.0_dp / pi) * atan(exp(((pi * x) / &
      & 2.0_dp))))))
  end function psecant

  pure elemental function varsecant(p, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((2.0_dp / pi) * log(tan(((pi * pp) / 2.0_dp))))
  end function varsecant

  pure elemental function essecant(p) result(res)
    real(dp), intent(in) :: p
    real(dp) :: res
    integer :: i
    real(dp) :: s_quad
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varsecant(p * s_quad**3)
    end do
  end function essecant

  pure elemental function dgenlogis(x, a, mu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, mu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((a_v / sigma_v) * exp(((-(x - mu_v)) / sigma_v))) * ((1.0_dp + exp(((-(x - &
      & mu_v)) / sigma_v))) ** ((-1.0_dp) - a_v)))
    if (((log_pdf_v))) res = (((log(a_v) - log(sigma_v)) - ((x - mu_v) / sigma_v)) - ((1.0_dp + a_v) * log((1.0_dp &
      & + exp(((-(x - mu_v)) / sigma_v))))))
  end function dgenlogis

  pure elemental function pgenlogis(x, a, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp + exp(((-(x - mu_v)) / sigma_v))) ** (-a_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp + exp(((-(x - mu_v)) / &
      & sigma_v))) ** (-a_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((-a_v) * log((1.0_dp + exp(((-(x - mu_v)) / sigma_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((1.0_dp + exp(((-(x - mu_v)) / sigma_v))) &
      & ** (-a_v))))
  end function pgenlogis

  pure elemental function vargenlogis(p, a, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
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
    res = (mu_v - (sigma_v * log(((pp ** ((-1.0_dp) / a_v)) - 1.0_dp))))
  end function vargenlogis

  pure elemental function esgenlogis(p, a, mu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, mu, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenlogis(p * s_quad**3, a_v, mu_v, sigma_v)
    end do
  end function esgenlogis

  pure elemental function dgenlogis3(x, alpha, mu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha, mu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((exp(((alpha_v * (x - mu_v)) / sigma_v)) * ((1.0_dp + exp(((x - mu_v) / &
      & sigma_v))) ** ((-1.0_dp) - alpha_v))) / (sigma_v * beta_fn(alpha_v, alpha_v)))
    if (((log_pdf_v))) res = (((((alpha_v * (x - mu_v)) / sigma_v) - ((1.0_dp + alpha_v) * log((1.0_dp + exp(((x - &
      & mu_v) / sigma_v)))))) - log(sigma_v)) - log_beta_fn(alpha_v, alpha_v))
  end function dgenlogis3

  pure elemental function pgenlogis3(x, alpha, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf((1.0_dp / (1.0_dp + exp(((mu_v - x) / sigma_v)))), shape1=alpha_v, shape2=alpha_v, &
      & log_p=log_p_v, lower_tail=lower_tail_v)
  end function pgenlogis3

  pure elemental function vargenlogis3(p, alpha, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
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
    res = (mu_v - (sigma_v * log(((1.0_dp / beta_quantile(pp, shape1=alpha_v, shape2=alpha_v)) - 1.0_dp))))
  end function vargenlogis3

  pure elemental function esgenlogis3(p, alpha, mu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha, mu, sigma
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenlogis3(p * s_quad**3, alpha_v, mu_v, sigma_v)
    end do
  end function esgenlogis3

  pure elemental function dgenlogis4(x, a, alpha, mu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, alpha, mu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((exp((((-alpha_v) * (x - mu_v)) / sigma_v)) * ((1.0_dp + exp(((-(x - mu_v)) / &
      & sigma_v))) ** ((-a_v) - alpha_v))) / (sigma_v * beta_fn(a_v, alpha_v)))
    if (((log_pdf_v))) res = ((((((-alpha_v) * (x - mu_v)) / sigma_v) - ((a_v + alpha_v) * log((1.0_dp + exp(((-(x &
      & - mu_v)) / sigma_v)))))) - log(sigma_v)) - log_beta_fn(a_v, alpha_v))
  end function dgenlogis4

  pure elemental function pgenlogis4(x, a, alpha, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, alpha, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf((1.0_dp / (1.0_dp + exp(((mu_v - x) / sigma_v)))), shape1=alpha_v, shape2=a_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pgenlogis4

  pure elemental function vargenlogis4(p, a, alpha, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, alpha, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
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
    res = (mu_v - (sigma_v * log(((1.0_dp / beta_quantile(pp, shape1=alpha_v, shape2=a_v)) - 1.0_dp))))
  end function vargenlogis4

  pure elemental function esgenlogis4(p, a, alpha, mu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, alpha, mu, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: alpha_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenlogis4(p * s_quad**3, a_v, alpha_v, mu_v, sigma_v)
    end do
  end function esgenlogis4

  pure elemental function dhalflogis(x, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((2.0_dp * lambda_v) * exp(((-lambda_v) * x))) * ((1.0_dp + exp(((-lambda_v) * &
      & x))) ** (-2.0_dp)))
    if (((log_pdf_v))) res = ((log((2.0_dp * lambda_v)) - (lambda_v * x)) - (2.0_dp * log((1.0_dp + &
      & exp(((-lambda_v) * x))))))
  end function dhalflogis

  pure elemental function phalflogis(x, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp - exp(((-lambda_v) * x))) / (1.0_dp + &
      & exp(((-lambda_v) * x))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) / &
      & (1.0_dp + exp(((-lambda_v) * x)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log(((1.0_dp - exp(((-lambda_v) * x))) / (1.0_dp + &
      & exp(((-lambda_v) * x)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) / &
      & (1.0_dp + exp(((-lambda_v) * x))))))
  end function phalflogis

  pure elemental function varhalflogis(p, lambda, log_p, lower_tail) result(res)
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
    res = ((-(1.0_dp / lambda_v)) * log(((1.0_dp - pp) / (1.0_dp + pp))))
  end function varhalflogis

  pure elemental function eshalflogis(p, lambda) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhalflogis(p * s_quad**3, lambda_v)
    end do
  end function eshalflogis

  pure elemental function dloglogis(x, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((b_v * (a_v ** b_v)) * (x ** (b_v - 1.0_dp))) * (((a_v ** b_v) + (x ** b_v)) &
      & ** (-2.0_dp)))
    if (((log_pdf_v))) res = (((log(b_v) + (b_v * log(a_v))) + ((b_v - 1.0_dp) * log(x))) - (2.0_dp * log(((a_v ** &
      & b_v) + (x ** b_v)))))
  end function dloglogis

  pure elemental function ploglogis(x, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((x ** b_v) / ((a_v ** b_v) + (x ** b_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((a_v ** b_v) / ((a_v ** b_v) + (x ** b_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((b_v * log(x)) - log(((a_v ** b_v) + (x ** b_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((b_v * log(a_v)) - log(((a_v ** b_v) + (x ** b_v))))
  end function ploglogis

  pure elemental function varloglogis(p, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
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
    res = (a_v * ((pp / (1.0_dp - pp)) ** (1.0_dp / b_v)))
  end function varloglogis

  pure elemental function esloglogis(p, a, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varloglogis(p * s_quad**3, a_v, b_v)
    end do
  end function esloglogis

  pure elemental function dkumloglogis(x, a, b, alpha, beta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, beta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: beta_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    beta_v = 1.0_dp
    if (present(beta)) beta_v = beta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((((a_v * b_v) * beta_v) * (alpha_v ** beta_v)) * (x ** ((a_v * beta_v) - &
      & 1.0_dp))) * (((alpha_v ** beta_v) + (x ** beta_v)) ** ((-a_v) - 1.0_dp))) * ((1.0_dp - ((x ** (a_v * &
      & beta_v)) * (((alpha_v ** beta_v) + (x ** beta_v)) ** (-a_v)))) ** (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((((log(((a_v * b_v) * beta_v)) + (beta_v * log(alpha_v))) + (((a_v * beta_v) - &
      & 1.0_dp) * log(x))) - ((a_v + 1.0_dp) * log(((alpha_v ** beta_v) + (x ** beta_v))))) + ((b_v - 1.0_dp) * &
      & log((1.0_dp - ((x ** (a_v * beta_v)) * (((alpha_v ** beta_v) + (x ** beta_v)) ** (-a_v)))))))
  end function dkumloglogis

  pure elemental function pkumloglogis(x, a, b, alpha, beta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, beta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: beta_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    beta_v = 1.0_dp
    if (present(beta)) beta_v = beta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - ((x ** (a_v * beta_v)) * (((alpha_v &
      & ** beta_v) + (x ** beta_v)) ** (-a_v)))) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - ((x ** (a_v * beta_v)) * (((alpha_v ** &
      & beta_v) + (x ** beta_v)) ** (-a_v)))) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - ((x ** (a_v * beta_v)) * (((alpha_v ** &
      & beta_v) + (x ** beta_v)) ** (-a_v)))) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - ((x ** (a_v * beta_v)) * (((alpha_v &
      & ** beta_v) + (x ** beta_v)) ** (-a_v))))))
  end function pkumloglogis

  pure elemental function varkumloglogis(p, a, b, alpha, beta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, beta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: beta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    beta_v = 1.0_dp
    if (present(beta)) beta_v = beta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    res = (alpha_v * ((((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v)) / (1.0_dp - ((1.0_dp - &
      & ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v)))) ** (1.0_dp / beta_v)))
  end function varkumloglogis

  pure elemental function eskumloglogis(p, a, b, alpha, beta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, beta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: beta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    beta_v = 1.0_dp
    if (present(beta)) beta_v = beta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumloglogis(p * s_quad**3, a_v, b_v, alpha_v, beta_v)
    end do
  end function eskumloglogis

end module vares_distributions_06
