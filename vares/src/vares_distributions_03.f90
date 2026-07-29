! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_03
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dloggamma, ploggamma, varloggamma, esloggamma, dinvgamma, pinvgamma
  public :: varinvgamma, esinvgamma, dstacygamma, pstacygamma, varstacygamma, esstacygamma
  public :: dbetadist, pbetadist, varbetadist, esbetadist, duniform, puniform
  public :: varuniform, esuniform, dgenunif, pgenunif, vargenunif, esgenunif
  public :: dpower1, ppower1, varpower1, espower1, dpower2, ppower2
  public :: varpower2, espower2, dlogbeta, plogbeta, varlogbeta, eslogbeta
  public :: dcompbeta, pcompbeta, varcompbeta, escompbeta, dlnbeta, plnbeta
  public :: varlnbeta, eslnbeta, dmrbeta, pmrbeta, varmrbeta, esmrbeta
contains
  pure elemental function dloggamma(x, a, r, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, r
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: r_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    r_v = 1.0_dp
    if (present(r)) r_v = r
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((1.0_dp / x) * gamma_pdf((-log(x)), shape=r_v, rate=a_v))
    if (((log_pdf_v))) res = (gamma_pdf((-log(x)), shape=r_v, rate=a_v, log_pdf=.true.) - log(x))
  end function dloggamma

  pure elemental function ploggamma(x, a, r, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, r
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: r_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    r_v = 1.0_dp
    if (present(r)) r_v = r
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - gamma_cdf(((-a_v) * log(x)), shape=r_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = gamma_cdf(((-a_v) * log(x)), shape=r_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - gamma_cdf(((-a_v) * log(x)), shape=r_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = gamma_cdf(((-a_v) * log(x)), shape=r_v, log_p=.true.)
  end function ploggamma

  pure elemental function varloggamma(p, a, r, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, r
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: r_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    r_v = 1.0_dp
    if (present(r)) r_v = r
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = exp(((-(1.0_dp / a_v)) * gamma_quantile((1.0_dp - pp), shape=r_v)))
  end function varloggamma

  pure elemental function esloggamma(p, a, r) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, r
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: r_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    r_v = 1.0_dp
    if (present(r)) r_v = r
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varloggamma(p * s_quad**3, a_v, r_v)
    end do
  end function esloggamma

  pure elemental function dinvgamma(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((1.0_dp / x) * gamma_pdf((1.0_dp / x), shape=a_v, rate=b_v))
    if (((log_pdf_v))) res = (gamma_pdf((1.0_dp / x), shape=a_v, rate=b_v, log_pdf=.true.) - log(x))
  end function dinvgamma

  pure elemental function pinvgamma(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - gamma_cdf((b_v / x), shape=a_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = gamma_cdf((b_v / x), shape=a_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - gamma_cdf((b_v / x), shape=a_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = gamma_cdf((b_v / x), shape=a_v, log_p=.true.)
  end function pinvgamma

  pure elemental function varinvgamma(p, a, b, log_p, lower_tail) result(res)
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
    res = (b_v / gamma_quantile((1.0_dp - pp), shape=a_v))
  end function varinvgamma

  pure elemental function esinvgamma(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varinvgamma(p * s_quad**3, a_v, b_v)
    end do
  end function esinvgamma

  pure elemental function dstacygamma(x, gamma, c, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: gamma, c, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: c_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    c_v = 1.0_dp
    if (present(c)) c_v = c
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((c_v * (x ** (c_v - 1.0_dp))) / (theta_v ** c_v)) * gamma_pdf(((x / theta_v) &
      & ** c_v), shape=gamma_v))
    if (((log_pdf_v))) res = (((log(c_v) + ((c_v - 1.0_dp) * log(x))) - (c_v * log(theta_v))) + gamma_pdf(((x / &
      & theta_v) ** c_v), shape=gamma_v, log_pdf=.true.))
  end function dstacygamma

  pure elemental function pstacygamma(x, gamma, c, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: gamma, c, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: c_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    c_v = 1.0_dp
    if (present(c)) c_v = c
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = gamma_cdf(((x / theta_v) ** c_v), shape=gamma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pstacygamma

  pure elemental function varstacygamma(p, gamma, c, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: gamma, c, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: c_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    c_v = 1.0_dp
    if (present(c)) c_v = c
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (theta_v * (gamma_quantile(pp, shape=gamma_v) ** (1.0_dp / c_v)))
  end function varstacygamma

  pure elemental function esstacygamma(p, gamma, c, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: gamma, c, theta
    real(dp) :: res
    real(dp) :: gamma_v
    real(dp) :: c_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    gamma_v = 1.0_dp
    if (present(gamma)) gamma_v = gamma
    c_v = 1.0_dp
    if (present(c)) c_v = c
    theta_v = 1.0_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varstacygamma(p * s_quad**3, gamma_v, c_v, theta_v)
    end do
  end function esstacygamma

  pure elemental function dbetadist(x, a, b, log_pdf) result(res)
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
    res = beta_pdf(x, shape1=a_v, shape2=b_v, log_pdf=log_pdf_v)
  end function dbetadist

  pure elemental function pbetadist(x, a, b, log_p, lower_tail) result(res)
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
    res = beta_cdf(x, shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pbetadist

  pure elemental function varbetadist(p, a, b, log_p, lower_tail) result(res)
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
    res = beta_quantile(pp, shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varbetadist

  pure elemental function esbetadist(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetadist(p * s_quad**3, a_v, b_v)
    end do
  end function esbetadist

  pure elemental function duniform(x, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = uniform_pdf(x, lower=a_v, upper=b_v, log_pdf=log_pdf_v)
  end function duniform

  pure elemental function puniform(x, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = uniform_cdf(x, lower=a_v, upper=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function puniform

  pure elemental function varuniform(p, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = uniform_quantile(pp, lower=a_v, upper=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varuniform

  pure elemental function esuniform(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varuniform(p * s_quad**3, a_v, b_v)
    end do
  end function esuniform

  pure elemental function dgenunif(x, a, c, h, k, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, c, h, k
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: h_v
    real(dp) :: k_v
    logical :: log_pdf_v
    a_v = 0.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    h_v = 1.0_dp
    if (present(h)) h_v = h
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((h_v * k_v) * c_v) * ((x - a_v) ** (c_v - 1.0_dp))) * ((1.0_dp - (k_v * ((x - &
      & a_v) ** c_v))) ** (h_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((log(((h_v * k_v) * c_v)) + ((c_v - 1.0_dp) * log((x - a_v)))) + ((h_v - 1.0_dp) * &
      & log((1.0_dp - (k_v * ((x - a_v) ** c_v))))))
  end function dgenunif

  pure elemental function pgenunif(x, a, c, h, k, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, c, h, k
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: h_v
    real(dp) :: k_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    h_v = 1.0_dp
    if (present(h)) h_v = h
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - (k_v * ((x - a_v) ** c_v))) ** h_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - (k_v * ((x - a_v) ** c_v))) ** h_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - (k_v * ((x - a_v) ** c_v))) ** h_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (h_v * log((1.0_dp - (k_v * ((x - a_v) ** c_v)))))
  end function pgenunif

  pure elemental function vargenunif(p, a, c, h, k, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, c, h, k
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: h_v
    real(dp) :: k_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 0.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    h_v = 1.0_dp
    if (present(h)) h_v = h
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (a_v + ((k_v ** ((-1.0_dp) / c_v)) * ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / h_v))) ** (1.0_dp / c_v))))
  end function vargenunif

  pure elemental function esgenunif(p, a, c, h, k) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, c, h, k
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: h_v
    real(dp) :: k_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    h_v = 1.0_dp
    if (present(h)) h_v = h
    k_v = 1.0_dp
    if (present(k)) k_v = k
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenunif(p * s_quad**3, a_v, c_v, h_v, k_v)
    end do
  end function esgenunif

  pure elemental function dpower1(x, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = beta_pdf(x, shape1=a_v, shape2=1.0_dp, log_pdf=log_pdf_v)
  end function dpower1

  pure elemental function ppower1(x, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(x, shape1=a_v, shape2=1.0_dp, log_p=log_p_v, lower_tail=lower_tail_v)
  end function ppower1

  pure elemental function varpower1(p, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = beta_quantile(pp, shape1=a_v, shape2=1.0_dp, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varpower1

  pure elemental function espower1(p, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a
    real(dp) :: res
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varpower1(p * s_quad**3, a_v)
    end do
  end function espower1

  pure elemental function dpower2(x, b, log_pdf) result(res)
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
    res = beta_pdf(x, shape1=1.0_dp, shape2=b_v, log_pdf=log_pdf_v)
  end function dpower2

  pure elemental function ppower2(x, b, log_p, lower_tail) result(res)
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
    res = beta_cdf(x, shape1=1.0_dp, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function ppower2

  pure elemental function varpower2(p, b, log_p, lower_tail) result(res)
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
    res = beta_quantile(pp, shape1=1.0_dp, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varpower2

  pure elemental function espower2(p, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varpower2(p * s_quad**3, b_v)
    end do
  end function espower2

  pure elemental function dlogbeta(x, a, b, c, d, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c, d
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 2.0_dp
    if (present(d)) d_v = d
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (beta_pdf((log((x / c_v)) / log((d_v / c_v))), shape1=a_v, shape2=b_v) / (x * &
      & log((d_v / c_v))))
    if (((log_pdf_v))) res = ((beta_pdf((log((x / c_v)) / log((d_v / c_v))), shape1=a_v, shape2=b_v, &
      & log_pdf=.true.) - log(x)) - log(log((d_v / c_v))))
  end function dlogbeta

  pure elemental function plogbeta(x, a, b, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 2.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf((log((x / c_v)) / log((d_v / c_v))), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function plogbeta

  pure elemental function varlogbeta(p, a, b, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 2.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = (c_v * ((d_v / c_v) ** beta_quantile(pp, shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)))
  end function varlogbeta

  pure elemental function eslogbeta(p, a, b, c, d) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c, d
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 2.0_dp
    if (present(d)) d_v = d
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlogbeta(p * s_quad**3, a_v, b_v, c_v, d_v)
    end do
  end function eslogbeta

  pure elemental function dcompbeta(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((beta_fn(a_v, b_v) * (beta_quantile(x, shape1=a_v, shape2=b_v) ** (1.0_dp - &
      & a_v))) * ((1.0_dp - beta_quantile(x, shape1=a_v, shape2=b_v)) ** (1.0_dp - b_v)))
    if (((log_pdf_v))) res = ((log_beta_fn(a_v, b_v) + ((1.0_dp - a_v) * log(beta_quantile(x, shape1=a_v, &
      & shape2=b_v)))) + ((1.0_dp - b_v) * log((1.0_dp - beta_quantile(x, shape1=a_v, shape2=b_v)))))
  end function dcompbeta

  pure elemental function pcompbeta(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = beta_quantile(x, shape1=a_v, shape2=b_v)
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - beta_quantile(x, shape1=a_v, shape2=b_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log(beta_quantile(x, shape1=a_v, shape2=b_v))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - beta_quantile(x, shape1=a_v, shape2=b_v)))
  end function pcompbeta

  pure elemental function varcompbeta(p, a, b, log_p, lower_tail) result(res)
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
    res = beta_cdf(pp, shape1=a_v, shape2=b_v)
  end function varcompbeta

  pure elemental function escompbeta(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varcompbeta(p * s_quad**3, a_v, b_v)
    end do
  end function escompbeta

  pure elemental function dlnbeta(x, lambda, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((lambda_v * ((1.0_dp + ((lambda_v - 1.0_dp) * x)) ** (-2.0_dp))) * &
      & beta_pdf(((lambda_v * x) / (1.0_dp + ((lambda_v - 1.0_dp) * x))), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = ((log(lambda_v) - (2.0_dp * log((1.0_dp + ((lambda_v - 1.0_dp) * x))))) + &
      & beta_pdf(((lambda_v * x) / (1.0_dp + ((lambda_v - 1.0_dp) * x))), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dlnbeta

  pure elemental function plnbeta(x, lambda, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(((lambda_v * x) / (1.0_dp + ((lambda_v - 1.0_dp) * x))), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function plnbeta

  pure elemental function varlnbeta(p, lambda, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
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
    res = (beta_quantile(pp, shape1=a_v, shape2=b_v) / (lambda_v - ((lambda_v - 1.0_dp) * beta_quantile(pp, &
      & shape1=a_v, shape2=b_v))))
  end function varlnbeta

  pure elemental function eslnbeta(p, lambda, a, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a, b
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlnbeta(p * s_quad**3, lambda_v, a_v, b_v)
    end do
  end function eslnbeta

  pure elemental function dmrbeta(x, a, b, r, q, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, r, q
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: r_v
    real(dp) :: q_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    r_v = 1.0_dp
    if (present(r)) r_v = r
    q_v = 1.0_dp
    if (present(q)) q_v = q
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((r_v / b_v) * ((x ** (r_v - 1.0_dp)) / (q_v ** r_v))) * beta_pdf(((1.0_dp / &
      & b_v) * ((x / q_v) ** r_v)), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = ((((log(r_v) - log(b_v)) - (r_v * log(q_v))) + ((r_v - 1.0_dp) * log(x))) + &
      & beta_pdf(((1.0_dp / b_v) * ((x / q_v) ** r_v)), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dmrbeta

  pure elemental function pmrbeta(x, a, b, r, q, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, r, q
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: r_v
    real(dp) :: q_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    r_v = 1.0_dp
    if (present(r)) r_v = r
    q_v = 1.0_dp
    if (present(q)) q_v = q
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(((1.0_dp / b_v) * ((x / q_v) ** r_v)), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pmrbeta

  pure elemental function varmrbeta(p, a, b, r, q, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, r, q
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: r_v
    real(dp) :: q_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    r_v = 1.0_dp
    if (present(r)) r_v = r
    q_v = 1.0_dp
    if (present(q)) q_v = q
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (((b_v ** (1.0_dp / r_v)) * q_v) * (beta_quantile(pp, shape1=a_v, shape2=b_v) ** (1.0_dp / r_v)))
  end function varmrbeta

  pure elemental function esmrbeta(p, a, b, r, q) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, r, q
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: r_v
    real(dp) :: q_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    r_v = 1.0_dp
    if (present(r)) r_v = r
    q_v = 1.0_dp
    if (present(q)) q_v = q
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varmrbeta(p * s_quad**3, a_v, b_v, r_v, q_v)
    end do
  end function esmrbeta

end module vares_distributions_03
