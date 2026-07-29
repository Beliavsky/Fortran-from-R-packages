! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_07
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dexplogis, pexplogis, varexplogis, esexplogis, dhlogis, phlogis
  public :: varhlogis, eshlogis, dlognorm, plognorm, varlognorm, eslognorm
  public :: dbetalognorm, pbetalognorm, varbetalognorm, esbetalognorm, dburr, pburr
  public :: varburr, esburr, dbetaburr, pbetaburr, varbetaburr, esbetaburr
  public :: dburr7, pburr7, varburr7, esburr7, dkumburr7, pkumburr7
  public :: varkumburr7, eskumburr7, dbetaburr7, pbetaburr7, varbetaburr7, esbetaburr7
  public :: ddagum, pdagum, vardagum, esdagum, dlomax, plomax
  public :: varlomax, eslomax, dbetalomax, pbetalomax, varbetalomax, esbetalomax
contains
  pure elemental function dexplogis(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((a_v / b_v) * exp(((-x) / b_v))) * ((1.0_dp + exp(((-x) / b_v))) ** ((-a_v) - &
      & 1.0_dp)))
    if (((log_pdf_v))) res = (((log(a_v) - log(b_v)) - (x / b_v)) - ((a_v + 1.0_dp) * log((1.0_dp + exp(((-x) / &
      & b_v))))))
  end function dexplogis

  pure elemental function pexplogis(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp + exp(((-x) / b_v))) ** (-a_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp + exp(((-x) / b_v))) ** (-a_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((-a_v) * log((1.0_dp + exp(((-x) / b_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((1.0_dp + exp(((-x) / b_v))) ** (-a_v))))
  end function pexplogis

  pure elemental function varexplogis(p, a, b, log_p, lower_tail) result(res)
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
    res = ((-b_v) * log(((pp ** ((-1.0_dp) / a_v)) - 1.0_dp)))
  end function varexplogis

  pure elemental function esexplogis(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexplogis(p * s_quad**3, a_v, b_v)
    end do
  end function esexplogis

  pure elemental function dhlogis(x, k, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: k_v
    logical :: log_pdf_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((1.0_dp - (k_v * x)) ** ((1.0_dp / k_v) - 1.0_dp)) * ((1.0_dp + ((1.0_dp - &
      & (k_v * x)) ** (1.0_dp / k_v))) ** (-2.0_dp)))
    if (((log_pdf_v))) res = ((((1.0_dp / k_v) - 1.0_dp) * log((1.0_dp - (k_v * x)))) - (2.0_dp * log((1.0_dp + &
      & ((1.0_dp - (k_v * x)) ** (1.0_dp / k_v))))))
  end function dhlogis

  pure elemental function phlogis(x, k, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    logical :: log_p_v
    logical :: lower_tail_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp + ((1.0_dp - (k_v * x)) ** (1.0_dp / k_v))) ** &
      & (-1.0_dp))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (((1.0_dp - (k_v * x)) ** (1.0_dp / k_v)) * &
      & ((1.0_dp + ((1.0_dp - (k_v * x)) ** (1.0_dp / k_v))) ** (-1.0_dp)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (-log((1.0_dp + ((1.0_dp - (k_v * x)) ** (1.0_dp / k_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (((1.0_dp / k_v) * log((1.0_dp - (k_v * x)))) - &
      & log((1.0_dp + ((1.0_dp - (k_v * x)) ** (1.0_dp / k_v)))))
  end function phlogis

  pure elemental function varhlogis(p, k, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / k_v) * (1.0_dp - (((1.0_dp / pp) - 1.0_dp) ** k_v)))
  end function varhlogis

  pure elemental function eshlogis(p, k) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k
    real(dp) :: res
    real(dp) :: k_v
    integer :: i
    real(dp) :: s_quad
    k_v = 1.0_dp
    if (present(k)) k_v = k
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhlogis(p * s_quad**3, k_v)
    end do
  end function eshlogis

  pure elemental function dlognorm(x, mu, sigma, log_pdf) result(res)
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
    res = lognormal_pdf(x, meanlog=mu_v, sdlog=sigma_v, log_pdf=log_pdf_v)
  end function dlognorm

  pure elemental function plognorm(x, mu, sigma, log_p, lower_tail) result(res)
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
    res = lognormal_cdf(x, meanlog=mu_v, sdlog=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function plognorm

  pure elemental function varlognorm(p, mu, sigma, log_p, lower_tail) result(res)
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
    res = lognormal_quantile(pp, meanlog=mu_v, sdlog=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varlognorm

  pure elemental function eslognorm(p, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlognorm(p * s_quad**3, mu_v, sigma_v)
    end do
  end function eslognorm

  pure elemental function dbetalognorm(x, a, b, mu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, mu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (lognormal_pdf(x, meanlog=mu_v, sdlog=sigma_v) * beta_pdf(lognormal_cdf(x, &
      & meanlog=mu_v, sdlog=sigma_v), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = (lognormal_pdf(x, meanlog=mu_v, sdlog=sigma_v, log_pdf=.true.) + &
      & beta_pdf(lognormal_cdf(x, meanlog=mu_v, sdlog=sigma_v), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dbetalognorm

  pure elemental function pbetalognorm(x, a, b, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(lognormal_cdf(x, meanlog=mu_v, sdlog=sigma_v), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetalognorm

  pure elemental function varbetalognorm(p, a, b, mu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, mu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
    res = exp((mu_v + (sigma_v * normal_quantile(beta_quantile(pp, shape1=a_v, shape2=b_v)))))
  end function varbetalognorm

  pure elemental function esbetalognorm(p, a, b, mu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, mu, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: mu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetalognorm(p * s_quad**3, a_v, b_v, mu_v, sigma_v)
    end do
  end function esbetalognorm

  pure elemental function dburr(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((b_v * (a_v ** b_v)) * (x ** ((-b_v) - 1.0_dp))) * ((1.0_dp + ((x / a_v) ** &
      & (-b_v))) ** (-2.0_dp)))
    if (((log_pdf_v))) res = (((log(b_v) + (b_v * log(a_v))) - ((b_v + 1.0_dp) * log(x))) - (2.0_dp * log((1.0_dp + &
      & ((x / a_v) ** (-b_v))))))
  end function dburr

  pure elemental function pburr(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp + ((x / a_v) ** (-b_v))) ** (-1.0_dp))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp + ((x / a_v) ** b_v)) ** (-1.0_dp))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (-log((1.0_dp + ((x / a_v) ** (-b_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (-log((1.0_dp + ((x / a_v) ** b_v))))
  end function pburr

  pure elemental function varburr(p, a, b, log_p, lower_tail) result(res)
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
    res = ((a_v * (pp ** (1.0_dp / b_v))) * ((1.0_dp - pp) ** ((-1.0_dp) / b_v)))
  end function varburr

  pure elemental function esburr(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varburr(p * s_quad**3, a_v, b_v)
    end do
  end function esburr

  pure elemental function dbetaburr(x, a, b, c, d, log_pdf) result(res)
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
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((b_v * (a_v ** b_v)) * (x ** ((-b_v) - 1.0_dp))) * ((1.0_dp + ((x / a_v) ** &
      & (-b_v))) ** (-2.0_dp))) * beta_pdf(((1.0_dp + ((x / a_v) ** (-b_v))) ** (-1.0_dp)), shape1=c_v, shape2=d_v))
    if (((log_pdf_v))) res = ((((log(b_v) + (b_v * log(a_v))) - ((b_v + 1.0_dp) * log(x))) - (2.0_dp * log((1.0_dp &
      & + ((x / a_v) ** (-b_v)))))) + beta_pdf(((1.0_dp + ((x / a_v) ** (-b_v))) ** (-1.0_dp)), shape1=c_v, &
      & shape2=d_v, log_pdf=.true.))
  end function dbetaburr

  pure elemental function pbetaburr(x, a, b, c, d, log_p, lower_tail) result(res)
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
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    res = beta_cdf(((1.0_dp + ((x / a_v) ** (-b_v))) ** (-1.0_dp)), shape1=c_v, shape2=d_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetaburr

  pure elemental function varbetaburr(p, a, b, c, d, log_p, lower_tail) result(res)
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
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (a_v * ((beta_quantile(pp, shape1=c_v, shape2=d_v) / (1.0_dp - beta_quantile(pp, shape1=c_v, &
      & shape2=d_v))) ** (1.0_dp / b_v)))
  end function varbetaburr

  pure elemental function esbetaburr(p, a, b, c, d) result(res)
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
    d_v = 1.0_dp
    if (present(d)) d_v = d
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetaburr(p * s_quad**3, a_v, b_v, c_v, d_v)
    end do
  end function esbetaburr

  pure elemental function dburr7(x, k, c, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k, c
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: c_v
    logical :: log_pdf_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((k_v * c_v) * (x ** (c_v - 1.0_dp))) * ((1.0_dp + (x ** c_v)) ** ((-k_v) - &
      & 1.0_dp)))
    if (((log_pdf_v))) res = ((log((k_v * c_v)) + ((c_v - 1.0_dp) * log(x))) - ((k_v + 1.0_dp) * log((1.0_dp + (x &
      & ** c_v)))))
  end function dburr7

  pure elemental function pburr7(x, k, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp + (x ** c_v)) ** (-k_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((-k_v) * log((1.0_dp + (x ** c_v))))
  end function pburr7

  pure elemental function varburr7(p, k, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((((1.0_dp - pp) ** ((-1.0_dp) / k_v)) - 1.0_dp) ** (1.0_dp / c_v))
  end function varburr7

  pure elemental function esburr7(p, k, c) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k, c
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: c_v
    integer :: i
    real(dp) :: s_quad
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varburr7(p * s_quad**3, k_v, c_v)
    end do
  end function esburr7

  pure elemental function dkumburr7(x, a, b, k, c, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, k, c
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: k_v
    real(dp) :: c_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((((((a_v * b_v) * k_v) * c_v) * (x ** (c_v - 1.0_dp))) * ((1.0_dp + (x ** &
      & c_v)) ** ((-k_v) - 1.0_dp))) * ((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))) ** (a_v - 1.0_dp))) * ((1.0_dp &
      & - ((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))) ** a_v)) ** (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((((log((((a_v * b_v) * k_v) * c_v)) + ((c_v - 1.0_dp) * log(x))) - ((k_v + 1.0_dp) * &
      & log((1.0_dp + (x ** c_v))))) + ((a_v - 1.0_dp) * log((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v)))))) + &
      & ((b_v - 1.0_dp) * log((1.0_dp - ((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))) ** a_v)))))
  end function dkumburr7

  pure elemental function pkumburr7(x, a, b, k, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, k, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: k_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - ((1.0_dp - ((1.0_dp + (x ** c_v)) ** &
      & (-k_v))) ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - ((1.0_dp - ((1.0_dp + (x ** c_v)) ** &
      & (-k_v))) ** a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - ((1.0_dp - ((1.0_dp + (x ** c_v)) ** &
      & (-k_v))) ** a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - ((1.0_dp - ((1.0_dp + (x ** c_v)) ** &
      & (-k_v))) ** a_v))))
  end function pkumburr7

  pure elemental function varkumburr7(p, a, b, k, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, k, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: k_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((((1.0_dp - ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v))) ** ((-1.0_dp) / k_v)) - &
      & 1.0_dp) ** (1.0_dp / c_v))
  end function varkumburr7

  pure elemental function eskumburr7(p, a, b, k, c) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, k, c
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: k_v
    real(dp) :: c_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    k_v = 1.0_dp
    if (present(k)) k_v = k
    c_v = 1.0_dp
    if (present(c)) c_v = c
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumburr7(p * s_quad**3, a_v, b_v, k_v, c_v)
    end do
  end function eskumburr7

  pure elemental function dbetaburr7(x, a, b, c, k, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c, k
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: k_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((c_v * k_v) * (x ** (c_v - 1.0_dp))) * ((1.0_dp + (x ** c_v)) ** ((-k_v) - &
      & 1.0_dp))) * beta_pdf((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = (((log((c_v * k_v)) + ((c_v - 1.0_dp) * log(x))) - ((k_v + 1.0_dp) * log((1.0_dp + (x &
      & ** c_v))))) + beta_pdf((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dbetaburr7

  pure elemental function pbetaburr7(x, a, b, c, k, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c, k
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: k_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    res = beta_cdf((1.0_dp - ((1.0_dp + (x ** c_v)) ** (-k_v))), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetaburr7

  pure elemental function varbetaburr7(p, a, b, c, k, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c, k
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: k_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    k_v = 1.0_dp
    if (present(k)) k_v = k
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((((1.0_dp - beta_quantile(pp, shape1=a_v, shape2=b_v)) ** ((-1.0_dp) / k_v)) - 1.0_dp) ** (1.0_dp / c_v))
  end function varbetaburr7

  pure elemental function esbetaburr7(p, a, b, c, k) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c, k
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: k_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    k_v = 1.0_dp
    if (present(k)) k_v = k
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetaburr7(p * s_quad**3, a_v, b_v, c_v, k_v)
    end do
  end function esbetaburr7

  pure elemental function ddagum(x, a, b, c, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * c_v) * (b_v ** a_v)) * (x ** ((a_v * c_v) - 1.0_dp))) * (((x ** a_v) + &
      & (b_v ** a_v)) ** ((-c_v) - 1.0_dp)))
    if (((log_pdf_v))) res = (((log((a_v * c_v)) + (a_v * log(b_v))) + (((a_v * c_v) - 1.0_dp) * log(x))) - ((c_v + &
      & 1.0_dp) * log(((x ** a_v) + (b_v ** a_v)))))
  end function ddagum

  pure elemental function pdagum(x, a, b, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
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
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp + ((b_v / x) ** a_v)) ** (-c_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp + ((b_v / x) ** a_v)) ** (-c_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((-c_v) * log((1.0_dp + ((b_v / x) ** a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((1.0_dp + ((b_v / x) ** a_v)) ** (-c_v))))
  end function pdagum

  pure elemental function vardagum(p, a, b, c, log_p, lower_tail) result(res)
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
    res = (b_v * (((pp ** ((-1.0_dp) / c_v)) - 1.0_dp) ** ((-1.0_dp) / a_v)))
  end function vardagum

  pure elemental function esdagum(p, a, b, c) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vardagum(p * s_quad**3, a_v, b_v, c_v)
    end do
  end function esdagum

  pure elemental function dlomax(x, a, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((a_v / lambda_v) * ((1.0_dp + (x / lambda_v)) ** ((-a_v) - 1.0_dp)))
    if (((log_pdf_v))) res = (log((a_v / lambda_v)) - ((a_v + 1.0_dp) * log((1.0_dp + (x / lambda_v)))))
  end function dlomax

  pure elemental function plomax(x, a, lambda, log_p, lower_tail) result(res)
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
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp + (x / lambda_v)) ** (-a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp + (x / lambda_v)) ** (-a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp + (x / lambda_v)) ** (-a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((-a_v) * log((1.0_dp + (x / lambda_v))))
  end function plomax

  pure elemental function varlomax(p, a, lambda, log_p, lower_tail) result(res)
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
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (lambda_v * (((1.0_dp - pp) ** ((-1.0_dp) / a_v)) - 1.0_dp))
  end function varlomax

  pure elemental function eslomax(p, a, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, lambda
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlomax(p * s_quad**3, a_v, lambda_v)
    end do
  end function eslomax

  pure elemental function dbetalomax(x, a, b, alpha, lambda, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, lambda
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: lambda_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((alpha_v / lambda_v) * ((1.0_dp + (x / lambda_v)) ** ((-alpha_v) - 1.0_dp))) * &
      & beta_pdf((1.0_dp - ((1.0_dp + (x / lambda_v)) ** (-alpha_v))), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = ((log((alpha_v / lambda_v)) - ((alpha_v + 1.0_dp) * log((1.0_dp + (x / lambda_v))))) + &
      & beta_pdf((1.0_dp - ((1.0_dp + (x / lambda_v)) ** (-alpha_v))), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dbetalomax

  pure elemental function pbetalomax(x, a, b, alpha, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf((1.0_dp - ((1.0_dp + (x / lambda_v)) ** (-alpha_v))), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetalomax

  pure elemental function varbetalomax(p, a, b, alpha, lambda, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, lambda
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: lambda_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (lambda_v * (((1.0_dp - beta_quantile(pp, shape1=a_v, shape2=b_v)) ** ((-1.0_dp) / alpha_v)) - 1.0_dp))
  end function varbetalomax

  pure elemental function esbetalomax(p, a, b, alpha, lambda) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, lambda
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: lambda_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetalomax(p * s_quad**3, a_v, b_v, alpha_v, lambda_v)
    end do
  end function esbetalomax

end module vares_distributions_07
