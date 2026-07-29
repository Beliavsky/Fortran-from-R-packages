! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_08
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dgumbel, pgumbel, vargumbel, esgumbel, dkumgumbel, pkumgumbel
  public :: varkumgumbel, eskumgumbel, dbetagumbel, pbetagumbel, varbetagumbel, esbetagumbel
  public :: dgumbel2, pgumbel2, vargumbel2, esgumbel2, dbetagumbel2, pbetagumbel2
  public :: varbetagumbel2, esbetagumbel2, dfrechet, pfrechet, varfrechet, esfrechet
  public :: dbetafrechet, pbetafrechet, varbetafrechet, esbetafrechet, dweibull, pweibull
  public :: varweibull, esweibull, dkumweibull, pkumweibull, varkumweibull, eskumweibull
  public :: dlogisrayleigh, plogisrayleigh, varlogisrayleigh, eslogisrayleigh, dmoweibull, pmoweibull
  public :: varmoweibull, esmoweibull, dbetaweibull, pbetaweibull, varbetaweibull, esbetaweibull
contains
  pure elemental function dgumbel(x, mu, sigma, log_pdf) result(res)
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
    res = x
    if (((.not. log_pdf_v))) res = (((1.0_dp / sigma_v) * exp(((-(x - mu_v)) / sigma_v))) * exp((-exp(((-(x - &
      & mu_v)) / sigma_v)))))
    if (((log_pdf_v))) res = (((-log(sigma_v)) - ((x - mu_v) / sigma_v)) - exp(((-(x - mu_v)) / sigma_v)))
  end function dgumbel

  pure elemental function pgumbel(x, mu, sigma, log_p, lower_tail) result(res)
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
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = exp((-exp(((-(x - mu_v)) / sigma_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - exp((-exp(((-(x - mu_v)) / sigma_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (-exp(((-(x - mu_v)) / sigma_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - exp((-exp(((-(x - mu_v)) / sigma_v))))))
  end function pgumbel

  pure elemental function vargumbel(p, mu, sigma, log_p, lower_tail) result(res)
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
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (mu_v - (sigma_v * log((-log(pp)))))
  end function vargumbel

  pure elemental function esgumbel(p, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargumbel(p * s_quad**3, mu_v, sigma_v)
    end do
  end function esgumbel

  pure elemental function dkumgumbel(x, a, b, mu, sigma, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((((a_v * b_v) / sigma_v) * exp(((-(x - mu_v)) / sigma_v))) * exp(((-a_v) * &
      & exp(((-(x - mu_v)) / sigma_v))))) * ((1.0_dp - exp(((-a_v) * exp(((-(x - mu_v)) / sigma_v))))) ** (b_v - &
      & 1.0_dp)))
    if (((log_pdf_v))) res = ((((-log(((a_v * b_v) / sigma_v))) - ((x - mu_v) / sigma_v)) - (a_v * exp(((-(x - &
      & mu_v)) / sigma_v)))) + ((b_v - 1.0_dp) * log((1.0_dp - exp(((-a_v) * exp(((-(x - mu_v)) / sigma_v))))))))
  end function dkumgumbel

  pure elemental function pkumgumbel(x, a, b, mu, sigma, log_p, lower_tail) result(res)
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
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - exp(((-a_v) * exp(((-(x - mu_v)) / &
      & sigma_v))))) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - exp(((-a_v) * exp(((-(x - mu_v)) / &
      & sigma_v))))) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - exp(((-a_v) * exp(((-(x - mu_v)) / &
      & sigma_v))))) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - exp(((-a_v) * exp(((-(x - mu_v)) / &
      & sigma_v)))))))
  end function pkumgumbel

  pure elemental function varkumgumbel(p, a, b, mu, sigma, log_p, lower_tail) result(res)
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
    res = (mu_v - (sigma_v * log(((-(1.0_dp / a_v)) * log((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))))))))
  end function varkumgumbel

  pure elemental function eskumgumbel(p, a, b, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumgumbel(p * s_quad**3, a_v, b_v, mu_v, sigma_v)
    end do
  end function eskumgumbel

  pure elemental function dbetagumbel(x, a, b, mu, sigma, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((((1.0_dp / sigma_v) * exp(((-(x - mu_v)) / sigma_v))) * exp((-exp(((-(x - &
      & mu_v)) / sigma_v))))) * beta_pdf(exp((-exp(((-(x - mu_v)) / sigma_v)))), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = ((((-log(sigma_v)) - ((x - mu_v) / sigma_v)) - exp(((-(x - mu_v)) / sigma_v))) + &
      & beta_pdf(exp((-exp(((-(x - mu_v)) / sigma_v)))), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dbetagumbel

  pure elemental function pbetagumbel(x, a, b, mu, sigma, log_p, lower_tail) result(res)
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
    res = beta_cdf(exp((-exp(((-(x - mu_v)) / sigma_v)))), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetagumbel

  pure elemental function varbetagumbel(p, a, b, mu, sigma, log_p, lower_tail) result(res)
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
    res = (mu_v - (sigma_v * log((-log(beta_quantile(pp, shape1=a_v, shape2=b_v))))))
  end function varbetagumbel

  pure elemental function esbetagumbel(p, a, b, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetagumbel(p * s_quad**3, a_v, b_v, mu_v, sigma_v)
    end do
  end function esbetagumbel

  pure elemental function dgumbel2(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((a_v * b_v) * (x ** ((-a_v) - 1.0_dp))) * exp(((-b_v) * (x ** (-a_v)))))
    if (((log_pdf_v))) res = ((log((a_v * b_v)) - ((a_v + 1.0_dp) * log(x))) - (b_v * (x ** (-a_v))))
  end function dgumbel2

  pure elemental function pgumbel2(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp(((-b_v) * (x ** (-a_v)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp(((-b_v) * (x ** (-a_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp(((-b_v) * (x ** (-a_v))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((-b_v) * (x ** (-a_v)))
  end function pgumbel2

  pure elemental function vargumbel2(p, a, b, log_p, lower_tail) result(res)
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
    res = ((b_v ** (1.0_dp / a_v)) * ((-log((1.0_dp - pp))) ** ((-1.0_dp) / a_v)))
  end function vargumbel2

  pure elemental function esgumbel2(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargumbel2(p * s_quad**3, a_v, b_v)
    end do
  end function esgumbel2

  pure elemental function dbetagumbel2(x, a, b, c, d, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((((a_v * b_v) * (x ** ((-a_v) - 1.0_dp))) * exp(((-b_v) * (x ** (-a_v))))) * &
      & beta_pdf((1.0_dp - exp(((-b_v) * (x ** (-a_v))))), shape1=c_v, shape2=d_v))
    if (((log_pdf_v))) res = (((log((a_v * b_v)) - ((a_v + 1.0_dp) * log(x))) - (b_v * (x ** (-a_v)))) + &
      & beta_pdf((1.0_dp - exp(((-b_v) * (x ** (-a_v))))), shape1=c_v, shape2=d_v, log_pdf=.true.))
  end function dbetagumbel2

  pure elemental function pbetagumbel2(x, a, b, c, d, log_p, lower_tail) result(res)
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
    res = beta_cdf((1.0_dp - exp(((-b_v) * (x ** (-a_v))))), shape1=c_v, shape2=d_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetagumbel2

  pure elemental function varbetagumbel2(p, a, b, c, d, log_p, lower_tail) result(res)
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
    res = ((b_v ** (1.0_dp / a_v)) * ((-log((1.0_dp - beta_quantile(pp, shape1=c_v, shape2=d_v)))) ** ((-1.0_dp) / &
      & a_v)))
  end function varbetagumbel2

  pure elemental function esbetagumbel2(p, a, b, c, d) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetagumbel2(p * s_quad**3, a_v, b_v, c_v, d_v)
    end do
  end function esbetagumbel2

  pure elemental function dfrechet(x, alpha, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((alpha_v * (sigma_v ** alpha_v)) * (x ** ((-alpha_v) - 1.0_dp))) * &
      & exp((-((sigma_v / x) ** alpha_v))))
    if (((log_pdf_v))) res = (((log(alpha_v) + (alpha_v * log(sigma_v))) - ((alpha_v + 1.0_dp) * log(x))) - &
      & ((sigma_v / x) ** alpha_v))
  end function dfrechet

  pure elemental function pfrechet(x, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = exp((-((sigma_v / x) ** alpha_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - exp((-((sigma_v / x) ** alpha_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (-((sigma_v / x) ** alpha_v))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - exp((-((sigma_v / x) ** alpha_v)))))
  end function pfrechet

  pure elemental function varfrechet(p, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
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
    res = (sigma_v * ((-log(pp)) ** ((-1.0_dp) / alpha_v)))
  end function varfrechet

  pure elemental function esfrechet(p, alpha, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha, sigma
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varfrechet(p * s_quad**3, alpha_v, sigma_v)
    end do
  end function esfrechet

  pure elemental function dbetafrechet(x, a, b, alpha, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((alpha_v * (sigma_v ** alpha_v)) * (x ** ((-alpha_v) - 1.0_dp))) * &
      & exp((-((sigma_v / x) ** alpha_v)))) * beta_pdf(exp((-((sigma_v / x) ** alpha_v))), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = ((((log(alpha_v) + (alpha_v * log(sigma_v))) - ((alpha_v + 1.0_dp) * log(x))) - &
      & ((sigma_v / x) ** alpha_v)) + beta_pdf(exp((-((sigma_v / x) ** alpha_v))), shape1=a_v, shape2=b_v, &
      & log_pdf=.true.))
  end function dbetafrechet

  pure elemental function pbetafrechet(x, a, b, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(exp((-((sigma_v / x) ** alpha_v))), shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pbetafrechet

  pure elemental function varbetafrechet(p, a, b, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
    res = (sigma_v * ((-log(beta_quantile(pp, shape1=a_v, shape2=b_v))) ** ((-1.0_dp) / alpha_v)))
  end function varbetafrechet

  pure elemental function esbetafrechet(p, a, b, alpha, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetafrechet(p * s_quad**3, a_v, b_v, alpha_v, sigma_v)
    end do
  end function esbetafrechet

  pure elemental function dweibull(x, alpha, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = weibull_pdf(x, shape=alpha_v, scale=sigma_v, log_pdf=log_pdf_v)
  end function dweibull

  pure elemental function pweibull(x, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = weibull_cdf(x, shape=alpha_v, scale=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pweibull

  pure elemental function varweibull(p, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = weibull_quantile(pp, shape=alpha_v, scale=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varweibull

  pure elemental function esweibull(p, alpha, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: alpha, sigma
    real(dp) :: res
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varweibull(p * s_quad**3, alpha_v, sigma_v)
    end do
  end function esweibull

  pure elemental function dkumweibull(x, a, b, alpha, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * b_v) * weibull_pdf(x, shape=alpha_v, scale=sigma_v)) * (weibull_cdf(x, &
      & shape=alpha_v, scale=sigma_v) ** (a_v - 1.0_dp))) * ((1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v)) ** (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = (((log((a_v * b_v)) + weibull_pdf(x, shape=alpha_v, scale=sigma_v, log_pdf=.true.)) + &
      & ((a_v - 1.0_dp) * log(weibull_cdf(x, shape=alpha_v, scale=sigma_v)))) + ((b_v - 1.0_dp) * log((1.0_dp - &
      & (weibull_cdf(x, shape=alpha_v, scale=sigma_v) ** a_v)))))
  end function dkumweibull

  pure elemental function pkumweibull(x, a, b, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - (weibull_cdf(x, shape=alpha_v, &
      & scale=sigma_v) ** a_v))))
  end function pkumweibull

  pure elemental function varkumweibull(p, a, b, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
    res = (sigma_v * ((-log((1.0_dp - ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v))))) ** &
      & (1.0_dp / alpha_v)))
  end function varkumweibull

  pure elemental function eskumweibull(p, a, b, alpha, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumweibull(p * s_quad**3, a_v, b_v, alpha_v, sigma_v)
    end do
  end function eskumweibull

  pure elemental function dlogisrayleigh(x, a, lambda, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((((a_v * lambda_v) * x) * exp((((lambda_v * x) * x) / 2.0_dp))) * &
      & ((exp((((lambda_v * x) * x) / 2.0_dp)) - 1.0_dp) ** (a_v - 1.0_dp))) * ((1.0_dp + ((exp((((lambda_v * x) * &
      & x) / 2.0_dp)) - 1.0_dp) ** a_v)) ** (-2.0_dp)))
    if (((log_pdf_v))) res = ((((log((a_v * lambda_v)) + log(x)) + (((lambda_v * x) * x) / 2.0_dp)) + ((a_v - &
      & 1.0_dp) * log((exp((((lambda_v * x) * x) / 2.0_dp)) - 1.0_dp)))) - (2.0_dp * log((1.0_dp + &
      & ((exp((((lambda_v * x) * x) / 2.0_dp)) - 1.0_dp) ** a_v)))))
  end function dlogisrayleigh

  pure elemental function plogisrayleigh(x, a, lambda, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - (1.0_dp / (1.0_dp + ((exp((((lambda_v * x) * x) &
      & / 2.0_dp)) - 1.0_dp) ** a_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp / (1.0_dp + ((exp((((lambda_v * x) * x) / &
      & 2.0_dp)) - 1.0_dp) ** a_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - (1.0_dp / (1.0_dp + ((exp((((lambda_v * x) * x) / &
      & 2.0_dp)) - 1.0_dp) ** a_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (-log((1.0_dp + ((exp((((lambda_v * x) * x) / 2.0_dp)) - &
      & 1.0_dp) ** a_v))))
  end function plogisrayleigh

  pure elemental function varlogisrayleigh(p, a, lambda, log_p, lower_tail) result(res)
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
    res = (sqrt((2.0_dp / lambda_v)) * sqrt(log((1.0_dp + ((pp / (1.0_dp - pp)) ** (1.0_dp / a_v))))))
  end function varlogisrayleigh

  pure elemental function eslogisrayleigh(p, a, lambda) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlogisrayleigh(p * s_quad**3, a_v, lambda_v)
    end do
  end function eslogisrayleigh

  pure elemental function dmoweibull(x, a, b, lambda, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((((b_v * (lambda_v ** b_v)) * (x ** (b_v - 1.0_dp))) * exp(((lambda_v * x) ** &
      & b_v))) * (((exp(((lambda_v * x) ** b_v)) - 1.0_dp) + a_v) ** (-2.0_dp)))
    if (((log_pdf_v))) res = ((((log(b_v) + (b_v * log(lambda_v))) + ((b_v - 1.0_dp) * log(x))) + ((lambda_v * x) &
      & ** b_v)) - (2.0_dp * log(((exp(((lambda_v * x) ** b_v)) - 1.0_dp) + a_v))))
  end function dmoweibull

  pure elemental function pmoweibull(x, a, b, lambda, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - (1.0_dp / ((exp(((lambda_v * x) ** b_v)) - &
      & 1.0_dp) + a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp / ((exp(((lambda_v * x) ** b_v)) - 1.0_dp) &
      & + a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - (1.0_dp / ((exp(((lambda_v * x) ** b_v)) - &
      & 1.0_dp) + a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (-log(((exp(((lambda_v * x) ** b_v)) - 1.0_dp) + a_v)))
  end function pmoweibull

  pure elemental function varmoweibull(p, a, b, lambda, log_p, lower_tail) result(res)
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
    res = ((1.0_dp / lambda_v) * (log((((1.0_dp / (1.0_dp - pp)) + 1.0_dp) - a_v)) ** (1.0_dp / b_v)))
  end function varmoweibull

  pure elemental function esmoweibull(p, a, b, lambda) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varmoweibull(p * s_quad**3, a_v, b_v, lambda_v)
    end do
  end function esmoweibull

  pure elemental function dbetaweibull(x, a, b, alpha, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (weibull_pdf(x, shape=alpha_v, scale=sigma_v) * beta_pdf(weibull_cdf(x, &
      & shape=alpha_v, scale=sigma_v), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = (weibull_pdf(x, shape=alpha_v, scale=sigma_v, log_pdf=.true.) + &
      & beta_pdf(weibull_cdf(x, shape=alpha_v, scale=sigma_v), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dbetaweibull

  pure elemental function pbetaweibull(x, a, b, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    alpha_v = 1.0_dp
    if (present(alpha)) alpha_v = alpha
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(weibull_cdf(x, shape=alpha_v, scale=sigma_v), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetaweibull

  pure elemental function varbetaweibull(p, a, b, alpha, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
    res = (sigma_v * ((-log((1.0_dp - beta_quantile(pp, shape1=a_v, shape2=b_v)))) ** (1.0_dp / alpha_v)))
  end function varbetaweibull

  pure elemental function esbetaweibull(p, a, b, alpha, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, alpha, sigma
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: alpha_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetaweibull(p * s_quad**3, a_v, b_v, alpha_v, sigma_v)
    end do
  end function esbetaweibull

end module vares_distributions_08
