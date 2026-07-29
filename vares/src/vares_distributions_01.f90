! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_01
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dexponential, pexponential, varexponential, esexponential, dkumexp, pkumexp
  public :: varkumexp, eskumexp, dexpexp, pexpexp, varexpexp, esexpexp
  public :: dinvexpexp, pinvexpexp, varinvexpexp, esinvexpexp, dbetaexp, pbetaexp
  public :: varbetaexp, esbetaexp, dlogisexp, plogisexp, varlogisexp, eslogisexp
  public :: dexpext, pexpext, varexpext, esexpext, dmoexp, pmoexp
  public :: varmoexp, esmoexp, dperks, pperks, varperks, esperks
  public :: dbeard, pbeard, varbeard, esbeard, dgompertz, pgompertz
  public :: vargompertz, esgompertz, dbetagompertz, pbetagompertz, varbetagompertz, esbetagompertz
contains
  pure elemental function dexponential(x, lambda, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (lambda_v * exp(((-lambda_v) * x)))
    if (((log_pdf_v))) res = (log(lambda_v) - (lambda_v * x))
  end function dexponential

  pure elemental function pexponential(x, lambda, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp(((-lambda_v) * x)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp(((-lambda_v) * x))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp(((-lambda_v) * x))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((-lambda_v) * x)
  end function pexponential

  pure elemental function varexponential(p, lambda, log_p, lower_tail) result(res)
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
    res = (((-1.0_dp) / lambda_v) * log((1.0_dp - pp)))
  end function varexponential

  pure elemental function esexponential(p, lambda) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexponential(p * s_quad**3, lambda_v)
    end do
  end function esexponential

  pure elemental function dkumexp(x, lambda, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((((a_v * b_v) * lambda_v) * exp(((-lambda_v) * x))) * ((1.0_dp - &
      & exp(((-lambda_v) * x))) ** (a_v - 1.0_dp))) * ((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) ** a_v)) ** &
      & (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = (((log(((a_v * b_v) * lambda_v)) - (lambda_v * x)) + ((a_v - 1.0_dp) * log((1.0_dp - &
      & exp(((-lambda_v) * x)))))) + ((b_v - 1.0_dp) * log((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) ** a_v)))))
  end function dkumexp

  pure elemental function pkumexp(x, lambda, a, b, log_p, lower_tail) result(res)
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
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) &
      & ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) ** &
      & a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) ** &
      & a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) &
      & ** a_v))))
  end function pkumexp

  pure elemental function varkumexp(p, lambda, a, b, log_p, lower_tail) result(res)
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
    res = (((-1.0_dp) / lambda_v) * log((1.0_dp - ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v)))))
  end function varkumexp

  pure elemental function eskumexp(p, lambda, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumexp(p * s_quad**3, lambda_v, a_v, b_v)
    end do
  end function eskumexp

  pure elemental function dexpexp(x, lambda, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((a_v * lambda_v) * exp(((-lambda_v) * x))) * ((1.0_dp - exp(((-lambda_v) * &
      & x))) ** (a_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((log((a_v * lambda_v)) - (lambda_v * x)) + ((a_v - 1.0_dp) * log((1.0_dp - &
      & exp(((-lambda_v) * x))))))
  end function dexpexp

  pure elemental function pexpexp(x, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp - exp(((-lambda_v) * x))) ** a_v)
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) ** a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (a_v * log((1.0_dp - exp(((-lambda_v) * x)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - exp(((-lambda_v) * x))) ** a_v)))
  end function pexpexp

  pure elemental function varexpexp(p, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (((-1.0_dp) / lambda_v) * log((1.0_dp - (pp ** (1.0_dp / a_v)))))
  end function varexpexp

  pure elemental function esexpexp(p, lambda, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexpexp(p * s_quad**3, lambda_v, a_v)
    end do
  end function esexpexp

  pure elemental function dinvexpexp(x, lambda, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * lambda_v) * (1.0_dp / (x * x))) * exp(((-lambda_v) / x))) * ((1.0_dp - &
      & exp(((-lambda_v) / x))) ** (a_v - 1.0_dp)))
    if (((log_pdf_v))) res = (((log((a_v * lambda_v)) - (2.0_dp * log(x))) - (lambda_v / x)) + ((a_v - 1.0_dp) * &
      & log((1.0_dp - exp(((-lambda_v) / x))))))
  end function dinvexpexp

  pure elemental function pinvexpexp(x, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - exp(((-lambda_v) / x))) ** a_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - exp(((-lambda_v) / x))) ** a_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - exp(((-lambda_v) / x))) ** a_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (a_v * log((1.0_dp - exp(((-lambda_v) / x)))))
  end function pinvexpexp

  pure elemental function varinvexpexp(p, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (lambda_v / (-log((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / a_v))))))
  end function varinvexpexp

  pure elemental function esinvexpexp(p, lambda, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varinvexpexp(p * s_quad**3, lambda_v, a_v)
    end do
  end function esinvexpexp

  pure elemental function dbetaexp(x, lambda, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((lambda_v * exp((((-b_v) * lambda_v) * x))) * ((1.0_dp - exp(((-lambda_v) * &
      & x))) ** (a_v - 1.0_dp))) / beta_fn(a_v, b_v))
    if (((log_pdf_v))) res = (((log(lambda_v) - ((b_v * lambda_v) * x)) + ((a_v - 1.0_dp) * log((1.0_dp - &
      & exp(((-lambda_v) * x)))))) - log_beta_fn(a_v, b_v))
  end function dbetaexp

  pure elemental function pbetaexp(x, lambda, a, b, log_p, lower_tail) result(res)
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
    res = beta_cdf((1.0_dp - exp(((-lambda_v) * x))), shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pbetaexp

  pure elemental function varbetaexp(p, lambda, a, b, log_p, lower_tail) result(res)
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
    res = (((-1.0_dp) / lambda_v) * log((1.0_dp - beta_quantile(pp, shape1=a_v, shape2=b_v))))
  end function varbetaexp

  pure elemental function esbetaexp(p, lambda, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetaexp(p * s_quad**3, lambda_v, a_v, b_v)
    end do
  end function esbetaexp

  pure elemental function dlogisexp(x, lambda, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * lambda_v) * exp((lambda_v * x))) * ((exp((lambda_v * x)) - 1.0_dp) ** &
      & (a_v - 1.0_dp))) / ((1.0_dp + ((exp((lambda_v * x)) - 1.0_dp) ** a_v)) ** 2.0_dp))
    if (((log_pdf_v))) res = (((log((a_v * lambda_v)) + (lambda_v * x)) + ((a_v - 1.0_dp) * log((exp((lambda_v * &
      & x)) - 1.0_dp)))) - (2.0_dp * log((1.0_dp + ((exp((lambda_v * x)) - 1.0_dp) ** a_v)))))
  end function dlogisexp

  pure elemental function plogisexp(x, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (((exp((lambda_v * x)) - 1.0_dp) ** a_v) / (1.0_dp + &
      & ((exp((lambda_v * x)) - 1.0_dp) ** a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp / (1.0_dp + ((exp((lambda_v * x)) - 1.0_dp) &
      & ** a_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((a_v * log((exp((lambda_v * x)) - 1.0_dp))) - log((1.0_dp + &
      & ((exp((lambda_v * x)) - 1.0_dp) ** a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (-log((1.0_dp + ((exp((lambda_v * x)) - 1.0_dp) ** a_v))))
  end function plogisexp

  pure elemental function varlogisexp(p, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / lambda_v) * log((1.0_dp + ((pp / (1.0_dp - pp)) ** (1.0_dp / a_v)))))
  end function varlogisexp

  pure elemental function eslogisexp(p, lambda, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlogisexp(p * s_quad**3, lambda_v, a_v)
    end do
  end function eslogisexp

  pure elemental function dexpext(x, lambda, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((a_v * lambda_v) * ((1.0_dp + (lambda_v * x)) ** (a_v - 1.0_dp))) * &
      & exp((1.0_dp - ((1.0_dp + (lambda_v * x)) ** a_v))))
    if (((log_pdf_v))) res = (((log((a_v * lambda_v)) + ((a_v - 1.0_dp) * log((1.0_dp + (lambda_v * x))))) + &
      & 1.0_dp) - ((1.0_dp + (lambda_v * x)) ** a_v))
  end function dexpext

  pure elemental function pexpext(x, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp((1.0_dp - ((1.0_dp + (lambda_v * x)) ** &
      & a_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp((1.0_dp - ((1.0_dp + (lambda_v * x)) ** a_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp((1.0_dp - ((1.0_dp + (lambda_v * x)) ** a_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp + (lambda_v * x)) ** a_v))
  end function pexpext

  pure elemental function varexpext(p, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / lambda_v) * (((1.0_dp - log((1.0_dp - pp))) ** (1.0_dp / a_v)) - 1.0_dp))
  end function varexpext

  pure elemental function esexpext(p, lambda, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexpext(p * s_quad**3, lambda_v, a_v)
    end do
  end function esexpext

  pure elemental function dmoexp(x, lambda, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((lambda_v * exp((lambda_v * x))) / (((exp((lambda_v * x)) - 1.0_dp) + a_v) ** &
      & 2.0_dp))
    if (((log_pdf_v))) res = ((log(lambda_v) + (lambda_v * x)) - (2.0_dp * log(((exp((lambda_v * x)) - 1.0_dp) + &
      & a_v))))
  end function dmoexp

  pure elemental function pmoexp(x, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (((exp((lambda_v * x)) - 2.0_dp) + a_v) / ((exp((lambda_v &
      & * x)) - 1.0_dp) + a_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp / ((exp((lambda_v * x)) - 1.0_dp) + a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = (log(((exp((lambda_v * x)) - 2.0_dp) + a_v)) - &
      & log(((exp((lambda_v * x)) - 1.0_dp) + a_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (-log(((exp((lambda_v * x)) - 1.0_dp) + a_v)))
  end function pmoexp

  pure elemental function varmoexp(p, lambda, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / lambda_v) * (log(((2.0_dp - a_v) - ((1.0_dp - a_v) * pp))) - log((1.0_dp - pp))))
  end function varmoexp

  pure elemental function esmoexp(p, lambda, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, a
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varmoexp(p * s_quad**3, lambda_v, a_v)
    end do
  end function esmoexp

  pure elemental function dperks(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((a_v * (a_v + 1.0_dp)) * exp((b_v * x))) / ((1.0_dp + (a_v * exp((b_v * x)))) &
      & ** 2.0_dp))
    if (((log_pdf_v))) res = ((log((a_v * (a_v + 1.0_dp))) + (b_v * x)) - (2.0_dp * log((1.0_dp + (a_v * exp((b_v * &
      & x)))))))
  end function dperks

  pure elemental function pperks(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((a_v * (exp((b_v * x)) - 1.0_dp)) / (1.0_dp + (a_v * &
      & exp((b_v * x)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp + a_v) / (1.0_dp + (a_v * exp((b_v * x)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = ((log(a_v) + log((exp((b_v * x)) - 1.0_dp))) - log((1.0_dp + &
      & (a_v * exp((b_v * x))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (log((1.0_dp + a_v)) - log((1.0_dp + (a_v * exp((b_v * &
      & x))))))
  end function pperks

  pure elemental function varperks(p, a, b, log_p, lower_tail) result(res)
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
    res = ((1.0_dp / b_v) * ((log((a_v + pp)) - log(a_v)) - log((1.0_dp - pp))))
  end function varperks

  pure elemental function esperks(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varperks(p * s_quad**3, a_v, b_v)
    end do
  end function esperks

  pure elemental function dbeard(x, a, b, rho, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, rho
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: rho_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    rho_v = 1.0_dp
    if (present(rho)) rho_v = rho
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((a_v * exp((b_v * x))) * ((1.0_dp + (a_v * rho_v)) ** (rho_v ** ((-1.0_dp) / &
      & b_v)))) / ((1.0_dp + ((a_v * rho_v) * exp((b_v * x)))) ** (1.0_dp + (rho_v ** ((-1.0_dp) / b_v)))))
    if (((log_pdf_v))) res = (((log(a_v) + (b_v * x)) + ((rho_v ** ((-1.0_dp) / b_v)) * log((1.0_dp + (a_v * &
      & rho_v))))) - ((1.0_dp + (rho_v ** ((-1.0_dp) / b_v))) * log((1.0_dp + ((a_v * rho_v) * exp((b_v * x)))))))
  end function dbeard

  pure elemental function pbeard(x, a, b, rho, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, rho
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: rho_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    rho_v = 1.0_dp
    if (present(rho)) rho_v = rho
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - (((1.0_dp + (a_v * rho_v)) / (1.0_dp + ((a_v * &
      & rho_v) * exp((b_v * x))))) ** (rho_v ** ((-1.0_dp) / b_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (((1.0_dp + (a_v * rho_v)) / (1.0_dp + ((a_v * &
      & rho_v) * exp((b_v * x))))) ** (rho_v ** ((-1.0_dp) / b_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - (((1.0_dp + (a_v * rho_v)) / (1.0_dp + ((a_v * &
      & rho_v) * exp((b_v * x))))) ** (rho_v ** ((-1.0_dp) / b_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((rho_v ** ((-1.0_dp) / b_v)) * (log((1.0_dp + (a_v * &
      & rho_v))) - log((1.0_dp + ((a_v * rho_v) * exp((b_v * x)))))))
  end function pbeard

  pure elemental function varbeard(p, a, b, rho, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, rho
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: rho_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    rho_v = 1.0_dp
    if (present(rho)) rho_v = rho
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / b_v) * (log((((1.0_dp + (a_v * rho_v)) / ((1.0_dp - pp) ** (rho_v ** (1.0_dp / b_v)))) - &
      & 1.0_dp)) - log((a_v * rho_v))))
  end function varbeard

  pure elemental function esbeard(p, a, b, rho) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, rho
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: rho_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    rho_v = 1.0_dp
    if (present(rho)) rho_v = rho
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbeard(p * s_quad**3, a_v, b_v, rho_v)
    end do
  end function esbeard

  pure elemental function dgompertz(x, b, eta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b, eta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: eta_v
    logical :: log_pdf_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((b_v * eta_v) * exp((b_v * x))) * exp((eta_v - (eta_v * exp((b_v * x))))))
    if (((log_pdf_v))) res = (((log((b_v * eta_v)) + (b_v * x)) + eta_v) - (eta_v * exp((b_v * x))))
  end function dgompertz

  pure elemental function pgompertz(x, b, eta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b, eta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: eta_v
    logical :: log_p_v
    logical :: lower_tail_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp((eta_v - (eta_v * exp((b_v * x))))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp((eta_v - (eta_v * exp((b_v * x)))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp((eta_v - (eta_v * exp((b_v * x)))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (eta_v - (eta_v * exp((b_v * x))))
  end function pgompertz

  pure elemental function vargompertz(p, b, eta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, eta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: eta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    b_v = 1.0_dp
    if (present(b)) b_v = b
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / b_v) * log((1.0_dp - ((1.0_dp / eta_v) * log((1.0_dp - pp))))))
  end function vargompertz

  pure elemental function esgompertz(p, b, eta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, eta
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: eta_v
    integer :: i
    real(dp) :: s_quad
    b_v = 1.0_dp
    if (present(b)) b_v = b
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargompertz(p * s_quad**3, b_v, eta_v)
    end do
  end function esgompertz

  pure elemental function dbetagompertz(x, b, c, d, eta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b, c, d, eta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    real(dp) :: eta_v
    logical :: log_pdf_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((b_v * eta_v) * exp((((d_v * eta_v) - ((d_v * eta_v) * exp((b_v * x)))) + (b_v &
      & * x)))) * ((1.0_dp - exp((eta_v - (eta_v * exp((b_v * x)))))) ** (c_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((((log((b_v * eta_v)) + (d_v * eta_v)) - ((d_v * eta_v) * exp((b_v * x)))) + (b_v * &
      & x)) + ((c_v - 1.0_dp) * log((1.0_dp - exp((eta_v - (eta_v * exp((b_v * x)))))))))
  end function dbetagompertz

  pure elemental function pbetagompertz(x, b, c, d, eta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: b, c, d, eta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    real(dp) :: eta_v
    logical :: log_p_v
    logical :: lower_tail_v
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf((1.0_dp - exp((eta_v - (eta_v * exp((b_v * x)))))), shape1=c_v, shape2=d_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetagompertz

  pure elemental function varbetagompertz(p, b, c, d, eta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, c, d, eta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    real(dp) :: eta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = ((1.0_dp / b_v) * log((1.0_dp - ((1.0_dp / eta_v) * log((1.0_dp - beta_quantile(pp, shape1=c_v, &
      & shape2=d_v)))))))
  end function varbetagompertz

  pure elemental function esbetagompertz(p, b, c, d, eta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: b, c, d, eta
    real(dp) :: res
    real(dp) :: b_v
    real(dp) :: c_v
    real(dp) :: d_v
    real(dp) :: eta_v
    integer :: i
    real(dp) :: s_quad
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    eta_v = 1.0_dp
    if (present(eta)) eta_v = eta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetagompertz(p * s_quad**3, b_v, c_v, d_v, eta_v)
    end do
  end function esbetagompertz

end module vares_distributions_01
