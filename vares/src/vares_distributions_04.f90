! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_04
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dgenbeta, pgenbeta, vargenbeta, esgenbeta, darcsine, parcsine
  public :: vararcsine, esarcsine, dtriangular, ptriangular, vartriangular, estriangular
  public :: dgenbeta2, pgenbeta2, vargenbeta2, esgenbeta2, dinvbeta, pinvbeta
  public :: varinvbeta, esinvbeta, dgeninvbeta, pgeninvbeta, vargeninvbeta, esgeninvbeta
  public :: dtsp, ptsp, vartsp, estsp, dkum, pkum
  public :: varkum, eskum, dnormal, pnormal, varnormal, esnormal
  public :: dkumnormal, pkumnormal, varkumnormal, eskumnormal, dexppower, pexppower
  public :: varexppower, esexppower, daep, paep, varaep, esaep
contains
  pure elemental function dgenbeta(x, a, b, c, d, log_pdf) result(res)
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
    c_v = 0.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((1.0_dp / (d_v - c_v)) * beta_pdf(((x - c_v) / (d_v - c_v)), shape1=a_v, &
      & shape2=b_v))
    if (((log_pdf_v))) res = (beta_pdf(((x - c_v) / (d_v - c_v)), shape1=a_v, shape2=b_v, log_pdf=.true.) - &
      & log((d_v - c_v)))
  end function dgenbeta

  pure elemental function pgenbeta(x, a, b, c, d, log_p, lower_tail) result(res)
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
    c_v = 0.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(((x - c_v) / (d_v - c_v)), shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pgenbeta

  pure elemental function vargenbeta(p, a, b, c, d, log_p, lower_tail) result(res)
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
    c_v = 0.0_dp
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
    res = (c_v + ((d_v - c_v) * beta_quantile(pp, shape1=a_v, shape2=b_v)))
  end function vargenbeta

  pure elemental function esgenbeta(p, a, b, c, d) result(res)
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
    c_v = 0.0_dp
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenbeta(p * s_quad**3, a_v, b_v, c_v, d_v)
    end do
  end function esgenbeta

  pure elemental function darcsine(x, a, b, log_pdf) result(res)
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
    res = dgenbeta(x, a=0.5_dp, b=0.5_dp, c=a_v, d=b_v, log_pdf=log_pdf_v)
  end function darcsine

  pure elemental function parcsine(x, a, b, log_p, lower_tail) result(res)
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
    res = pgenbeta(x, a=0.5_dp, b=0.5_dp, c=a_v, d=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function parcsine

  pure elemental function vararcsine(p, a, b, log_p, lower_tail) result(res)
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
    res = vargenbeta(pp, a=0.5_dp, b=0.5_dp, c=a_v, d=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function vararcsine

  pure elemental function esarcsine(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vararcsine(p * s_quad**3, a_v, b_v)
    end do
  end function esarcsine

  pure elemental function dtriangular(x, a, b, c, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_pdf_v
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 2.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x < a_v)))) res = 0.0_dp
    if ((((.not. log_pdf_v)) .and. ((x >= a_v)) .and. ((x <= c_v)))) res = ((2.0_dp * (x - a_v)) / ((b_v - a_v) * &
      & (c_v - a_v)))
    if ((((.not. log_pdf_v)) .and. ((x > c_v)) .and. ((x <= b_v)))) res = ((2.0_dp * (b_v - x)) / ((b_v - a_v) * &
      & (b_v - c_v)))
    if ((((.not. log_pdf_v)) .and. ((x > b_v)))) res = 1.0_dp
    if ((((log_pdf_v)) .and. ((x < a_v)))) res = (-huge(1.0_dp))
    if ((((log_pdf_v)) .and. ((x >= a_v)) .and. ((x <= c_v)))) res = log(((2.0_dp * (x - a_v)) / ((b_v - a_v) * &
      & (c_v - a_v))))
    if ((((log_pdf_v)) .and. ((x > c_v)) .and. ((x <= b_v)))) res = log(((2.0_dp * (b_v - x)) / ((b_v - a_v) * (b_v &
      & - c_v))))
    if ((((log_pdf_v)) .and. ((x > b_v)))) res = 0.0_dp
  end function dtriangular

  pure elemental function ptriangular(x, a, b, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 2.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x < a_v)))) res = 0.0_dp
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x >= a_v)) .and. ((x <= c_v)))) res = (((x - a_v) ** &
      & 2.0_dp) / ((b_v - a_v) * (c_v - a_v)))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > c_v)) .and. ((x <= b_v)))) res = (1.0_dp - (((b_v - &
      & x) ** 2.0_dp) / ((b_v - a_v) * (b_v - c_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > b_v)))) res = 1.0_dp
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x < a_v)))) res = 1.0_dp
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x >= a_v)) .and. ((x <= c_v)))) res = (1.0_dp - &
      & (((x - a_v) ** 2.0_dp) / ((b_v - a_v) * (c_v - a_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > c_v)) .and. ((x <= b_v)))) res = (((b_v - x) ** &
      & 2.0_dp) / ((b_v - a_v) * (b_v - c_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > b_v)))) res = 0.0_dp
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x < a_v)))) res = (-huge(1.0_dp))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x >= a_v)) .and. ((x <= c_v)))) res = log((((x - a_v) ** &
      & 2.0_dp) / ((b_v - a_v) * (c_v - a_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > c_v)) .and. ((x <= b_v)))) res = log((1.0_dp - (((b_v - x) &
      & ** 2.0_dp) / ((b_v - a_v) * (b_v - c_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > b_v)))) res = 0.0_dp
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x < a_v)))) res = 0.0_dp
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x >= a_v)) .and. ((x <= c_v)))) res = log((1.0_dp - (((x &
      & - a_v) ** 2.0_dp) / ((b_v - a_v) * (c_v - a_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > c_v)) .and. ((x <= b_v)))) res = log((((b_v - x) ** &
      & 2.0_dp) / ((b_v - a_v) * (b_v - c_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > b_v)))) res = (-huge(1.0_dp))
  end function ptriangular

  pure elemental function vartriangular(p, a, b, c, log_p, lower_tail) result(res)
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
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 2.0_dp
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
    res = pp
    if (((pp < ((c_v - a_v) / (b_v - a_v))))) res = (a_v + sqrt(((pp * (b_v - a_v)) * (c_v - a_v))))
    if (((pp >= ((c_v - a_v) / (b_v - a_v))))) res = (b_v - sqrt((((1.0_dp - pp) * (b_v - a_v)) * (b_v - c_v))))
  end function vartriangular

  pure elemental function estriangular(p, a, b, c) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, c
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.0_dp
    if (present(a)) a_v = a
    b_v = 2.0_dp
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vartriangular(p * s_quad**3, a_v, b_v, c_v)
    end do
  end function estriangular

  pure elemental function dgenbeta2(x, a, b, c, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((c_v * (x ** (c_v - 1.0_dp))) * beta_pdf((x ** c_v), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = ((beta_pdf((x ** c_v), shape1=a_v, shape2=b_v, log_pdf=.true.) + log(c_v)) + ((c_v - &
      & 1.0_dp) * log(x)))
  end function dgenbeta2

  pure elemental function pgenbeta2(x, a, b, c, log_p, lower_tail) result(res)
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
    res = beta_cdf((x ** c_v), shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pgenbeta2

  pure elemental function vargenbeta2(p, a, b, c, log_p, lower_tail) result(res)
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
    res = (beta_quantile(pp, shape1=a_v, shape2=b_v) ** (1.0_dp / c_v))
  end function vargenbeta2

  pure elemental function esgenbeta2(p, a, b, c) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenbeta2(p * s_quad**3, a_v, b_v, c_v)
    end do
  end function esgenbeta2

  pure elemental function dinvbeta(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((1.0_dp + x) ** (-2.0_dp)) * beta_pdf((x / (1.0_dp + x)), shape1=a_v, &
      & shape2=b_v))
    if (((log_pdf_v))) res = (beta_pdf((x / (1.0_dp + x)), shape1=a_v, shape2=b_v, log_pdf=.true.) - (2.0_dp * &
      & log((1.0_dp + x))))
  end function dinvbeta

  pure elemental function pinvbeta(x, a, b, log_p, lower_tail) result(res)
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
    res = beta_cdf((x / (1.0_dp + x)), shape1=a_v, shape2=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pinvbeta

  pure elemental function varinvbeta(p, a, b, log_p, lower_tail) result(res)
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
    res = (beta_quantile(pp, shape1=a_v, shape2=b_v) / (1.0_dp - beta_quantile(pp, shape1=a_v, shape2=b_v)))
  end function varinvbeta

  pure elemental function esinvbeta(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varinvbeta(p * s_quad**3, a_v, b_v)
    end do
  end function esinvbeta

  pure elemental function dgeninvbeta(x, a, c, d, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, c, d
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((a_v * (x ** (a_v - 1.0_dp))) * ((1.0_dp + (x ** a_v)) ** (-2.0_dp))) * &
      & beta_pdf(((x ** a_v) / (1.0_dp + (x ** a_v))), shape1=c_v, shape2=d_v))
    if (((log_pdf_v))) res = (((log(a_v) + ((a_v - 1.0_dp) * log(x))) - (2.0_dp * log((1.0_dp + (x ** a_v))))) + &
      & beta_pdf(((x ** a_v) / (1.0_dp + (x ** a_v))), shape1=c_v, shape2=d_v, log_pdf=.true.))
  end function dgeninvbeta

  pure elemental function pgeninvbeta(x, a, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(((x ** a_v) / (1.0_dp + (x ** a_v))), shape1=c_v, shape2=d_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pgeninvbeta

  pure elemental function vargeninvbeta(p, a, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
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
    res = ((beta_quantile(pp, shape1=c_v, shape2=d_v) / (1.0_dp - beta_quantile(pp, shape1=c_v, shape2=d_v))) ** &
      & (1.0_dp / a_v))
  end function vargeninvbeta

  pure elemental function esgeninvbeta(p, a, c, d) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, c, d
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargeninvbeta(p * s_quad**3, a_v, c_v, d_v)
    end do
  end function esgeninvbeta

  pure elemental function dtsp(x, a, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= theta_v)))) res = (a_v * ((x / theta_v) ** (a_v - 1.0_dp)))
    if ((((.not. log_pdf_v)) .and. ((x > theta_v)))) res = (a_v * (((1.0_dp - x) / (1.0_dp - theta_v)) ** (a_v - &
      & 1.0_dp)))
    if ((((log_pdf_v)) .and. ((x <= theta_v)))) res = (log(a_v) + ((a_v - 1.0_dp) * log((x / theta_v))))
    if ((((log_pdf_v)) .and. ((x > theta_v)))) res = (log(a_v) + ((a_v - 1.0_dp) * log(((1.0_dp - x) / (1.0_dp - &
      & theta_v)))))
  end function dtsp

  pure elemental function ptsp(x, a, theta, log_p, lower_tail) result(res)
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
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (theta_v * ((x / theta_v) ** a_v))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = (1.0_dp - ((1.0_dp - theta_v) * &
      & (((1.0_dp - x) / (1.0_dp - theta_v)) ** a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = (1.0_dp - (theta_v * ((x / &
      & theta_v) ** a_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = ((1.0_dp - theta_v) * &
      & (((1.0_dp - x) / (1.0_dp - theta_v)) ** a_v))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (((1.0_dp - a_v) * log(theta_v)) + (a_v &
      & * log(x)))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = log((1.0_dp - ((1.0_dp - theta_v) * &
      & (((1.0_dp - x) / (1.0_dp - theta_v)) ** a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = log((1.0_dp - (theta_v * ((x / &
      & theta_v) ** a_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = (((1.0_dp - a_v) * log((1.0_dp - &
      & theta_v))) + (a_v * log((1.0_dp - x))))
  end function ptsp

  pure elemental function vartsp(p, a, theta, log_p, lower_tail) result(res)
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
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= theta_v))) res = (theta_v * ((pp / theta_v) ** (1.0_dp / a_v)))
    if (((pp > theta_v))) res = (1.0_dp - ((1.0_dp - theta_v) * (((1.0_dp - pp) / (1.0_dp - theta_v)) ** (1.0_dp / &
      & a_v))))
  end function vartsp

  pure elemental function estsp(p, a, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 0.5_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vartsp(p * s_quad**3, a_v, theta_v)
    end do
  end function estsp

  pure elemental function dkum(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (((a_v * b_v) * (x ** (a_v - 1.0_dp))) * ((1.0_dp - (x ** a_v)) ** (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((log((a_v * b_v)) + ((a_v - 1.0_dp) * log(x))) + ((b_v - 1.0_dp) * log((1.0_dp - (x &
      & ** a_v)))))
  end function dkum

  pure elemental function pkum(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - (x ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - (x ** a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - (x ** a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - (x ** a_v))))
  end function pkum

  pure elemental function varkum(p, a, b, log_p, lower_tail) result(res)
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
    res = ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v))
  end function varkum

  pure elemental function eskum(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkum(p * s_quad**3, a_v, b_v)
    end do
  end function eskum

  pure elemental function dnormal(x, mu, sigma, log_pdf) result(res)
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
    res = normal_pdf(x, mean=mu_v, sd=sigma_v, log_pdf=log_pdf_v)
  end function dnormal

  pure elemental function pnormal(x, mu, sigma, log_p, lower_tail) result(res)
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
    res = normal_cdf(x, mean=mu_v, sd=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pnormal

  pure elemental function varnormal(p, mu, sigma, log_p, lower_tail) result(res)
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
    res = normal_quantile(pp, mean=mu_v, sd=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varnormal

  pure elemental function esnormal(p, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varnormal(p * s_quad**3, mu_v, sigma_v)
    end do
  end function esnormal

  pure elemental function dkumnormal(x, mu, sigma, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma, a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * b_v) * normal_pdf(x, mean=mu_v, sd=sigma_v)) * (normal_cdf(x, &
      & mean=mu_v, sd=sigma_v) ** (a_v - 1.0_dp))) * ((1.0_dp - (normal_cdf(x, mean=mu_v, sd=sigma_v) ** a_v)) ** &
      & (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = (((log((a_v * b_v)) + normal_pdf(x, mean=mu_v, sd=sigma_v, log_pdf=.true.)) + ((a_v - &
      & 1.0_dp) * normal_cdf(x, mean=mu_v, sd=sigma_v, log_p=.true.))) + ((b_v - 1.0_dp) * log((1.0_dp - &
      & (normal_cdf(x, mean=mu_v, sd=sigma_v) ** a_v)))))
  end function dkumnormal

  pure elemental function pkumnormal(x, mu, sigma, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma, a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - (normal_cdf(x, mean=mu_v, &
      & sd=sigma_v) ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - (normal_cdf(x, mean=mu_v, sd=sigma_v) ** &
      & a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - (normal_cdf(x, mean=mu_v, sd=sigma_v) &
      & ** a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - (normal_cdf(x, mean=mu_v, &
      & sd=sigma_v) ** a_v))))
  end function pkumnormal

  pure elemental function varkumnormal(p, mu, sigma, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma, a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
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
    res = (mu_v + (sigma_v * normal_quantile(((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v)))))
  end function varkumnormal

  pure elemental function eskumnormal(p, mu, sigma, a, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma, a, b
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumnormal(p * s_quad**3, mu_v, sigma_v, a_v, b_v)
    end do
  end function eskumnormal

  pure elemental function dexppower(x, mu, sigma, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    logical :: log_pdf_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (exp((((-(abs((x - mu_v)) ** a_v)) * (sigma_v ** (-a_v))) / a_v)) / (((2.0_dp * &
      & (a_v ** (1.0_dp / a_v))) * sigma_v) * gamma_fn((1.0_dp + (1.0_dp / a_v)))))
    if (((log_pdf_v))) res = (((((((-(abs((x - mu_v)) ** a_v)) * (sigma_v ** (-a_v))) / a_v) - log(2.0_dp)) - &
      & ((1.0_dp / a_v) * log(a_v))) - log(sigma_v)) - log_gamma_fn((1.0_dp + (1.0_dp / a_v))))
  end function dexppower

  pure elemental function pexppower(x, mu, sigma, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: mu, sigma, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= mu_v)))) res = (0.5_dp - (0.5_dp * gamma_cdf((((mu_v &
      & - x) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > mu_v)))) res = (0.5_dp + (0.5_dp * gamma_cdf((((x - &
      & mu_v) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= mu_v)))) res = (0.5_dp + (0.5_dp * &
      & gamma_cdf((((mu_v - x) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > mu_v)))) res = (0.5_dp - (0.5_dp * &
      & gamma_cdf((((x - mu_v) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= mu_v)))) res = log((0.5_dp - (0.5_dp * gamma_cdf((((mu_v - &
      & x) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > mu_v)))) res = log((0.5_dp + (0.5_dp * gamma_cdf((((x - &
      & mu_v) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= mu_v)))) res = log((0.5_dp + (0.5_dp * &
      & gamma_cdf((((mu_v - x) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > mu_v)))) res = log((0.5_dp - (0.5_dp * gamma_cdf((((x &
      & - mu_v) ** a_v) / (a_v * (sigma_v ** a_v))), shape=(1.0_dp / a_v)))))
  end function pexppower

  pure elemental function varexppower(p, mu, sigma, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= 0.5_dp))) res = (mu_v - (((a_v ** (1.0_dp / a_v)) * sigma_v) * (gamma_quantile((1.0_dp - (2.0_dp * &
      & pp)), shape=(1.0_dp / a_v)) ** (1.0_dp / a_v))))
    if (((pp > 0.5_dp))) res = (mu_v + (((a_v ** (1.0_dp / a_v)) * sigma_v) * (gamma_quantile(((2.0_dp * pp) - &
      & 1.0_dp), shape=(1.0_dp / a_v)) ** (1.0_dp / a_v))))
  end function varexppower

  pure elemental function esexppower(p, mu, sigma, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: mu, sigma, a
    real(dp) :: res
    real(dp) :: mu_v
    real(dp) :: sigma_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    mu_v = 0.0_dp
    if (present(mu)) mu_v = mu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varexppower(p * s_quad**3, mu_v, sigma_v, a_v)
    end do
  end function esexppower

  pure elemental function daep(x, q1, q2, alpha, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: q1, q2, alpha
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: q1_v
    real(dp) :: q2_v
    real(dp) :: alpha_v
    logical :: log_pdf_v
    real(dp) :: k1_local
    real(dp) :: k2_local
    real(dp) :: alphastar_local
    q1_v = 1.0_dp
    if (present(q1)) q1_v = q1
    q2_v = 1.0_dp
    if (present(q2)) q2_v = q2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    k1_local = (1.0_dp / ((2.0_dp * (q1_v ** (1.0_dp / q1_v))) * gamma_fn((1.0_dp + (1.0_dp / q1_v)))))
    k2_local = (1.0_dp / ((2.0_dp * (q2_v ** (1.0_dp / q2_v))) * gamma_fn((1.0_dp + (1.0_dp / q2_v)))))
    alphastar_local = ((alpha_v * k1_local) / ((alpha_v * k1_local) + ((1.0_dp - alpha_v) * k2_local)))
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= 0.0_dp)))) res = ((alpha_v * exp(((-(1.0_dp / q1_v)) * (abs((x / (2.0_dp &
      & * alphastar_local))) ** q1_v)))) / (((2.0_dp * alphastar_local) * (q1_v ** (1.0_dp / q1_v))) * &
      & gamma_fn((1.0_dp + (1.0_dp / q1_v)))))
    if ((((.not. log_pdf_v)) .and. ((x > 0.0_dp)))) res = (((1.0_dp - alpha_v) * exp(((-(1.0_dp / q2_v)) * (abs((x &
      & / (2.0_dp - (2.0_dp * alphastar_local)))) ** q2_v)))) / (((2.0_dp * (1.0_dp - alphastar_local)) * (q2_v ** &
      & (1.0_dp / q2_v))) * gamma_fn((1.0_dp + (1.0_dp / q2_v)))))
    if ((((log_pdf_v)) .and. ((x <= 0.0_dp)))) res = (((((log(alpha_v) - ((1.0_dp / q1_v) * (abs((x / (2.0_dp * &
      & alphastar_local))) ** q1_v))) - log(2.0_dp)) - log(alphastar_local)) - ((1.0_dp / q1_v) * log(q1_v))) - &
      & log_gamma_fn((1.0_dp + (1.0_dp / q1_v))))
    if ((((log_pdf_v)) .and. ((x > 0.0_dp)))) res = (((((log((1.0_dp - alpha_v)) - ((1.0_dp / q2_v) * (abs((x / &
      & (2.0_dp - (2.0_dp * alphastar_local)))) ** q2_v))) - log(2.0_dp)) - log((1.0_dp - alphastar_local))) - &
      & ((1.0_dp / q2_v) * log(q2_v))) - log_gamma_fn((1.0_dp + (1.0_dp / q2_v))))
  end function daep

  pure elemental function paep(x, q1, q2, alpha, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: q1, q2, alpha
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: q1_v
    real(dp) :: q2_v
    real(dp) :: alpha_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: k1_local
    real(dp) :: k2_local
    real(dp) :: alphastar_local
    q1_v = 1.0_dp
    if (present(q1)) q1_v = q1
    q2_v = 1.0_dp
    if (present(q2)) q2_v = q2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    k1_local = (1.0_dp / ((2.0_dp * (q1_v ** (1.0_dp / q1_v))) * gamma_fn((1.0_dp + (1.0_dp / q1_v)))))
    k2_local = (1.0_dp / ((2.0_dp * (q2_v ** (1.0_dp / q2_v))) * gamma_fn((1.0_dp + (1.0_dp / q2_v)))))
    alphastar_local = ((alpha_v * k1_local) / ((alpha_v * k1_local) + ((1.0_dp - alpha_v) * k2_local)))
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= 0.0_dp)))) res = (alpha_v - (alpha_v * &
      & gamma_cdf(((1.0_dp / q1_v) * (((-x) / (2.0_dp * alphastar_local)) ** q1_v)), shape=(1.0_dp / q1_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > 0.0_dp)))) res = (alpha_v + ((1.0_dp - alpha_v) * &
      & gamma_cdf(((1.0_dp / q2_v) * (((-x) / (2.0_dp - (2.0_dp * alphastar_local))) ** q2_v)), shape=(1.0_dp / &
      & q2_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= 0.0_dp)))) res = ((1.0_dp - alpha_v) + &
      & (alpha_v * gamma_cdf(((1.0_dp / q1_v) * (((-x) / (2.0_dp * alphastar_local)) ** q1_v)), shape=(1.0_dp / &
      & q1_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > 0.0_dp)))) res = ((1.0_dp - alpha_v) - ((1.0_dp &
      & - alpha_v) * gamma_cdf(((1.0_dp / q2_v) * (((-x) / (2.0_dp - (2.0_dp * alphastar_local))) ** q2_v)), &
      & shape=(1.0_dp / q2_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= 0.0_dp)))) res = log((alpha_v - (alpha_v * &
      & gamma_cdf(((1.0_dp / q1_v) * (((-x) / (2.0_dp * alphastar_local)) ** q1_v)), shape=(1.0_dp / q1_v)))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > 0.0_dp)))) res = log((alpha_v + ((1.0_dp - alpha_v) * &
      & gamma_cdf(((1.0_dp / q2_v) * (((-x) / (2.0_dp - (2.0_dp * alphastar_local))) ** q2_v)), shape=(1.0_dp / &
      & q2_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= 0.0_dp)))) res = log(((1.0_dp - alpha_v) + (alpha_v &
      & * gamma_cdf(((1.0_dp / q1_v) * (((-x) / (2.0_dp * alphastar_local)) ** q1_v)), shape=(1.0_dp / q1_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > 0.0_dp)))) res = log(((1.0_dp - alpha_v) - ((1.0_dp - &
      & alpha_v) * gamma_cdf(((1.0_dp / q2_v) * (((-x) / (2.0_dp - (2.0_dp * alphastar_local))) ** q2_v)), &
      & shape=(1.0_dp / q2_v)))))
  end function paep

  pure elemental function varaep(p, q1, q2, alpha, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: q1, q2, alpha
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: q1_v
    real(dp) :: q2_v
    real(dp) :: alpha_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    real(dp) :: k1_local
    real(dp) :: k2_local
    real(dp) :: alphastar_local
    q1_v = 1.0_dp
    if (present(q1)) q1_v = q1
    q2_v = 1.0_dp
    if (present(q2)) q2_v = q2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    k1_local = (1.0_dp / ((2.0_dp * (q1_v ** (1.0_dp / q1_v))) * gamma_fn((1.0_dp + (1.0_dp / q1_v)))))
    k2_local = (1.0_dp / ((2.0_dp * (q2_v ** (1.0_dp / q2_v))) * gamma_fn((1.0_dp + (1.0_dp / q2_v)))))
    alphastar_local = ((alpha_v * k1_local) / ((alpha_v * k1_local) + ((1.0_dp - alpha_v) * k2_local)))
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= alpha_v))) res = (((-2.0_dp) * alphastar_local) * ((q1_v * gamma_quantile((1.0_dp - (pp / &
      & alpha_v)), shape=(1.0_dp / q1_v))) ** (1.0_dp / q1_v)))
    if (((pp > alpha_v))) res = (((-2.0_dp) * (1.0_dp - alphastar_local)) * ((q2_v * gamma_quantile((1.0_dp - &
      & ((1.0_dp - pp) / (1.0_dp - alpha_v))), shape=(1.0_dp / q2_v))) ** (1.0_dp / q2_v)))
  end function varaep

  pure elemental function esaep(p, q1, q2, alpha) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: q1, q2, alpha
    real(dp) :: res
    real(dp) :: q1_v
    real(dp) :: q2_v
    real(dp) :: alpha_v
    integer :: i
    real(dp) :: s_quad
    q1_v = 1.0_dp
    if (present(q1)) q1_v = q1
    q2_v = 1.0_dp
    if (present(q2)) q2_v = q2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varaep(p * s_quad**3, q1_v, q2_v, alpha_v)
    end do
  end function esaep

end module vares_distributions_04
