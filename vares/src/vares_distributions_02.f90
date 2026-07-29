! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_02
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dlfr, plfr, varlfr, eslfr, dpareto, ppareto
  public :: varpareto, espareto, dkumpareto, pkumpareto, varkumpareto, eskumpareto
  public :: df, pf, varf, esf, dgenpareto, pgenpareto
  public :: vargenpareto, esgenpareto, dbetapareto, pbetapareto, varbetapareto, esbetapareto
  public :: dparetostable, pparetostable, varparetostable, esparetostable, dgamma, pgamma
  public :: vargamma, esgamma, dkumgamma, pkumgamma, varkumgamma, eskumgamma
  public :: dnakagami, pnakagami, varnakagami, esnakagami, drgamma, prgamma
  public :: varrgamma, esrgamma, dclg, pclg, varclg, esclg
contains
  pure elemental function dlfr(x, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((a_v + (b_v * x)) * exp((((-a_v) * x) - (((b_v * x) * x) / 2.0_dp))))
    if (((log_pdf_v))) res = ((log((a_v + (b_v * x))) - (a_v * x)) - (((b_v * x) * x) / 2.0_dp))
  end function dlfr

  pure elemental function plfr(x, a, b, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp((((-a_v) * x) - (((b_v * x) * x) / 2.0_dp))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp((((-a_v) * x) - (((b_v * x) * x) / 2.0_dp)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp((((-a_v) * x) - (((b_v * x) * x) / 2.0_dp)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (((-a_v) * x) - (((b_v * x) * x) / 2.0_dp))
  end function plfr

  pure elemental function varlfr(p, a, b, log_p, lower_tail) result(res)
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
    res = ((1.0_dp / b_v) * (sqrt(((a_v * a_v) - ((2.0_dp * b_v) * log((1.0_dp - pp))))) - a_v))
  end function varlfr

  pure elemental function eslfr(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlfr(p * s_quad**3, a_v, b_v)
    end do
  end function eslfr

  pure elemental function dpareto(x, k, c, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((c_v * (k_v ** c_v)) * (x ** ((-c_v) - 1.0_dp)))
    if (((log_pdf_v))) res = ((log(c_v) + (c_v * log(k_v))) - ((c_v + 1.0_dp) * log(x)))
  end function dpareto

  pure elemental function ppareto(x, k, c, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((k_v ** c_v) * (x ** (-c_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((k_v ** c_v) * (x ** (-c_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((k_v ** c_v) * (x ** (-c_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (c_v * (log(k_v) - log(x)))
  end function ppareto

  pure elemental function varpareto(p, k, c, log_p, lower_tail) result(res)
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
    res = (k_v * ((1.0_dp - pp) ** ((-1.0_dp) / c_v)))
  end function varpareto

  pure elemental function espareto(p, k, c) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varpareto(p * s_quad**3, k_v, c_v)
    end do
  end function espareto

  pure elemental function dkumpareto(x, k, a, b, c, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k, a, b, c
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_pdf_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    c_v = 1.0_dp
    if (present(c)) c_v = c
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((((a_v * b_v) * c_v) * (k_v ** c_v)) * (x ** ((-c_v) - 1.0_dp))) * ((1.0_dp - &
      & ((k_v ** c_v) * (x ** (-c_v)))) ** (a_v - 1.0_dp))) * ((1.0_dp - ((1.0_dp - ((k_v ** c_v) * (x ** (-c_v)))) &
      & ** a_v)) ** (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = ((((log(((a_v * b_v) * c_v)) + (c_v * log(k_v))) - ((c_v + 1.0_dp) * log(x))) + ((a_v &
      & - 1.0_dp) * log((1.0_dp - ((k_v ** c_v) * (x ** (-c_v))))))) + ((b_v - 1.0_dp) * log((1.0_dp - ((1.0_dp - &
      & ((k_v ** c_v) * (x ** (-c_v)))) ** a_v)))))
  end function dkumpareto

  pure elemental function pkumpareto(x, k, a, b, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k, a, b, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - ((1.0_dp - ((k_v ** c_v) * (x ** &
      & (-c_v)))) ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - ((1.0_dp - ((k_v ** c_v) * (x ** &
      & (-c_v)))) ** a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - ((1.0_dp - ((k_v ** c_v) * (x ** &
      & (-c_v)))) ** a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - ((1.0_dp - ((k_v ** c_v) * (x ** &
      & (-c_v)))) ** a_v))))
  end function pkumpareto

  pure elemental function varkumpareto(p, k, a, b, c, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k, a, b, c
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    k_v = 1.0_dp
    if (present(k)) k_v = k
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
    res = (k_v * ((1.0_dp - ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / a_v))) ** ((-1.0_dp) / c_v)))
  end function varkumpareto

  pure elemental function eskumpareto(p, k, a, b, c) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k, a, b, c
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: c_v
    integer :: i
    real(dp) :: s_quad
    k_v = 1.0_dp
    if (present(k)) k_v = k
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumpareto(p * s_quad**3, k_v, a_v, b_v, c_v)
    end do
  end function eskumpareto

  pure elemental function df(x, d1, d2, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: d1, d2
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: d1_v
    real(dp) :: d2_v
    logical :: log_pdf_v
    d1_v = 1.0_dp
    if (present(d1)) d1_v = d1
    d2_v = 1.0_dp
    if (present(d2)) d2_v = d2
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = f_pdf(x, df1=d1_v, df2=d2_v, log_pdf=log_pdf_v)
  end function df

  pure elemental function pf(x, d1, d2, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: d1, d2
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: d1_v
    real(dp) :: d2_v
    logical :: log_p_v
    logical :: lower_tail_v
    d1_v = 1.0_dp
    if (present(d1)) d1_v = d1
    d2_v = 1.0_dp
    if (present(d2)) d2_v = d2
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = f_cdf(x, df1=d1_v, df2=d2_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pf

  pure elemental function varf(p, d1, d2, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: d1, d2
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: d1_v
    real(dp) :: d2_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    d1_v = 1.0_dp
    if (present(d1)) d1_v = d1
    d2_v = 1.0_dp
    if (present(d2)) d2_v = d2
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = f_quantile(pp, df1=d1_v, df2=d2_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varf

  pure elemental function esf(p, d1, d2) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: d1, d2
    real(dp) :: res
    real(dp) :: d1_v
    real(dp) :: d2_v
    integer :: i
    real(dp) :: s_quad
    d1_v = 1.0_dp
    if (present(d1)) d1_v = d1
    d2_v = 1.0_dp
    if (present(d2)) d2_v = d2
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varf(p * s_quad**3, d1_v, d2_v)
    end do
  end function esf

  pure elemental function dgenpareto(x, k, c, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((1.0_dp / k_v) * ((1.0_dp - ((c_v * x) / k_v)) ** ((1.0_dp / c_v) - 1.0_dp)))
    if (((log_pdf_v))) res = ((-log(k_v)) + (((1.0_dp / c_v) - 1.0_dp) * log((1.0_dp - ((c_v * x) / k_v)))))
  end function dgenpareto

  pure elemental function pgenpareto(x, k, c, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - ((c_v * x) / k_v)) ** (1.0_dp / c_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - ((c_v * x) / k_v)) ** (1.0_dp / c_v))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - ((c_v * x) / k_v)) ** (1.0_dp / c_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp / c_v) * log((1.0_dp - ((c_v * x) / k_v))))
  end function pgenpareto

  pure elemental function vargenpareto(p, k, c, log_p, lower_tail) result(res)
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
    res = ((k_v / c_v) * (1.0_dp - ((1.0_dp - pp) ** c_v)))
  end function vargenpareto

  pure elemental function esgenpareto(p, k, c) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargenpareto(p * s_quad**3, k_v, c_v)
    end do
  end function esgenpareto

  pure elemental function dbetapareto(x, k, a, c, d, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k, a, c, d
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_pdf_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
    a_v = 1.0_dp
    if (present(a)) a_v = a
    c_v = 1.0_dp
    if (present(c)) c_v = c
    d_v = 1.0_dp
    if (present(d)) d_v = d
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((a_v * (k_v ** (a_v * d_v))) * (x ** (((-a_v) * d_v) - 1.0_dp))) * ((1.0_dp - &
      & ((k_v / x) ** a_v)) ** (c_v - 1.0_dp))) / beta_fn(c_v, d_v))
    if (((log_pdf_v))) res = ((((log(a_v) + ((a_v * d_v) * log(k_v))) - (((a_v * d_v) + 1.0_dp) * log(x))) + ((c_v &
      & - 1.0_dp) * log((1.0_dp - ((k_v / x) ** a_v))))) - log_beta_fn(c_v, d_v))
  end function dbetapareto

  pure elemental function pbetapareto(x, k, a, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: k, a, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    k_v = 1.0_dp
    if (present(k)) k_v = k
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
    res = beta_cdf((1.0_dp - ((k_v / x) ** a_v)), shape1=c_v, shape2=d_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pbetapareto

  pure elemental function varbetapareto(p, k, a, c, d, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k, a, c, d
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    k_v = 1.0_dp
    if (present(k)) k_v = k
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
    res = (k_v * ((1.0_dp - beta_quantile(pp, shape1=c_v, shape2=d_v)) ** ((-1.0_dp) / a_v)))
  end function varbetapareto

  pure elemental function esbetapareto(p, k, a, c, d) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: k, a, c, d
    real(dp) :: res
    real(dp) :: k_v
    real(dp) :: a_v
    real(dp) :: c_v
    real(dp) :: d_v
    integer :: i
    real(dp) :: s_quad
    k_v = 1.0_dp
    if (present(k)) k_v = k
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetapareto(p * s_quad**3, k_v, a_v, c_v, d_v)
    end do
  end function esbetapareto

  pure elemental function dparetostable(x, lambda, nu, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, nu, sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: nu_v
    real(dp) :: sigma_v
    logical :: log_pdf_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    nu_v = 1.0_dp
    if (present(nu)) nu_v = nu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((nu_v * lambda_v) * (1.0_dp / x)) * (log((x / sigma_v)) ** (nu_v - 1.0_dp))) &
      & * exp(((-lambda_v) * (log((x / sigma_v)) ** nu_v))))
    if (((log_pdf_v))) res = (((log((nu_v * lambda_v)) - log(x)) + ((nu_v - 1.0_dp) * log(log((x / sigma_v))))) - &
      & (lambda_v * (log((x / sigma_v)) ** nu_v)))
  end function dparetostable

  pure elemental function pparetostable(x, lambda, nu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: lambda, nu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: nu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    nu_v = 1.0_dp
    if (present(nu)) nu_v = nu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - exp(((-lambda_v) * (log((x / sigma_v)) ** nu_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = exp(((-lambda_v) * (log((x / sigma_v)) ** nu_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - exp(((-lambda_v) * (log((x / sigma_v)) ** nu_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((-lambda_v) * (log((x / sigma_v)) ** nu_v))
  end function pparetostable

  pure elemental function varparetostable(p, lambda, nu, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, nu, sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: nu_v
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    nu_v = 1.0_dp
    if (present(nu)) nu_v = nu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (sigma_v * exp(((((-1.0_dp) / lambda_v) * log((1.0_dp - pp))) ** (1.0_dp / nu_v))))
  end function varparetostable

  pure elemental function esparetostable(p, lambda, nu, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: lambda, nu, sigma
    real(dp) :: res
    real(dp) :: lambda_v
    real(dp) :: nu_v
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    lambda_v = 1.0_dp
    if (present(lambda)) lambda_v = lambda
    nu_v = 1.0_dp
    if (present(nu)) nu_v = nu
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varparetostable(p * s_quad**3, lambda_v, nu_v, sigma_v)
    end do
  end function esparetostable

  pure elemental function dgamma(x, a, b, log_pdf) result(res)
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
    res = gamma_pdf(x, shape=a_v, rate=b_v, log_pdf=log_pdf_v)
  end function dgamma

  pure elemental function pgamma(x, a, b, log_p, lower_tail) result(res)
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
    res = gamma_cdf(x, shape=a_v, rate=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pgamma

  pure elemental function vargamma(p, a, b, log_p, lower_tail) result(res)
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
    res = gamma_quantile(pp, shape=a_v, rate=b_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function vargamma

  pure elemental function esgamma(p, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vargamma(p * s_quad**3, a_v, b_v)
    end do
  end function esgamma

  pure elemental function dkumgamma(x, a, b, c, d, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((((c_v * d_v) * gamma_pdf(x, shape=a_v, rate=b_v)) * (gamma_cdf(x, shape=a_v, &
      & rate=b_v) ** (c_v - 1.0_dp))) * ((1.0_dp - (gamma_cdf(x, shape=a_v, rate=b_v) ** c_v)) ** (d_v - 1.0_dp)))
    if (((log_pdf_v))) res = (((log((c_v * d_v)) + gamma_pdf(x, shape=a_v, rate=b_v, log_pdf=.true.)) + ((c_v - &
      & 1.0_dp) * gamma_cdf(x, shape=a_v, rate=b_v, log_p=.true.))) + ((d_v - 1.0_dp) * log((1.0_dp - (gamma_cdf(x, &
      & shape=a_v, rate=b_v) ** c_v)))))
  end function dkumgamma

  pure elemental function pkumgamma(x, a, b, c, d, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - (gamma_cdf(x, shape=a_v, rate=b_v) &
      & ** c_v)) ** d_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - (gamma_cdf(x, shape=a_v, rate=b_v) ** &
      & c_v)) ** d_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - (gamma_cdf(x, shape=a_v, rate=b_v) ** &
      & c_v)) ** d_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (d_v * log((1.0_dp - (gamma_cdf(x, shape=a_v, rate=b_v) &
      & ** c_v))))
  end function pkumgamma

  pure elemental function varkumgamma(p, a, b, c, d, log_p, lower_tail) result(res)
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
    res = ((1.0_dp / b_v) * gamma_quantile(((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / d_v))) ** (1.0_dp / c_v)), &
      & shape=a_v))
  end function varkumgamma

  pure elemental function eskumgamma(p, a, b, c, d) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumgamma(p * s_quad**3, a_v, b_v, c_v, d_v)
    end do
  end function eskumgamma

  pure elemental function dnakagami(x, m, a, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: m, a
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: m_v
    real(dp) :: a_v
    logical :: log_pdf_v
    m_v = 1.0_dp
    if (present(m)) m_v = m
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((((2.0_dp * ((m_v / a_v) ** m_v)) * (x ** ((2.0_dp * m_v) - 1.0_dp))) * &
      & exp(((((-m_v) * x) * x) / a_v))) / gamma_fn(m_v))
    if (((log_pdf_v))) res = ((((log(2.0_dp) + (m_v * log((m_v / a_v)))) + (((2.0_dp * m_v) - 1.0_dp) * log(x))) - &
      & (((m_v * x) * x) / a_v)) - log_gamma_fn(m_v))
  end function dnakagami

  pure elemental function pnakagami(x, m, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: m, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: m_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    m_v = 1.0_dp
    if (present(m)) m_v = m
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = gamma_cdf((((m_v * x) * x) / a_v), shape=m_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pnakagami

  pure elemental function varnakagami(p, m, a, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: m, a
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: m_v
    real(dp) :: a_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    m_v = 1.0_dp
    if (present(m)) m_v = m
    a_v = 1.0_dp
    if (present(a)) a_v = a
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (sqrt((a_v / m_v)) * sqrt(gamma_quantile(pp, shape=m_v)))
  end function varnakagami

  pure elemental function esnakagami(p, m, a) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: m, a
    real(dp) :: res
    real(dp) :: m_v
    real(dp) :: a_v
    integer :: i
    real(dp) :: s_quad
    m_v = 1.0_dp
    if (present(m)) m_v = m
    a_v = 1.0_dp
    if (present(a)) a_v = a
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varnakagami(p * s_quad**3, m_v, a_v)
    end do
  end function esnakagami

  pure elemental function drgamma(x, a, theta, phi, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta, phi
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = ((1.0_dp / (2.0_dp * phi_v)) * gamma_pdf((abs((x - theta_v)) / phi_v), shape=a_v))
    if (((log_pdf_v))) res = (gamma_pdf((abs((x - theta_v)) / phi_v), shape=a_v, log_pdf=.true.) - log((2.0_dp * &
      & phi_v)))
  end function drgamma

  pure elemental function prgamma(x, a, theta, phi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta, phi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((x <= theta_v)) .and. ((.not. log_p_v)) .and. ((lower_tail_v)))) res = (0.5_dp - (0.5_dp * &
      & gamma_cdf(((theta_v - x) / phi_v), shape=a_v)))
    if ((((x > theta_v)) .and. ((.not. log_p_v)) .and. ((lower_tail_v)))) res = (0.5_dp + (0.5_dp * gamma_cdf(((x - &
      & theta_v) / phi_v), shape=a_v)))
    if ((((x <= theta_v)) .and. ((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (0.5_dp + (0.5_dp * &
      & gamma_cdf(((theta_v - x) / phi_v), shape=a_v)))
    if ((((x > theta_v)) .and. ((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (0.5_dp - (0.5_dp * &
      & gamma_cdf(((x - theta_v) / phi_v), shape=a_v)))
    if ((((x <= theta_v)) .and. ((log_p_v)) .and. ((lower_tail_v)))) res = log((0.5_dp - (0.5_dp * &
      & gamma_cdf(((theta_v - x) / phi_v), shape=a_v))))
    if ((((x > theta_v)) .and. ((log_p_v)) .and. ((lower_tail_v)))) res = log((0.5_dp + (0.5_dp * gamma_cdf(((x - &
      & theta_v) / phi_v), shape=a_v))))
    if ((((x <= theta_v)) .and. ((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((0.5_dp + (0.5_dp * &
      & gamma_cdf(((theta_v - x) / phi_v), shape=a_v))))
    if ((((x > theta_v)) .and. ((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((0.5_dp - (0.5_dp * &
      & gamma_cdf(((x - theta_v) / phi_v), shape=a_v))))
  end function prgamma

  pure elemental function varrgamma(p, a, theta, phi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta, phi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    if (((pp <= 0.5_dp))) res = (theta_v - (phi_v * gamma_quantile((1.0_dp - (2.0_dp * pp)), shape=a_v)))
    if (((pp > 0.5_dp))) res = (theta_v + (phi_v * gamma_quantile(((2.0_dp * pp) - 1.0_dp), shape=a_v)))
  end function varrgamma

  pure elemental function esrgamma(p, a, theta, phi) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta, phi
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varrgamma(p * s_quad**3, a_v, theta_v, phi_v)
    end do
  end function esrgamma

  pure elemental function dclg(x, a, b, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((0.5_dp * a_v) * b_v) * ((1.0_dp + (b_v * abs((x - theta_v)))) ** ((-a_v) - &
      & 1.0_dp)))
    if (((log_pdf_v))) res = ((log((a_v * b_v)) - log(2.0_dp)) - ((a_v + 1.0_dp) * log((1.0_dp + (b_v * abs((x - &
      & theta_v)))))))
  end function dclg

  pure elemental function pclg(x, a, b, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, b, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((x <= theta_v)) .and. ((.not. log_p_v)) .and. ((lower_tail_v)))) res = (0.5_dp * ((1.0_dp + (b_v * abs((x &
      & - theta_v)))) ** (-a_v)))
    if ((((x > theta_v)) .and. ((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - (0.5_dp * ((1.0_dp + (b_v &
      & * abs((x - theta_v)))) ** (-a_v))))
    if ((((x <= theta_v)) .and. ((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - (0.5_dp * ((1.0_dp &
      & + (b_v * abs((x - theta_v)))) ** (-a_v))))
    if ((((x > theta_v)) .and. ((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (0.5_dp * ((1.0_dp + (b_v * &
      & abs((x - theta_v)))) ** (-a_v)))
    if ((((x <= theta_v)) .and. ((log_p_v)) .and. ((lower_tail_v)))) res = ((-log(2.0_dp)) - (a_v * log((1.0_dp + &
      & (b_v * abs((x - theta_v)))))))
    if ((((x > theta_v)) .and. ((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - (0.5_dp * ((1.0_dp + (b_v * &
      & abs((x - theta_v)))) ** (-a_v)))))
    if ((((x <= theta_v)) .and. ((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - (0.5_dp * ((1.0_dp + &
      & (b_v * abs((x - theta_v)))) ** (-a_v)))))
    if ((((x > theta_v)) .and. ((log_p_v)) .and. ((.not. lower_tail_v)))) res = ((-log(2.0_dp)) - (a_v * &
      & log((1.0_dp + (b_v * abs((x - theta_v)))))))
  end function pclg

  pure elemental function varclg(p, a, b, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
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
    if (((pp <= 0.5_dp))) res = ((theta_v + (1.0_dp / b_v)) - ((1.0_dp / b_v) * ((2.0_dp * pp) ** ((-1.0_dp) / a_v))))
    if (((pp > 0.5_dp))) res = ((theta_v - (1.0_dp / b_v)) + ((1.0_dp / b_v) * ((2.0_dp * (1.0_dp - pp)) ** &
      & ((-1.0_dp) / a_v))))
  end function varclg

  pure elemental function esclg(p, a, b, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, b, theta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: b_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varclg(p * s_quad**3, a_v, b_v, theta_v)
    end do
  end function esclg

end module vares_distributions_02
