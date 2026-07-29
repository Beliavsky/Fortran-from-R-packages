! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of VaRES 1.0.2.
! See NOTICE.md and original/ for attribution and provenance.
module vares_distributions_05
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use vares_kinds, only : dp, pi
  use vares_special
  use vares_quadrature, only : gl_n, gl_x, gl_w
  implicit none
  private
  public :: dbetanorm, pbetanorm, varbetanorm, esbetanorm, dhalfnorm, phalfnorm
  public :: varhalfnorm, eshalfnorm, dkumhalfnorm, pkumhalfnorm, varkumhalfnorm, eskumhalfnorm
  public :: dt, pt, vart, est, dast, past
  public :: varast, esast, dhalft, phalft, varhalft, eshalft
  public :: dcauchy, pcauchy, varcauchy, escauchy, dlogcauchy, plogcauchy
  public :: varlogcauchy, eslogcauchy, dhalfcauchy, phalfcauchy, varhalfcauchy, eshalfcauchy
  public :: dlaplace, plaplace, varlaplace, eslaplace, dpctalaplace, ppctalaplace
  public :: varpctalaplace, espctalaplace, dhblaplace, phblaplace, varhblaplace, eshblaplace
contains
  pure elemental function dbetanorm(x, mu, sigma, a, b, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = (normal_pdf(x, mean=mu_v, sd=sigma_v) * beta_pdf(normal_cdf(x, mean=mu_v, &
      & sd=sigma_v), shape1=a_v, shape2=b_v))
    if (((log_pdf_v))) res = (normal_pdf(x, mean=mu_v, sd=sigma_v, log_pdf=.true.) + beta_pdf(normal_cdf(x, &
      & mean=mu_v, sd=sigma_v), shape1=a_v, shape2=b_v, log_pdf=.true.))
  end function dbetanorm

  pure elemental function pbetanorm(x, mu, sigma, a, b, log_p, lower_tail) result(res)
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
    res = beta_cdf(normal_cdf(x, mean=mu_v, sd=sigma_v), shape1=a_v, shape2=b_v, log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function pbetanorm

  pure elemental function varbetanorm(p, mu, sigma, a, b, log_p, lower_tail) result(res)
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
    res = (mu_v + (sigma_v * normal_quantile(beta_quantile(pp, shape1=a_v, shape2=b_v))))
  end function varbetanorm

  pure elemental function esbetanorm(p, mu, sigma, a, b) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varbetanorm(p * s_quad**3, mu_v, sigma_v, a_v, b_v)
    end do
  end function esbetanorm

  pure elemental function dhalfnorm(x, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: sigma_v
    logical :: log_pdf_v
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (2.0_dp * normal_pdf(x, sd=sigma_v))
    if (((log_pdf_v))) res = (log(2.0_dp) + normal_pdf(x, sd=sigma_v, log_pdf=.true.))
  end function dhalfnorm

  pure elemental function phalfnorm(x, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((2.0_dp * normal_cdf(x, sd=sigma_v)) - 1.0_dp)
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (2.0_dp - (2.0_dp * normal_cdf(x, sd=sigma_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log(((2.0_dp * normal_cdf(x, sd=sigma_v)) - 1.0_dp))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((2.0_dp - (2.0_dp * normal_cdf(x, sd=sigma_v))))
  end function phalfnorm

  pure elemental function varhalfnorm(p, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (sigma_v * normal_quantile((0.5_dp * (1.0_dp + pp))))
  end function varhalfnorm

  pure elemental function eshalfnorm(p, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: sigma
    real(dp) :: res
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhalfnorm(p * s_quad**3, sigma_v)
    end do
  end function eshalfnorm

  pure elemental function dkumhalfnorm(x, sigma, a, b, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: sigma, a, b
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_pdf_v
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    a_v = 1.0_dp
    if (present(a)) a_v = a
    b_v = 1.0_dp
    if (present(b)) b_v = b
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (((((2.0_dp * a_v) * b_v) * normal_pdf(x, sd=sigma_v)) * (((2.0_dp * &
      & normal_cdf(x, sd=sigma_v)) - 1.0_dp) ** (a_v - 1.0_dp))) * ((1.0_dp - (((2.0_dp * normal_cdf(x, &
      & sd=sigma_v)) - 1.0_dp) ** a_v)) ** (b_v - 1.0_dp)))
    if (((log_pdf_v))) res = (((log(((2.0_dp * a_v) * b_v)) + normal_pdf(x, sd=sigma_v, log_pdf=.true.)) + ((a_v - &
      & 1.0_dp) * log(((2.0_dp * normal_cdf(x, sd=sigma_v)) - 1.0_dp)))) + ((b_v - 1.0_dp) * log((1.0_dp - &
      & (((2.0_dp * normal_cdf(x, sd=sigma_v)) - 1.0_dp) ** a_v)))))
  end function dkumhalfnorm

  pure elemental function pkumhalfnorm(x, sigma, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: sigma, a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (1.0_dp - ((1.0_dp - (((2.0_dp * normal_cdf(x, &
      & sd=sigma_v)) - 1.0_dp) ** a_v)) ** b_v))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((1.0_dp - (((2.0_dp * normal_cdf(x, sd=sigma_v)) - &
      & 1.0_dp) ** a_v)) ** b_v)
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((1.0_dp - ((1.0_dp - (((2.0_dp * normal_cdf(x, sd=sigma_v)) &
      & - 1.0_dp) ** a_v)) ** b_v)))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = (b_v * log((1.0_dp - (((2.0_dp * normal_cdf(x, &
      & sd=sigma_v)) - 1.0_dp) ** a_v))))
  end function pkumhalfnorm

  pure elemental function varkumhalfnorm(p, sigma, a, b, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: sigma, a, b
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
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
    res = (sigma_v * normal_quantile((0.5_dp + (0.5_dp * ((1.0_dp - ((1.0_dp - pp) ** (1.0_dp / b_v))) ** (1.0_dp / &
      & a_v))))))
  end function varkumhalfnorm

  pure elemental function eskumhalfnorm(p, sigma, a, b) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: sigma, a, b
    real(dp) :: res
    real(dp) :: sigma_v
    real(dp) :: a_v
    real(dp) :: b_v
    integer :: i
    real(dp) :: s_quad
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varkumhalfnorm(p * s_quad**3, sigma_v, a_v, b_v)
    end do
  end function eskumhalfnorm

  pure elemental function dt(x, n, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: n
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: n_v
    logical :: log_pdf_v
    n_v = 1.0_dp
    if (present(n)) n_v = n
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = student_t_pdf(x, df=n_v, log_pdf=log_pdf_v)
  end function dt

  pure elemental function pt(x, n, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: n
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: n_v
    logical :: log_p_v
    logical :: lower_tail_v
    n_v = 1.0_dp
    if (present(n)) n_v = n
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = student_t_cdf(x, df=n_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pt

  pure elemental function vart(p, n, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: n
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: n_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    n_v = 1.0_dp
    if (present(n)) n_v = n
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    res = student_t_quantile(pp, df=n_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function vart

  pure elemental function est(p, n) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: n
    real(dp) :: res
    real(dp) :: n_v
    integer :: i
    real(dp) :: s_quad
    n_v = 1.0_dp
    if (present(n)) n_v = n
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * vart(p * s_quad**3, n_v)
    end do
  end function est

  pure elemental function dast(x, nu1, nu2, alpha, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: nu1, nu2, alpha
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: nu1_v
    real(dp) :: nu2_v
    real(dp) :: alpha_v
    logical :: log_pdf_v
    real(dp) :: k1_local
    real(dp) :: k2_local
    real(dp) :: alphastar_local
    nu1_v = 1.0_dp
    if (present(nu1)) nu1_v = nu1
    nu2_v = 1.0_dp
    if (present(nu2)) nu2_v = nu2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    k1_local = (gamma_fn(((nu1_v + 1.0_dp) / 2.0_dp)) / (sqrt((nu1_v / 2.0_dp)) * gamma_fn((nu1_v / 2.0_dp))))
    k2_local = (gamma_fn(((nu2_v + 1.0_dp) / 2.0_dp)) / (sqrt((nu2_v / 2.0_dp)) * gamma_fn((nu2_v / 2.0_dp))))
    alphastar_local = ((alpha_v * k1_local) / ((alpha_v * k1_local) + ((1.0_dp - alpha_v) * k2_local)))
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= 0.0_dp)))) res = ((alpha_v / alphastar_local) * student_t_pdf((x / &
      & (2.0_dp * alphastar_local)), df=nu1_v))
    if ((((.not. log_pdf_v)) .and. ((x > 0.0_dp)))) res = (((1.0_dp - alpha_v) / (1.0_dp - alphastar_local)) * &
      & student_t_pdf((x / (2.0_dp * (1.0_dp - alphastar_local))), df=nu2_v))
    if ((((log_pdf_v)) .and. ((x <= 0.0_dp)))) res = (log((alpha_v / alphastar_local)) + student_t_pdf((x / (2.0_dp &
      & * alphastar_local)), df=nu1_v, log_pdf=.true.))
    if ((((log_pdf_v)) .and. ((x > 0.0_dp)))) res = (log(((1.0_dp - alpha_v) / (1.0_dp - alphastar_local))) + &
      & student_t_pdf((x / (2.0_dp * (1.0_dp - alphastar_local))), df=nu2_v, log_pdf=.true.))
  end function dast

  pure elemental function past(x, nu1, nu2, alpha, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: nu1, nu2, alpha
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: nu1_v
    real(dp) :: nu2_v
    real(dp) :: alpha_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: k1_local
    real(dp) :: k2_local
    real(dp) :: alphastar_local
    real(dp) :: minp_local
    real(dp) :: maxp_local
    nu1_v = 1.0_dp
    if (present(nu1)) nu1_v = nu1
    nu2_v = 1.0_dp
    if (present(nu2)) nu2_v = nu2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    k1_local = (gamma_fn(((nu1_v + 1.0_dp) / 2.0_dp)) / (sqrt((nu1_v / 2.0_dp)) * gamma_fn((nu1_v / 2.0_dp))))
    k2_local = (gamma_fn(((nu2_v + 1.0_dp) / 2.0_dp)) / (sqrt((nu2_v / 2.0_dp)) * gamma_fn((nu2_v / 2.0_dp))))
    alphastar_local = ((alpha_v * k1_local) / ((alpha_v * k1_local) + ((1.0_dp - alpha_v) * k2_local)))
    res = x
    minp_local = x
    maxp_local = x
    minp_local = min(x, 0.0_dp)
    maxp_local = max(x, 0.0_dp)
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = (((((2.0_dp * alpha_v) * student_t_cdf((minp_local / &
      & (2.0_dp * alphastar_local)), df=nu1_v)) - 1.0_dp) + alpha_v) + ((2.0_dp * (1.0_dp - alpha_v)) * &
      & student_t_cdf((maxp_local / (2.0_dp - (2.0_dp * alphastar_local))), df=nu2_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = ((((1.0_dp - ((2.0_dp * alpha_v) * &
      & student_t_cdf((minp_local / (2.0_dp * alphastar_local)), df=nu1_v))) + 1.0_dp) - alpha_v) - ((2.0_dp * &
      & (1.0_dp - alpha_v)) * student_t_cdf((maxp_local / (2.0_dp - (2.0_dp * alphastar_local))), df=nu2_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log((((((2.0_dp * alpha_v) * student_t_cdf((minp_local / &
      & (2.0_dp * alphastar_local)), df=nu1_v)) - 1.0_dp) + alpha_v) + ((2.0_dp * (1.0_dp - alpha_v)) * &
      & student_t_cdf((maxp_local / (2.0_dp - (2.0_dp * alphastar_local))), df=nu2_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log(((((1.0_dp - ((2.0_dp * alpha_v) * &
      & student_t_cdf((minp_local / (2.0_dp * alphastar_local)), df=nu1_v))) + 1.0_dp) - alpha_v) - ((2.0_dp * &
      & (1.0_dp - alpha_v)) * student_t_cdf((maxp_local / (2.0_dp - (2.0_dp * alphastar_local))), df=nu2_v))))
  end function past

  pure elemental function varast(p, nu1, nu2, alpha, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: nu1, nu2, alpha
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: nu1_v
    real(dp) :: nu2_v
    real(dp) :: alpha_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    real(dp) :: k1_local
    real(dp) :: k2_local
    real(dp) :: alphastar_local
    real(dp) :: minp_local
    real(dp) :: maxp_local
    nu1_v = 1.0_dp
    if (present(nu1)) nu1_v = nu1
    nu2_v = 1.0_dp
    if (present(nu2)) nu2_v = nu2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    k1_local = (gamma_fn(((nu1_v + 1.0_dp) / 2.0_dp)) / (sqrt((nu1_v / 2.0_dp)) * gamma_fn((nu1_v / 2.0_dp))))
    k2_local = (gamma_fn(((nu2_v + 1.0_dp) / 2.0_dp)) / (sqrt((nu2_v / 2.0_dp)) * gamma_fn((nu2_v / 2.0_dp))))
    alphastar_local = ((alpha_v * k1_local) / ((alpha_v * k1_local) + ((1.0_dp - alpha_v) * k2_local)))
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = pp
    minp_local = pp
    maxp_local = pp
    minp_local = min(pp, alpha_v)
    maxp_local = max(pp, alpha_v)
    res = (((2.0_dp * alphastar_local) * student_t_quantile((minp_local / (2.0_dp * alpha_v)), df=nu1_v)) + &
      & ((2.0_dp * (1.0_dp - alphastar_local)) * student_t_quantile((((maxp_local + 1.0_dp) - (2.0_dp * alpha_v)) / &
      & (2.0_dp - (2.0_dp * alpha_v))), df=nu2_v)))
  end function varast

  pure elemental function esast(p, nu1, nu2, alpha) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: nu1, nu2, alpha
    real(dp) :: res
    real(dp) :: nu1_v
    real(dp) :: nu2_v
    real(dp) :: alpha_v
    integer :: i
    real(dp) :: s_quad
    nu1_v = 1.0_dp
    if (present(nu1)) nu1_v = nu1
    nu2_v = 1.0_dp
    if (present(nu2)) nu2_v = nu2
    alpha_v = 0.5_dp
    if (present(alpha)) alpha_v = alpha
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varast(p * s_quad**3, nu1_v, nu2_v, alpha_v)
    end do
  end function esast

  pure elemental function dhalft(x, n, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: n
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: n_v
    logical :: log_pdf_v
    n_v = 1.0_dp
    if (present(n)) n_v = n
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (2.0_dp * student_t_pdf(x, df=n_v))
    if (((log_pdf_v))) res = (log(2.0_dp) + student_t_pdf(x, df=n_v, log_pdf=.true.))
  end function dhalft

  pure elemental function phalft(x, n, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: n
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: n_v
    logical :: log_p_v
    logical :: lower_tail_v
    n_v = 1.0_dp
    if (present(n)) n_v = n
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = beta_cdf(((x * x) / (n_v + (x * x))), shape1=0.5_dp, shape2=(n_v / 2.0_dp), log_p=log_p_v, &
      & lower_tail=lower_tail_v)
  end function phalft

  pure elemental function varhalft(p, n, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: n
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: n_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    n_v = 1.0_dp
    if (present(n)) n_v = n
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = sqrt(((n_v * beta_quantile(pp, shape1=0.5_dp, shape2=(n_v / 2.0_dp))) / (1.0_dp - beta_quantile(pp, &
      & shape1=0.5_dp, shape2=(n_v / 2.0_dp)))))
  end function varhalft

  pure elemental function eshalft(p, n) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: n
    real(dp) :: res
    real(dp) :: n_v
    integer :: i
    real(dp) :: s_quad
    n_v = 1.0_dp
    if (present(n)) n_v = n
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhalft(p * s_quad**3, n_v)
    end do
  end function eshalft

  pure elemental function dcauchy(x, mu, sigma, log_pdf) result(res)
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
    res = cauchy_pdf(x, location=mu_v, scale=sigma_v, log_pdf=log_pdf_v)
  end function dcauchy

  pure elemental function pcauchy(x, mu, sigma, log_p, lower_tail) result(res)
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
    res = cauchy_cdf(x, location=mu_v, scale=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function pcauchy

  pure elemental function varcauchy(p, mu, sigma, log_p, lower_tail) result(res)
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
    res = cauchy_quantile(pp, location=mu_v, scale=sigma_v, log_p=log_p_v, lower_tail=lower_tail_v)
  end function varcauchy

  pure elemental function escauchy(p, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varcauchy(p * s_quad**3, mu_v, sigma_v)
    end do
  end function escauchy

  pure elemental function dlogcauchy(x, mu, sigma, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((1.0_dp / x) * cauchy_pdf(log(x), location=mu_v, scale=sigma_v))
    if (((log_pdf_v))) res = (cauchy_pdf(log(x), location=mu_v, scale=sigma_v, log_pdf=.true.) - log(x))
  end function dlogcauchy

  pure elemental function plogcauchy(x, mu, sigma, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((1.0_dp / pi) * atan(((log(x) - mu_v) / sigma_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((1.0_dp / pi) * atan(((log(x) - mu_v) / &
      & sigma_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log(((1.0_dp / pi) * atan(((log(x) - mu_v) / sigma_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((1.0_dp / pi) * atan(((log(x) - mu_v) / &
      & sigma_v)))))
  end function plogcauchy

  pure elemental function varlogcauchy(p, mu, sigma, log_p, lower_tail) result(res)
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
    res = exp((mu_v + (sigma_v * tan((pi * pp)))))
  end function varlogcauchy

  pure elemental function eslogcauchy(p, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlogcauchy(p * s_quad**3, mu_v, sigma_v)
    end do
  end function eslogcauchy

  pure elemental function dhalfcauchy(x, sigma, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: sigma
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: sigma_v
    logical :: log_pdf_v
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if (((.not. log_pdf_v))) res = (2.0_dp * cauchy_pdf(x, scale=sigma_v))
    if (((log_pdf_v))) res = (log(2.0_dp) + cauchy_pdf(x, scale=sigma_v, log_pdf=.true.))
  end function dhalfcauchy

  pure elemental function phalfcauchy(x, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)))) res = ((2.0_dp / pi) * atan((x / sigma_v)))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)))) res = (1.0_dp - ((2.0_dp / pi) * atan((x / sigma_v))))
    if ((((log_p_v)) .and. ((lower_tail_v)))) res = log(((2.0_dp / pi) * atan((x / sigma_v))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)))) res = log((1.0_dp - ((2.0_dp / pi) * atan((x / sigma_v)))))
  end function phalfcauchy

  pure elemental function varhalfcauchy(p, sigma, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: sigma
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: sigma_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    pp = p
    if (log_p_v) pp = exp(pp)
    if (.not. lower_tail_v) pp = 1.0_dp - pp
    res = (sigma_v * tan(((pi * pp) / 2.0_dp)))
  end function varhalfcauchy

  pure elemental function eshalfcauchy(p, sigma) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: sigma
    real(dp) :: res
    real(dp) :: sigma_v
    integer :: i
    real(dp) :: s_quad
    sigma_v = 1.0_dp
    if (present(sigma)) sigma_v = sigma
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhalfcauchy(p * s_quad**3, sigma_v)
    end do
  end function eshalfcauchy

  pure elemental function dlaplace(x, mu, sigma, log_pdf) result(res)
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
    if (((.not. log_pdf_v))) res = ((1.0_dp / (2.0_dp * sigma_v)) * exp(((-abs((x - mu_v))) / sigma_v)))
    if (((log_pdf_v))) res = (((-abs((x - mu_v))) / sigma_v) - log((2.0_dp * sigma_v)))
  end function dlaplace

  pure elemental function plaplace(x, mu, sigma, log_p, lower_tail) result(res)
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x < mu_v)))) res = (0.5_dp * exp(((x - mu_v) / sigma_v)))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x >= mu_v)))) res = (1.0_dp - (0.5_dp * exp(((mu_v - x) / &
      & sigma_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x < mu_v)))) res = (1.0_dp - (0.5_dp * exp(((x - &
      & mu_v) / sigma_v))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x >= mu_v)))) res = (0.5_dp * exp(((mu_v - x) / &
      & sigma_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x < mu_v)))) res = (((x - mu_v) / sigma_v) - log(2.0_dp))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x >= mu_v)))) res = log((1.0_dp - (0.5_dp * exp(((mu_v - x) / &
      & sigma_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x < mu_v)))) res = log((1.0_dp - (0.5_dp * exp(((x - &
      & mu_v) / sigma_v)))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x >= mu_v)))) res = (((mu_v - x) / sigma_v) - log(2.0_dp))
  end function plaplace

  pure elemental function varlaplace(p, mu, sigma, log_p, lower_tail) result(res)
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
    res = pp
    if (((pp < 0.5_dp))) res = (mu_v + (sigma_v * log((2.0_dp * pp))))
    if (((pp >= 0.5_dp))) res = (mu_v - (sigma_v * log((2.0_dp * (1.0_dp - pp)))))
  end function varlaplace

  pure elemental function eslaplace(p, mu, sigma) result(res)
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varlaplace(p * s_quad**3, mu_v, sigma_v)
    end do
  end function eslaplace

  pure elemental function dpctalaplace(x, a, theta, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_pdf_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= theta_v)))) res = ((a_v * (1.0_dp - a_v)) * exp(((1.0_dp - a_v) * (x - &
      & theta_v))))
    if ((((.not. log_pdf_v)) .and. ((x > theta_v)))) res = ((a_v * (1.0_dp - a_v)) * exp((a_v * (theta_v - x))))
    if ((((log_pdf_v)) .and. ((x <= theta_v)))) res = (log((a_v * (1.0_dp - a_v))) + ((1.0_dp - a_v) * (x - theta_v)))
    if ((((log_pdf_v)) .and. ((x > theta_v)))) res = (log((a_v * (1.0_dp - a_v))) + (a_v * (theta_v - x)))
  end function dpctalaplace

  pure elemental function ppctalaplace(x, a, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    log_p_v = .false.
    if (present(log_p)) log_p_v = log_p
    lower_tail_v = .true.
    if (present(lower_tail)) lower_tail_v = lower_tail
    res = x
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (a_v * exp(((1.0_dp - a_v) * (x - &
      & theta_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = (1.0_dp - ((1.0_dp - a_v) * &
      & exp((a_v * (theta_v - x)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = (1.0_dp - (a_v * &
      & exp(((1.0_dp - a_v) * (x - theta_v)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = ((1.0_dp - a_v) * exp((a_v * &
      & (theta_v - x))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (log(a_v) + ((1.0_dp - a_v) * (x - &
      & theta_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = log((1.0_dp - ((1.0_dp - a_v) * exp((a_v &
      & * (theta_v - x))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = log((1.0_dp - (a_v * exp(((1.0_dp &
      & - a_v) * (x - theta_v))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = (log((1.0_dp - a_v)) + (a_v * &
      & (theta_v - x)))
  end function ppctalaplace

  pure elemental function varpctalaplace(p, a, theta, log_p, lower_tail) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    logical :: log_p_v
    logical :: lower_tail_v
    real(dp) :: pp
    a_v = 0.5_dp
    if (present(a)) a_v = a
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
    if (((pp <= a_v))) res = (theta_v + (log((pp / a_v)) / (1.0_dp - a_v)))
    if (((pp > a_v))) res = (theta_v - (log(((1.0_dp - pp) / (1.0_dp - a_v))) / a_v))
  end function varpctalaplace

  pure elemental function espctalaplace(p, a, theta) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.5_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    if (p <= 0.0_dp .or. p > 1.0_dp) then
      res = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    res = 0.0_dp
    do i = 1, gl_n
      s_quad = 0.5_dp * (gl_x(i) + 1.0_dp)
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varpctalaplace(p * s_quad**3, a_v, theta_v)
    end do
  end function espctalaplace

  pure elemental function dhblaplace(x, a, theta, phi, log_pdf) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta, phi
    logical, intent(in), optional :: log_pdf
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    logical :: log_pdf_v
    a_v = 0.5_dp
    if (present(a)) a_v = a
    theta_v = 0.0_dp
    if (present(theta)) theta_v = theta
    phi_v = 1.0_dp
    if (present(phi)) phi_v = phi
    log_pdf_v = .false.
    if (present(log_pdf)) log_pdf_v = log_pdf
    res = x
    if ((((.not. log_pdf_v)) .and. ((x <= theta_v)))) res = ((a_v * phi_v) * exp((phi_v * (x - theta_v))))
    if ((((.not. log_pdf_v)) .and. ((x > theta_v)))) res = (((1.0_dp - a_v) * phi_v) * exp((phi_v * (theta_v - x))))
    if ((((log_pdf_v)) .and. ((x <= theta_v)))) res = (log((a_v * phi_v)) + (phi_v * (x - theta_v)))
    if ((((log_pdf_v)) .and. ((x > theta_v)))) res = (log(((1.0_dp - a_v) * phi_v)) + (phi_v * (theta_v - x)))
  end function dhblaplace

  pure elemental function phblaplace(x, a, theta, phi, log_p, lower_tail) result(res)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: a, theta, phi
    logical, intent(in), optional :: log_p, lower_tail
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    logical :: log_p_v
    logical :: lower_tail_v
    a_v = 0.5_dp
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
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (a_v * exp((phi_v * (x - theta_v))))
    if ((((.not. log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = (1.0_dp - ((1.0_dp - a_v) * &
      & exp((phi_v * (theta_v - x)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = (1.0_dp - (a_v * exp((phi_v &
      & * (x - theta_v)))))
    if ((((.not. log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = ((1.0_dp - a_v) * exp((phi_v &
      & * (theta_v - x))))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x <= theta_v)))) res = (log(a_v) + (phi_v * (x - theta_v)))
    if ((((log_p_v)) .and. ((lower_tail_v)) .and. ((x > theta_v)))) res = log((1.0_dp - ((1.0_dp - a_v) * &
      & exp((phi_v * (theta_v - x))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x <= theta_v)))) res = log((1.0_dp - (a_v * exp((phi_v * &
      & (x - theta_v))))))
    if ((((log_p_v)) .and. ((.not. lower_tail_v)) .and. ((x > theta_v)))) res = (log((1.0_dp - a_v)) + (phi_v * &
      & (theta_v - x)))
  end function phblaplace

  pure elemental function varhblaplace(p, a, theta, phi, log_p, lower_tail) result(res)
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
    a_v = 0.5_dp
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
    if (((pp <= a_v))) res = (theta_v + ((1.0_dp / phi_v) * log((pp / a_v))))
    if (((pp > a_v))) res = (theta_v - ((1.0_dp / phi_v) * log(((1.0_dp - pp) / (1.0_dp - a_v)))))
  end function varhblaplace

  pure elemental function eshblaplace(p, a, theta, phi) result(res)
    real(dp), intent(in) :: p
    real(dp), intent(in), optional :: a, theta, phi
    real(dp) :: res
    real(dp) :: a_v
    real(dp) :: theta_v
    real(dp) :: phi_v
    integer :: i
    real(dp) :: s_quad
    a_v = 0.5_dp
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
      res = res + 1.5_dp * gl_w(i) * s_quad**2 * varhblaplace(p * s_quad**3, a_v, theta_v, phi_v)
    end do
  end function eshblaplace

end module vares_distributions_05
