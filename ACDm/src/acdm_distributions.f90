! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_distributions
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use acdm_kinds, only : dp, pi, tiny_pos
  use acdm_math, only : rng_state, random_uniform, random_normal, random_gamma, &
                        normal_cdf, gamma_p, beta_i, gamma_quantile, &
                        beta_quantile
  implicit none
  private

  integer, parameter, public :: DIST_EXPONENTIAL = 1
  integer, parameter, public :: DIST_WEIBULL = 2
  integer, parameter, public :: DIST_BURR = 3
  integer, parameter, public :: DIST_GENGAMMA = 4
  integer, parameter, public :: DIST_GENF = 5
  integer, parameter, public :: DIST_QWEIBULL = 6
  integer, parameter, public :: DIST_MIXQWE = 7
  integer, parameter, public :: DIST_MIXQWW = 8
  integer, parameter, public :: DIST_MIXINVGAUSS = 9
  integer, parameter, public :: DIST_BIRNBAUM_SAUNDERS = 10

  public :: distribution_parameter_count, distribution_pdf
  public :: distribution_logpdf, distribution_cdf, distribution_quantile
  public :: sample_distribution, forced_scale
  public :: dburr, pburr, qburr, rburr, burr_expectation
  public :: dgenf, pgenf, qgenf, rgenf, genf_hazard
  public :: dgengamma, pgengamma, qgengamma, rgengamma, gengamma_hazard
  public :: dqweibull, pqweibull, qqweibull, rqweibull
  public :: qweibull_expectation, qweibull_hazard
  public :: dmixqwe, pmixqwe, qmixqwe, rmixqwe, mixqwe_hazard
  public :: dmixqww, pmixqww, qmixqww, rmixqww, mixqww_hazard
  public :: dmixinvgauss, pmixinvgauss, qmixinvgauss
  public :: rmixinvgauss, mixinvgauss_hazard
  public :: dbirnbaum_saunders, pbirnbaum_saunders
  public :: qbirnbaum_saunders, rbirnbaum_saunders

contains

  pure integer function distribution_parameter_count(code) result(n)
    integer, intent(in) :: code

    select case (code)
    case (DIST_EXPONENTIAL)
      n = 0
    case (DIST_WEIBULL)
      n = 1
    case (DIST_BURR, DIST_GENGAMMA, DIST_QWEIBULL)
      n = 2
    case (DIST_GENF, DIST_MIXINVGAUSS)
      n = 3
    case (DIST_MIXQWE)
      n = 4
    case (DIST_MIXQWW)
      n = 5
    case (DIST_BIRNBAUM_SAUNDERS)
      n = 1
    case default
      n = -1
    end select
  end function distribution_parameter_count

  pure function nan_value() result(x)
    real(dp) :: x
    x = ieee_value(x, ieee_quiet_nan)
  end function nan_value

  pure function weibull_scale(shape, force_mean) result(scale)
    real(dp), intent(in) :: shape
    logical, intent(in) :: force_mean
    real(dp) :: scale

    if (shape <= 0.0_dp) then
      scale = nan_value()
    else if (force_mean) then
      scale = 1.0_dp / gamma(1.0_dp + 1.0_dp / shape)
    else
      scale = 1.0_dp
    end if
  end function weibull_scale

  pure function burr_theta(kappa, sig2, force_mean) result(theta)
    real(dp), intent(in) :: kappa, sig2
    logical, intent(in) :: force_mean
    real(dp) :: theta, lm

    if (kappa <= 0.0_dp .or. sig2 <= 0.0_dp) then
      theta = nan_value()
      return
    end if
    if (.not. force_mean) then
      theta = 1.0_dp
      return
    end if
    if (1.0_dp / sig2 <= 1.0_dp / kappa) then
      theta = nan_value()
      return
    end if
    lm = log_gamma(1.0_dp + 1.0_dp / kappa) + &
         log_gamma(1.0_dp / sig2 - 1.0_dp / kappa) - &
         (1.0_dp + 1.0_dp / kappa) * log(sig2) - &
         log_gamma(1.0_dp / sig2 + 1.0_dp)
    theta = exp(kappa * lm)
  end function burr_theta

  pure function gengamma_lambda(kappa, gamma_shape, force_mean) result(lambda)
    real(dp), intent(in) :: kappa, gamma_shape
    logical, intent(in) :: force_mean
    real(dp) :: lambda

    if (kappa <= 0.0_dp .or. gamma_shape <= 0.0_dp) then
      lambda = nan_value()
    else if (force_mean) then
      lambda = exp(log_gamma(kappa) - &
                   log_gamma(kappa + 1.0_dp / gamma_shape))
    else
      lambda = 1.0_dp
    end if
  end function gengamma_lambda

  pure function genf_lambda(kappa, eta, gamma_shape, force_mean) result(lambda)
    real(dp), intent(in) :: kappa, eta, gamma_shape
    logical, intent(in) :: force_mean
    real(dp) :: lambda

    if (min(kappa, eta, gamma_shape) <= 0.0_dp) then
      lambda = nan_value()
    else if (.not. force_mean) then
      lambda = 1.0_dp
    else if (eta <= 1.0_dp / gamma_shape) then
      lambda = nan_value()
    else
      lambda = exp(log_gamma(kappa) + log_gamma(eta) - &
                   log_gamma(kappa + 1.0_dp / gamma_shape) - &
                   log_gamma(eta - 1.0_dp / gamma_shape) - &
                   log(eta) / gamma_shape)
    end if
  end function genf_lambda

  pure function qweibull_mean_factor(a, qdist) result(factor)
    real(dp), intent(in) :: a, qdist
    real(dp) :: factor, m

    if (a <= 0.0_dp .or. qdist >= 2.0_dp) then
      factor = nan_value()
      return
    end if
    if (abs(qdist - 1.0_dp) < 1.0e-10_dp) then
      factor = gamma(1.0_dp + 1.0_dp / a)
    else if (qdist < 1.0_dp) then
      m = (2.0_dp - qdist) / (1.0_dp - qdist)
      factor = exp(log_gamma(1.0_dp + 1.0_dp / a) + &
                   log_gamma(m + 1.0_dp) - &
                   log_gamma(m + 1.0_dp + 1.0_dp / a)) / &
               (1.0_dp - qdist)**(1.0_dp / a)
    else
      m = (2.0_dp - qdist) / (qdist - 1.0_dp)
      if (m <= 1.0_dp / a) then
        factor = nan_value()
      else
        factor = exp(log_gamma(1.0_dp + 1.0_dp / a) + &
                     log_gamma(m - 1.0_dp / a) - log_gamma(m)) / &
                 (qdist - 1.0_dp)**(1.0_dp / a)
      end if
    end if
  end function qweibull_mean_factor

  pure function qweibull_b(a, qdist, force_mean) result(b)
    real(dp), intent(in) :: a, qdist
    logical, intent(in) :: force_mean
    real(dp) :: b, factor

    if (.not. force_mean) then
      b = 1.0_dp
    else
      factor = qweibull_mean_factor(a, qdist)
      if (.not. (factor > 0.0_dp)) then
        b = nan_value()
      else
        b = 1.0_dp / factor
      end if
    end if
  end function qweibull_b

  pure function mixqwe_b(pdist, a, qdist, lambda, force_mean) result(b)
    real(dp), intent(in) :: pdist, a, qdist, lambda
    logical, intent(in) :: force_mean
    real(dp) :: b, f, numerator

    if (.not. force_mean) then
      b = 1.0_dp
      return
    end if
    f = qweibull_mean_factor(a, qdist)
    numerator = 1.0_dp - (1.0_dp - pdist) * lambda
    if (pdist <= 0.0_dp .or. pdist >= 1.0_dp .or. &
        f <= 0.0_dp .or. numerator <= 0.0_dp) then
      b = nan_value()
    else
      b = numerator / (pdist * f)
    end if
  end function mixqwe_b

  pure function mixqww_b(pdist, a, qdist, theta, gamma_shape, &
                         force_mean) result(b)
    real(dp), intent(in) :: pdist, a, qdist, theta, gamma_shape
    logical, intent(in) :: force_mean
    real(dp) :: b, f, wmean, numerator

    if (.not. force_mean) then
      b = 1.0_dp
      return
    end if
    f = qweibull_mean_factor(a, qdist)
    if (theta <= 0.0_dp .or. gamma_shape <= 0.0_dp) then
      b = nan_value()
      return
    end if
    wmean = theta**(-1.0_dp / gamma_shape) * &
            gamma(1.0_dp + 1.0_dp / gamma_shape)
    numerator = 1.0_dp - (1.0_dp - pdist) * wmean
    if (pdist <= 0.0_dp .or. pdist >= 1.0_dp .or. &
        f <= 0.0_dp .or. numerator <= 0.0_dp) then
      b = nan_value()
    else
      b = numerator / (pdist * f)
    end if
  end function mixqww_b

  pure function mixinvgauss_phi(theta, lambda, gamma_mix, &
                                force_mean) result(phi)
    real(dp), intent(in) :: theta, lambda, gamma_mix
    logical, intent(in) :: force_mean
    real(dp) :: phi

    if (theta <= 0.0_dp .or. lambda <= 0.0_dp .or. gamma_mix < 0.0_dp) then
      phi = nan_value()
    else if (force_mean) then
      phi = theta * (1.0_dp + theta**2 / lambda / (1.0_dp + gamma_mix))
    else
      phi = 1.0_dp
    end if
  end function mixinvgauss_phi

  pure function forced_scale(code, para) result(scale)
    integer, intent(in) :: code
    real(dp), intent(in) :: para(:)
    real(dp) :: scale

    select case (code)
    case (DIST_EXPONENTIAL)
      scale = 1.0_dp
    case (DIST_WEIBULL)
      scale = weibull_scale(para(1), .true.)
    case (DIST_BURR)
      scale = burr_theta(para(1), para(2), .true.)
    case (DIST_GENGAMMA)
      scale = gengamma_lambda(para(1), para(2), .true.)
    case (DIST_GENF)
      scale = genf_lambda(para(1), para(2), para(3), .true.)
    case (DIST_QWEIBULL)
      scale = qweibull_b(para(1), para(2), .true.)
    case (DIST_MIXQWE)
      scale = mixqwe_b(para(1), para(2), para(3), para(4), .true.)
    case (DIST_MIXQWW)
      scale = mixqww_b(para(1), para(2), para(3), para(4), &
                       para(5), .true.)
    case (DIST_MIXINVGAUSS)
      scale = mixinvgauss_phi(para(1), para(2), para(3), .true.)
    case (DIST_BIRNBAUM_SAUNDERS)
      scale = 1.0_dp / (1.0_dp + 0.5_dp * para(1)**2)
    case default
      scale = nan_value()
    end select
  end function forced_scale

  pure elemental function dburr(x, theta, kappa, sig2) result(f)
    real(dp), intent(in) :: x, theta, kappa, sig2
    real(dp) :: f

    if (x <= 0.0_dp .or. min(theta, kappa, sig2) <= 0.0_dp) then
      f = 0.0_dp
    else
      f = theta * kappa * x**(kappa - 1.0_dp) / &
          (1.0_dp + sig2 * theta * x**kappa)**(1.0_dp / sig2 + 1.0_dp)
    end if
  end function dburr

  pure elemental function pburr(x, theta, kappa, sig2) result(p)
    real(dp), intent(in) :: x, theta, kappa, sig2
    real(dp) :: p

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (min(theta, kappa, sig2) <= 0.0_dp) then
      p = nan_value()
    else
      p = 1.0_dp - (1.0_dp + sig2 * theta * x**kappa)**(-1.0_dp / sig2)
    end if
  end function pburr

  pure elemental function qburr(p, theta, kappa, sig2) result(x)
    real(dp), intent(in) :: p, theta, kappa, sig2
    real(dp) :: x

    if (p <= 0.0_dp) then
      x = 0.0_dp
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (min(theta, kappa, sig2) <= 0.0_dp) then
      x = nan_value()
    else
      x = (((1.0_dp - p)**(-sig2) - 1.0_dp) / &
           (sig2 * theta))**(1.0_dp / kappa)
    end if
  end function qburr

  pure function burr_expectation(theta, kappa, sig2) result(value)
    real(dp), intent(in) :: theta, kappa, sig2
    real(dp) :: value

    if (min(theta, kappa, sig2) <= 0.0_dp .or. &
        1.0_dp / sig2 <= 1.0_dp / kappa) then
      value = nan_value()
    else
      value = exp(-log(theta) / kappa + &
                  log_gamma(1.0_dp + 1.0_dp / kappa) + &
                  log_gamma(1.0_dp / sig2 - 1.0_dp / kappa) - &
                  (1.0_dp + 1.0_dp / kappa) * log(sig2) - &
                  log_gamma(1.0_dp / sig2 + 1.0_dp))
    end if
  end function burr_expectation

  function rburr(rng, theta, kappa, sig2) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: theta, kappa, sig2
    real(dp) :: x
    x = qburr(random_uniform(rng), theta, kappa, sig2)
  end function rburr

  pure elemental function dgengamma(x, gamma_shape, kappa, lambda) result(f)
    real(dp), intent(in) :: x, gamma_shape, kappa, lambda
    real(dp) :: f, logf

    if (x <= 0.0_dp .or. min(gamma_shape, kappa, lambda) <= 0.0_dp) then
      f = 0.0_dp
    else
      logf = (kappa * gamma_shape - 1.0_dp) * log(x) - &
             kappa * gamma_shape * log(lambda) - log_gamma(kappa) + &
             log(gamma_shape) - (x / lambda)**gamma_shape
      f = exp(logf)
    end if
  end function dgengamma

  pure elemental function pgengamma(x, gamma_shape, kappa, lambda) result(p)
    real(dp), intent(in) :: x, gamma_shape, kappa, lambda
    real(dp) :: p

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (min(gamma_shape, kappa, lambda) <= 0.0_dp) then
      p = nan_value()
    else
      p = gamma_p(kappa, (x / lambda)**gamma_shape)
    end if
  end function pgengamma

  function qgengamma(p, gamma_shape, kappa, lambda) result(x)
    real(dp), intent(in) :: p, gamma_shape, kappa, lambda
    real(dp) :: x

    if (min(gamma_shape, kappa, lambda) <= 0.0_dp) then
      x = nan_value()
    else
      x = lambda * gamma_quantile(p, kappa)**(1.0_dp / gamma_shape)
    end if
  end function qgengamma

  function rgengamma(rng, gamma_shape, kappa, lambda) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: gamma_shape, kappa, lambda
    real(dp) :: x

    if (min(gamma_shape, kappa, lambda) <= 0.0_dp) then
      x = nan_value()
    else
      x = lambda * random_gamma(rng, kappa)**(1.0_dp / gamma_shape)
    end if
  end function rgengamma

  pure elemental function gengamma_hazard(x, gamma_shape, kappa, lambda) &
      result(h)
    real(dp), intent(in) :: x, gamma_shape, kappa, lambda
    real(dp) :: h, s

    s = 1.0_dp - pgengamma(x, gamma_shape, kappa, lambda)
    if (s <= 0.0_dp) then
      h = huge(1.0_dp)
    else
      h = dgengamma(x, gamma_shape, kappa, lambda) / s
    end if
  end function gengamma_hazard

  pure elemental function dgenf(x, kappa, eta, gamma_shape, lambda) result(f)
    real(dp), intent(in) :: x, kappa, eta, gamma_shape, lambda
    real(dp) :: f, logf

    if (x <= 0.0_dp .or. min(kappa, eta, gamma_shape, lambda) <= 0.0_dp) then
      f = 0.0_dp
    else
      logf = (kappa * gamma_shape - 1.0_dp) * log(x) - &
             (eta + kappa) * log(eta + (x / lambda)**gamma_shape) + &
             eta * log(eta) - &
             (log_gamma(kappa) + log_gamma(eta) - &
              log_gamma(kappa + eta)) - &
             kappa * gamma_shape * log(lambda) + log(gamma_shape)
      f = exp(logf)
    end if
  end function dgenf

  pure elemental function pgenf(x, kappa, eta, gamma_shape, lambda) result(p)
    real(dp), intent(in) :: x, kappa, eta, gamma_shape, lambda
    real(dp) :: p, t, y

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (min(kappa, eta, gamma_shape, lambda) <= 0.0_dp) then
      p = nan_value()
    else
      t = (x / lambda)**gamma_shape
      y = t / (eta + t)
      p = beta_i(y, kappa, eta)
    end if
  end function pgenf

  function qgenf(p, kappa, eta, gamma_shape, lambda) result(x)
    real(dp), intent(in) :: p, kappa, eta, gamma_shape, lambda
    real(dp) :: x, y

    if (min(kappa, eta, gamma_shape, lambda) <= 0.0_dp) then
      x = nan_value()
      return
    end if
    y = beta_quantile(p, kappa, eta)
    if (y >= 1.0_dp) then
      x = huge(1.0_dp)
    else
      x = lambda * (eta * y / max(tiny_pos, 1.0_dp - y))** &
          (1.0_dp / gamma_shape)
    end if
  end function qgenf

  function rgenf(rng, kappa, eta, gamma_shape, lambda) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: kappa, eta, gamma_shape, lambda
    real(dp) :: x, g1, g2

    g1 = random_gamma(rng, kappa)
    g2 = random_gamma(rng, eta)
    x = lambda * (eta * g1 / g2)**(1.0_dp / gamma_shape)
  end function rgenf

  pure elemental function genf_hazard(x, kappa, eta, gamma_shape, lambda) &
      result(h)
    real(dp), intent(in) :: x, kappa, eta, gamma_shape, lambda
    real(dp) :: h, s

    s = 1.0_dp - pgenf(x, kappa, eta, gamma_shape, lambda)
    if (s <= 0.0_dp) then
      h = huge(1.0_dp)
    else
      h = dgenf(x, kappa, eta, gamma_shape, lambda) / s
    end if
  end function genf_hazard

  pure elemental function dqweibull(x, a, qdist, b) result(f)
    real(dp), intent(in) :: x, a, qdist, b
    real(dp) :: f, base

    if (x <= 0.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp .or. &
        qdist >= 2.0_dp) then
      f = 0.0_dp
      return
    end if
    if (abs(qdist - 1.0_dp) < 1.0e-10_dp) then
      f = a / b * (x / b)**(a - 1.0_dp) * exp(-(x / b)**a)
      return
    end if
    base = 1.0_dp - (1.0_dp - qdist) * (x / b)**a
    if (base <= 0.0_dp) then
      f = 0.0_dp
    else
      f = (2.0_dp - qdist) * a / b**a * x**(a - 1.0_dp) * &
          base**(1.0_dp / (1.0_dp - qdist))
    end if
  end function dqweibull

  pure elemental function pqweibull(x, a, qdist, b) result(p)
    real(dp), intent(in) :: x, a, qdist, b
    real(dp) :: p, base

    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    else if (a <= 0.0_dp .or. b <= 0.0_dp .or. qdist >= 2.0_dp) then
      p = nan_value()
      return
    end if
    if (abs(qdist - 1.0_dp) < 1.0e-10_dp) then
      p = 1.0_dp - exp(-(x / b)**a)
      return
    end if
    base = 1.0_dp - (1.0_dp - qdist) * (x / b)**a
    if (base <= 0.0_dp) then
      p = 1.0_dp
    else
      p = 1.0_dp - base**((2.0_dp - qdist) / (1.0_dp - qdist))
      p = max(0.0_dp, min(1.0_dp, p))
    end if
  end function pqweibull

  pure elemental function qqweibull(p, a, qdist, b) result(x)
    real(dp), intent(in) :: p, a, qdist, b
    real(dp) :: x, inside

    if (p <= 0.0_dp) then
      x = 0.0_dp
    else if (p >= 1.0_dp) then
      if (qdist < 1.0_dp) then
        x = b / (1.0_dp - qdist)**(1.0_dp / a)
      else
        x = huge(1.0_dp)
      end if
    else if (a <= 0.0_dp .or. b <= 0.0_dp .or. qdist >= 2.0_dp) then
      x = nan_value()
    else if (abs(qdist - 1.0_dp) < 1.0e-10_dp) then
      x = b * (-log(1.0_dp - p))**(1.0_dp / a)
    else
      inside = (1.0_dp - (1.0_dp - p)** &
               ((1.0_dp - qdist) / (2.0_dp - qdist))) / &
               (1.0_dp - qdist)
      x = b * max(0.0_dp, inside)**(1.0_dp / a)
    end if
  end function qqweibull

  function rqweibull(rng, a, qdist, b) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: a, qdist, b
    real(dp) :: x
    x = qqweibull(random_uniform(rng), a, qdist, b)
  end function rqweibull

  pure function qweibull_expectation(a, qdist, b) result(value)
    real(dp), intent(in) :: a, qdist, b
    real(dp) :: value, f

    f = qweibull_mean_factor(a, qdist)
    if (b <= 0.0_dp .or. .not. (f > 0.0_dp)) then
      value = nan_value()
    else
      value = b * f
    end if
  end function qweibull_expectation

  pure elemental function qweibull_hazard(x, a, qdist, b) result(h)
    real(dp), intent(in) :: x, a, qdist, b
    real(dp) :: h, s

    s = 1.0_dp - pqweibull(x, a, qdist, b)
    if (s <= 0.0_dp) then
      h = huge(1.0_dp)
    else
      h = dqweibull(x, a, qdist, b) / s
    end if
  end function qweibull_hazard

  pure elemental function dmixqwe(x, pdist, a, qdist, lambda, b) result(f)
    real(dp), intent(in) :: x, pdist, a, qdist, lambda, b
    real(dp) :: f

    if (pdist < 0.0_dp .or. pdist > 1.0_dp .or. lambda <= 0.0_dp) then
      f = nan_value()
    else
      f = pdist * dqweibull(x, a, qdist, b)
      if (x >= 0.0_dp) f = f + (1.0_dp - pdist) / lambda * exp(-x / lambda)
    end if
  end function dmixqwe

  pure elemental function pmixqwe(x, pdist, a, qdist, lambda, b) result(p)
    real(dp), intent(in) :: x, pdist, a, qdist, lambda, b
    real(dp) :: p

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = pdist * pqweibull(x, a, qdist, b) + &
          (1.0_dp - pdist) * (1.0_dp - exp(-x / lambda))
    end if
  end function pmixqwe

  function qmixqwe(p, pdist, a, qdist, lambda, b) result(x)
    real(dp), intent(in) :: p, pdist, a, qdist, lambda, b
    real(dp) :: x, lo, hi, mid
    integer :: i

    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = 0.0_dp
    hi = max(lambda, b)
    do while (pmixqwe(hi, pdist, a, qdist, lambda, b) < p)
      hi = 2.0_dp * hi
      if (hi > 1.0e100_dp) exit
    end do
    do i = 1, 140
      mid = 0.5_dp * (lo + hi)
      if (pmixqwe(mid, pdist, a, qdist, lambda, b) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function qmixqwe

  function rmixqwe(rng, pdist, a, qdist, lambda, b) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: pdist, a, qdist, lambda, b
    real(dp) :: x

    if (random_uniform(rng) < pdist) then
      x = rqweibull(rng, a, qdist, b)
    else
      x = -lambda * log(random_uniform(rng))
    end if
  end function rmixqwe

  pure elemental function mixqwe_hazard(x, pdist, a, qdist, lambda, b) &
      result(h)
    real(dp), intent(in) :: x, pdist, a, qdist, lambda, b
    real(dp) :: h, s

    s = 1.0_dp - pmixqwe(x, pdist, a, qdist, lambda, b)
    if (s <= 0.0_dp) then
      h = huge(1.0_dp)
    else
      h = dmixqwe(x, pdist, a, qdist, lambda, b) / s
    end if
  end function mixqwe_hazard

  pure elemental function dmixqww(x, pdist, a, qdist, theta, &
                                  gamma_shape, b) result(f)
    real(dp), intent(in) :: x, pdist, a, qdist, theta, gamma_shape, b
    real(dp) :: f

    if (pdist < 0.0_dp .or. pdist > 1.0_dp .or. &
        theta <= 0.0_dp .or. gamma_shape <= 0.0_dp) then
      f = nan_value()
    else
      f = pdist * dqweibull(x, a, qdist, b)
      if (x > 0.0_dp) then
        f = f + (1.0_dp - pdist) * theta * gamma_shape * &
            x**(gamma_shape - 1.0_dp) * exp(-theta * x**gamma_shape)
      end if
    end if
  end function dmixqww

  pure elemental function pmixqww(x, pdist, a, qdist, theta, &
                                  gamma_shape, b) result(p)
    real(dp), intent(in) :: x, pdist, a, qdist, theta, gamma_shape, b
    real(dp) :: p

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else
      p = pdist * pqweibull(x, a, qdist, b) + &
          (1.0_dp - pdist) * (1.0_dp - exp(-theta * x**gamma_shape))
    end if
  end function pmixqww

  function qmixqww(p, pdist, a, qdist, theta, gamma_shape, b) result(x)
    real(dp), intent(in) :: p, pdist, a, qdist, theta, gamma_shape, b
    real(dp) :: x, lo, hi, mid
    integer :: i

    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = 0.0_dp
    hi = max(b, theta**(-1.0_dp / gamma_shape))
    do while (pmixqww(hi, pdist, a, qdist, theta, gamma_shape, b) < p)
      hi = 2.0_dp * hi
      if (hi > 1.0e100_dp) exit
    end do
    do i = 1, 140
      mid = 0.5_dp * (lo + hi)
      if (pmixqww(mid, pdist, a, qdist, theta, gamma_shape, b) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function qmixqww

  function rmixqww(rng, pdist, a, qdist, theta, gamma_shape, b) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: pdist, a, qdist, theta, gamma_shape, b
    real(dp) :: x

    if (random_uniform(rng) < pdist) then
      x = rqweibull(rng, a, qdist, b)
    else
      x = (-log(random_uniform(rng)) / theta)**(1.0_dp / gamma_shape)
    end if
  end function rmixqww

  pure elemental function mixqww_hazard(x, pdist, a, qdist, theta, &
                                        gamma_shape, b) result(h)
    real(dp), intent(in) :: x, pdist, a, qdist, theta, gamma_shape, b
    real(dp) :: h, s

    s = 1.0_dp - pmixqww(x, pdist, a, qdist, theta, gamma_shape, b)
    if (s <= 0.0_dp) then
      h = huge(1.0_dp)
    else
      h = dmixqww(x, pdist, a, qdist, theta, gamma_shape, b) / s
    end if
  end function mixqww_hazard

  pure elemental function dmixinvgauss(x, theta, lambda, gamma_mix, phi) &
      result(f)
    real(dp), intent(in) :: x, theta, lambda, gamma_mix, phi
    real(dp) :: f

    if (x <= 0.0_dp .or. min(theta, lambda, phi) <= 0.0_dp .or. &
        gamma_mix < 0.0_dp) then
      f = 0.0_dp
    else
      f = (gamma_mix + phi * x) / (gamma_mix + theta) * &
          sqrt(lambda / (2.0_dp * pi * x**3 * phi)) * &
          exp(-lambda * (phi * x - theta)**2 / &
              (2.0_dp * phi * x * theta**2))
    end if
  end function dmixinvgauss

  pure elemental function pmixinvgauss(x, theta, lambda, gamma_mix, phi) &
      result(p)
    real(dp), intent(in) :: x, theta, lambda, gamma_mix, phi
    real(dp) :: p, t1, t2, a1, a2

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (min(theta, lambda, phi) <= 0.0_dp .or. gamma_mix < 0.0_dp) then
      p = nan_value()
    else
      t1 = x / theta - 1.0_dp
      t2 = -x / theta - 1.0_dp
      a1 = t1 * sqrt(lambda / (x * phi))
      a2 = t2 * sqrt(lambda / (x * phi))
      p = normal_cdf(a1) + (gamma_mix - theta) / &
          (theta + gamma_mix) * normal_cdf(a2) * exp(2.0_dp * lambda / theta)
      p = max(0.0_dp, min(1.0_dp, p))
    end if
  end function pmixinvgauss

  function qmixinvgauss(p, theta, lambda, gamma_mix, phi) result(x)
    real(dp), intent(in) :: p, theta, lambda, gamma_mix, phi
    real(dp) :: x, lo, hi, mid
    integer :: i

    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    lo = 0.0_dp
    hi = max(theta / phi, 1.0_dp)
    do while (pmixinvgauss(hi, theta, lambda, gamma_mix, phi) < p)
      hi = 2.0_dp * hi
      if (hi > 1.0e100_dp) exit
    end do
    do i = 1, 140
      mid = 0.5_dp * (lo + hi)
      if (pmixinvgauss(mid, theta, lambda, gamma_mix, phi) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp * (lo + hi)
  end function qmixinvgauss

  function rmixinvgauss(rng, theta, lambda, gamma_mix, phi) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: theta, lambda, gamma_mix, phi
    real(dp) :: x
    x = qmixinvgauss(random_uniform(rng), theta, lambda, gamma_mix, phi)
  end function rmixinvgauss

  pure elemental function mixinvgauss_hazard(x, theta, lambda, gamma_mix, &
                                             phi) result(h)
    real(dp), intent(in) :: x, theta, lambda, gamma_mix, phi
    real(dp) :: h, s

    s = 1.0_dp - pmixinvgauss(x, theta, lambda, gamma_mix, phi)
    if (s <= 0.0_dp) then
      h = huge(1.0_dp)
    else
      h = dmixinvgauss(x, theta, lambda, gamma_mix, phi) / s
    end if
  end function mixinvgauss_hazard

  pure elemental function dbirnbaum_saunders(x, kappa, sigma) result(f)
    real(dp), intent(in) :: x, kappa, sigma
    real(dp) :: f

    if (x <= 0.0_dp .or. kappa <= 0.0_dp .or. sigma <= 0.0_dp) then
      f = 0.0_dp
    else
      f = ((sigma / x)**0.5_dp + (sigma / x)**1.5_dp) / &
          (2.0_dp * kappa * sigma * sqrt(2.0_dp * pi)) * &
          exp(-(x / sigma + sigma / x - 2.0_dp) / &
              (2.0_dp * kappa**2))
    end if
  end function dbirnbaum_saunders

  pure elemental function pbirnbaum_saunders(x, kappa, sigma) result(p)
    real(dp), intent(in) :: x, kappa, sigma
    real(dp) :: p

    if (x <= 0.0_dp) then
      p = 0.0_dp
    else if (kappa <= 0.0_dp .or. sigma <= 0.0_dp) then
      p = nan_value()
    else
      p = normal_cdf(((x / sigma)**0.5_dp - &
                      (sigma / x)**0.5_dp) / kappa)
    end if
  end function pbirnbaum_saunders

  pure elemental function qbirnbaum_saunders(p, kappa, sigma) result(x)
    real(dp), intent(in) :: p, kappa, sigma
    real(dp) :: x, z, t

    if (p <= 0.0_dp) then
      x = 0.0_dp
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else
      z = inverse_normal_local(p)
      t = 0.5_dp * kappa * z
      x = sigma * (t + sqrt(1.0_dp + t * t))**2
    end if
  end function qbirnbaum_saunders

  function rbirnbaum_saunders(rng, kappa, sigma) result(x)
    type(rng_state), intent(inout) :: rng
    real(dp), intent(in) :: kappa, sigma
    real(dp) :: x, z, t

    z = random_normal(rng)
    t = 0.5_dp * kappa * z
    x = sigma * (t + sqrt(1.0_dp + t * t))**2
  end function rbirnbaum_saunders

  pure elemental function inverse_normal_local(p) result(z)
    use acdm_math, only : normal_quantile
    real(dp), intent(in) :: p
    real(dp) :: z
    z = normal_quantile(p)
  end function inverse_normal_local

  pure function distribution_pdf(x, code, para, force_mean) result(f)
    real(dp), intent(in) :: x
    integer, intent(in) :: code
    real(dp), intent(in) :: para(:)
    logical, intent(in) :: force_mean
    real(dp) :: f, scale

    select case (code)
    case (DIST_EXPONENTIAL)
      if (x < 0.0_dp) then
        f = 0.0_dp
      else
        f = exp(-x)
      end if
    case (DIST_WEIBULL)
      scale = weibull_scale(para(1), force_mean)
      if (x <= 0.0_dp .or. .not. (scale > 0.0_dp)) then
        f = 0.0_dp
      else
        f = para(1) / scale * (x / scale)**(para(1) - 1.0_dp) * &
            exp(-(x / scale)**para(1))
      end if
    case (DIST_BURR)
      scale = burr_theta(para(1), para(2), force_mean)
      f = dburr(x, scale, para(1), para(2))
    case (DIST_GENGAMMA)
      scale = gengamma_lambda(para(1), para(2), force_mean)
      f = dgengamma(x, para(2), para(1), scale)
    case (DIST_GENF)
      scale = genf_lambda(para(1), para(2), para(3), force_mean)
      f = dgenf(x, para(1), para(2), para(3), scale)
    case (DIST_QWEIBULL)
      scale = qweibull_b(para(1), para(2), force_mean)
      f = dqweibull(x, para(1), para(2), scale)
    case (DIST_MIXQWE)
      scale = mixqwe_b(para(1), para(2), para(3), para(4), force_mean)
      f = dmixqwe(x, para(1), para(2), para(3), para(4), scale)
    case (DIST_MIXQWW)
      scale = mixqww_b(para(1), para(2), para(3), para(4), &
                       para(5), force_mean)
      f = dmixqww(x, para(1), para(2), para(3), para(4), para(5), scale)
    case (DIST_MIXINVGAUSS)
      scale = mixinvgauss_phi(para(1), para(2), para(3), force_mean)
      f = dmixinvgauss(x, para(1), para(2), para(3), scale)
    case (DIST_BIRNBAUM_SAUNDERS)
      scale = 1.0_dp
      if (force_mean) scale = 1.0_dp / (1.0_dp + 0.5_dp * para(1)**2)
      f = dbirnbaum_saunders(x, para(1), scale)
    case default
      f = nan_value()
    end select
  end function distribution_pdf

  pure function distribution_logpdf(x, code, para, force_mean) result(logf)
    real(dp), intent(in) :: x
    integer, intent(in) :: code
    real(dp), intent(in) :: para(:)
    logical, intent(in) :: force_mean
    real(dp) :: logf, f

    f = distribution_pdf(x, code, para, force_mean)
    if (f > 0.0_dp) then
      logf = log(f)
    else
      logf = -huge(1.0_dp)
    end if
  end function distribution_logpdf

  pure function distribution_cdf(x, code, para, force_mean) result(p)
    real(dp), intent(in) :: x
    integer, intent(in) :: code
    real(dp), intent(in) :: para(:)
    logical, intent(in) :: force_mean
    real(dp) :: p, scale

    select case (code)
    case (DIST_EXPONENTIAL)
      if (x <= 0.0_dp) then
        p = 0.0_dp
      else
        p = 1.0_dp - exp(-x)
      end if
    case (DIST_WEIBULL)
      scale = weibull_scale(para(1), force_mean)
      if (x <= 0.0_dp) then
        p = 0.0_dp
      else
        p = 1.0_dp - exp(-(x / scale)**para(1))
      end if
    case (DIST_BURR)
      scale = burr_theta(para(1), para(2), force_mean)
      p = pburr(x, scale, para(1), para(2))
    case (DIST_GENGAMMA)
      scale = gengamma_lambda(para(1), para(2), force_mean)
      p = pgengamma(x, para(2), para(1), scale)
    case (DIST_GENF)
      scale = genf_lambda(para(1), para(2), para(3), force_mean)
      p = pgenf(x, para(1), para(2), para(3), scale)
    case (DIST_QWEIBULL)
      scale = qweibull_b(para(1), para(2), force_mean)
      p = pqweibull(x, para(1), para(2), scale)
    case (DIST_MIXQWE)
      scale = mixqwe_b(para(1), para(2), para(3), para(4), force_mean)
      p = pmixqwe(x, para(1), para(2), para(3), para(4), scale)
    case (DIST_MIXQWW)
      scale = mixqww_b(para(1), para(2), para(3), para(4), &
                       para(5), force_mean)
      p = pmixqww(x, para(1), para(2), para(3), para(4), para(5), scale)
    case (DIST_MIXINVGAUSS)
      scale = mixinvgauss_phi(para(1), para(2), para(3), force_mean)
      p = pmixinvgauss(x, para(1), para(2), para(3), scale)
    case (DIST_BIRNBAUM_SAUNDERS)
      scale = 1.0_dp
      if (force_mean) scale = 1.0_dp / (1.0_dp + 0.5_dp * para(1)**2)
      p = pbirnbaum_saunders(x, para(1), scale)
    case default
      p = nan_value()
    end select
  end function distribution_cdf

  function distribution_quantile(p, code, para, force_mean) result(x)
    real(dp), intent(in) :: p
    integer, intent(in) :: code
    real(dp), intent(in) :: para(:)
    logical, intent(in) :: force_mean
    real(dp) :: x, scale

    select case (code)
    case (DIST_EXPONENTIAL)
      x = -log(max(tiny_pos, 1.0_dp - p))
    case (DIST_WEIBULL)
      scale = weibull_scale(para(1), force_mean)
      x = scale * (-log(max(tiny_pos, 1.0_dp - p)))**(1.0_dp / para(1))
    case (DIST_BURR)
      scale = burr_theta(para(1), para(2), force_mean)
      x = qburr(p, scale, para(1), para(2))
    case (DIST_GENGAMMA)
      scale = gengamma_lambda(para(1), para(2), force_mean)
      x = qgengamma(p, para(2), para(1), scale)
    case (DIST_GENF)
      scale = genf_lambda(para(1), para(2), para(3), force_mean)
      x = qgenf(p, para(1), para(2), para(3), scale)
    case (DIST_QWEIBULL)
      scale = qweibull_b(para(1), para(2), force_mean)
      x = qqweibull(p, para(1), para(2), scale)
    case (DIST_MIXQWE)
      scale = mixqwe_b(para(1), para(2), para(3), para(4), force_mean)
      x = qmixqwe(p, para(1), para(2), para(3), para(4), scale)
    case (DIST_MIXQWW)
      scale = mixqww_b(para(1), para(2), para(3), para(4), &
                       para(5), force_mean)
      x = qmixqww(p, para(1), para(2), para(3), para(4), para(5), scale)
    case (DIST_MIXINVGAUSS)
      scale = mixinvgauss_phi(para(1), para(2), para(3), force_mean)
      x = qmixinvgauss(p, para(1), para(2), para(3), scale)
    case (DIST_BIRNBAUM_SAUNDERS)
      scale = 1.0_dp
      if (force_mean) scale = 1.0_dp / (1.0_dp + 0.5_dp * para(1)**2)
      x = qbirnbaum_saunders(p, para(1), scale)
    case default
      x = nan_value()
    end select
  end function distribution_quantile

  function sample_distribution(rng, code, para, force_mean) result(x)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: code
    real(dp), intent(in) :: para(:)
    logical, intent(in) :: force_mean
    real(dp) :: x, scale

    select case (code)
    case (DIST_EXPONENTIAL)
      x = -log(random_uniform(rng))
    case (DIST_WEIBULL)
      scale = weibull_scale(para(1), force_mean)
      x = scale * (-log(random_uniform(rng)))**(1.0_dp / para(1))
    case (DIST_BURR)
      scale = burr_theta(para(1), para(2), force_mean)
      x = rburr(rng, scale, para(1), para(2))
    case (DIST_GENGAMMA)
      scale = gengamma_lambda(para(1), para(2), force_mean)
      x = rgengamma(rng, para(2), para(1), scale)
    case (DIST_GENF)
      scale = genf_lambda(para(1), para(2), para(3), force_mean)
      x = rgenf(rng, para(1), para(2), para(3), scale)
    case (DIST_QWEIBULL)
      scale = qweibull_b(para(1), para(2), force_mean)
      x = rqweibull(rng, para(1), para(2), scale)
    case (DIST_MIXQWE)
      scale = mixqwe_b(para(1), para(2), para(3), para(4), force_mean)
      x = rmixqwe(rng, para(1), para(2), para(3), para(4), scale)
    case (DIST_MIXQWW)
      scale = mixqww_b(para(1), para(2), para(3), para(4), &
                       para(5), force_mean)
      x = rmixqww(rng, para(1), para(2), para(3), para(4), para(5), scale)
    case (DIST_MIXINVGAUSS)
      scale = mixinvgauss_phi(para(1), para(2), para(3), force_mean)
      x = rmixinvgauss(rng, para(1), para(2), para(3), scale)
    case (DIST_BIRNBAUM_SAUNDERS)
      scale = 1.0_dp
      if (force_mean) scale = 1.0_dp / (1.0_dp + 0.5_dp * para(1)**2)
      x = rbirnbaum_saunders(rng, para(1), scale)
    case default
      x = nan_value()
    end select
  end function sample_distribution

end module acdm_distributions
