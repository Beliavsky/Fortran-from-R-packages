! SPDX-License-Identifier: GPL-2.0-only
module propcis_bayesian
  !! Deterministic posterior integration and inversion, replacing simulation starts.
  use propcis_common
  use r_special, only: r_regularized_beta, r_log_beta, r_digamma, r_trigamma
  use rcppnumerical_integration_1d, only: integrate_1d, integration_result_t, gk61
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: diffci_bayes, rrci_bayes, orci_bayes, diffci_bayes_hpd, bayes_cdf
  integer, parameter, public :: bayes_difference = 1, bayes_risk_ratio = 2, bayes_odds_ratio = 3
  type :: integration_context
    real(dp) :: shapes(4), threshold, location, scale, log_beta2
    integer :: distribution
    logical :: upper
  end type integration_context
contains
  function diffci_bayes(x1, n1, x2, n2, a, b, c, d, conf_level) result(ci)
    integer, intent(in) :: x1 !! Successes in the first sample, between zero and n1.
    integer, intent(in) :: n1 !! Nonnegative first-sample trial count; zero permits prior-only inference.
    integer, intent(in) :: x2 !! Successes in the second sample, between zero and n2.
    integer, intent(in) :: n2 !! Nonnegative second-sample trial count; zero permits prior-only inference.
    real(dp), intent(in) :: a !! Positive finite first-sample beta prior success shape.
    real(dp), intent(in) :: b !! Positive finite first-sample beta prior failure shape.
    real(dp), intent(in) :: c !! Positive finite second-sample beta prior success shape.
    real(dp), intent(in) :: d !! Positive finite second-sample beta prior failure shape.
    real(dp), intent(in), optional :: conf_level !! Posterior interval probability in zero to one, default 0.95.
    type(ci_result) :: ci
    real(dp) :: shapes(4), level, alpha
    integer :: status1, status2
    level = 0.95_dp
    if (present(conf_level)) level = conf_level
    call prepare_posterior(x1, n1, x2, n2, a, b, c, d, level, shapes, ci)
    if (ci%status /= ci_ok) return
    alpha = 0.5_dp * (1.0_dp - level)
    ci%lower = posterior_quantile(bayes_difference, shapes, alpha, .false., status1)
    ci%upper = posterior_quantile(bayes_difference, shapes, alpha, .true., status2)
    ci%status = max(status1, status2)
  end function diffci_bayes

  function rrci_bayes(x1, n1, x2, n2, a, b, c, d, conf_level) result(ci)
    integer, intent(in) :: x1 !! Successes in the first sample, between zero and n1.
    integer, intent(in) :: n1 !! Nonnegative first-sample trial count; zero permits prior-only inference.
    integer, intent(in) :: x2 !! Successes in the second sample, between zero and n2.
    integer, intent(in) :: n2 !! Nonnegative second-sample trial count; zero permits prior-only inference.
    real(dp), intent(in) :: a !! Positive finite first-sample beta prior success shape.
    real(dp), intent(in) :: b !! Positive finite first-sample beta prior failure shape.
    real(dp), intent(in) :: c !! Positive finite second-sample beta prior success shape.
    real(dp), intent(in) :: d !! Positive finite second-sample beta prior failure shape.
    real(dp), intent(in), optional :: conf_level !! Posterior interval probability in zero to one, default 0.95.
    type(ci_result) :: ci
    real(dp) :: shapes(4), level, alpha
    integer :: status1, status2
    level = 0.95_dp
    if (present(conf_level)) level = conf_level
    call prepare_posterior(x1, n1, x2, n2, a, b, c, d, level, shapes, ci)
    if (ci%status /= ci_ok) return
    alpha = 0.5_dp * (1.0_dp - level)
    ci%lower = posterior_quantile(bayes_risk_ratio, shapes, alpha, .false., status1)
    ci%upper = posterior_quantile(bayes_risk_ratio, shapes, alpha, .true., status2)
    ci%status = max(status1, status2)
  end function rrci_bayes

  function orci_bayes(x1, n1, x2, n2, a, b, c, d, conf_level) result(ci)
    integer, intent(in) :: x1 !! Successes in the first sample, between zero and n1.
    integer, intent(in) :: n1 !! Nonnegative first-sample trial count; zero permits prior-only inference.
    integer, intent(in) :: x2 !! Successes in the second sample, between zero and n2.
    integer, intent(in) :: n2 !! Nonnegative second-sample trial count; zero permits prior-only inference.
    real(dp), intent(in) :: a !! Positive finite first-sample beta prior success shape.
    real(dp), intent(in) :: b !! Positive finite first-sample beta prior failure shape.
    real(dp), intent(in) :: c !! Positive finite second-sample beta prior success shape.
    real(dp), intent(in) :: d !! Positive finite second-sample beta prior failure shape.
    real(dp), intent(in), optional :: conf_level !! Posterior interval probability in zero to one, default 0.95.
    type(ci_result) :: ci
    real(dp) :: shapes(4), level, alpha
    integer :: status1, status2
    level = 0.95_dp
    if (present(conf_level)) level = conf_level
    call prepare_posterior(x1, n1, x2, n2, a, b, c, d, level, shapes, ci)
    if (ci%status /= ci_ok) return
    alpha = 0.5_dp * (1.0_dp - level)
    ci%lower = posterior_quantile(bayes_odds_ratio, shapes, alpha, .false., status1)
    ci%upper = posterior_quantile(bayes_odds_ratio, shapes, alpha, .true., status2)
    ci%status = max(status1, status2)
  end function orci_bayes

  function diffci_bayes_hpd(x1, n1, x2, n2, a, b, c, d, conf_level) result(ci)
    !! Shortest connected credible interval; not a disconnected HPD set.
    integer, intent(in) :: x1 !! Successes in the first sample, between zero and n1.
    integer, intent(in) :: n1 !! Nonnegative first-sample trial count; zero permits prior-only inference.
    integer, intent(in) :: x2 !! Successes in the second sample, between zero and n2.
    integer, intent(in) :: n2 !! Nonnegative second-sample trial count; zero permits prior-only inference.
    real(dp), intent(in) :: a !! Positive finite first-sample beta prior success shape.
    real(dp), intent(in) :: b !! Positive finite first-sample beta prior failure shape.
    real(dp), intent(in) :: c !! Positive finite second-sample beta prior success shape.
    real(dp), intent(in) :: d !! Positive finite second-sample beta prior failure shape.
    real(dp), intent(in), optional :: conf_level !! Posterior interval probability in zero to one, default 0.95.
    type(ci_result) :: ci
    real(dp) :: shapes(4), level, left, right, u, v, fu, fv, best_p, best_width, width, fraction
    real(dp), parameter :: golden = 0.6180339887498948482_dp
    integer :: iter, status, k, best_k
    level = 0.95_dp
    if (present(conf_level)) level = conf_level
    call prepare_posterior(x1, n1, x2, n2, a, b, c, d, level, shapes, ci)
    if (ci%status /= ci_ok) return
    if (all(shapes >= 1.0_dp) .and. shapes(1) == shapes(3) .and. shapes(2) == shapes(4)) then
      ci = diffci_bayes(x1, n1, x2, n2, a, b, c, d, level)
      return
    end if
    ! A coarse global scan also includes boundary intervals and both asymmetric modes.
    best_width = huge(1.0_dp)
    best_k = 0
    best_p = 0.0_dp
    do k = 0, 20
      fraction = real(k, dp) / 20.0_dp
      width = interval_width(fraction * (1.0_dp - level))
      if (width < best_width) then
        best_width = width
        best_k = k
        best_p = fraction * (1.0_dp - level)
      end if
    end do
    left = real(max(0, best_k - 1), dp) * (1.0_dp - level) / 20.0_dp
    right = real(min(20, best_k + 1), dp) * (1.0_dp - level) / 20.0_dp
    u = right - golden * (right - left)
    v = left + golden * (right - left)
    fu = interval_width(u)
    fv = interval_width(v)
    do iter = 1, 45
      if (fu < fv) then
        right = v
        v = u
        fv = fu
        u = right - golden * (right - left)
        fu = interval_width(u)
      else
        left = u
        u = v
        fu = fv
        v = left + golden * (right - left)
        fv = interval_width(v)
      end if
      if (right - left < 2.0e-8_dp * (1.0_dp - level)) exit
    end do
    if (fu < best_width) then
      best_width = fu
      best_p = u
    end if
    if (fv < best_width) best_p = v
    ci%lower = posterior_quantile(bayes_difference, shapes, best_p, .false., status)
    ci%status = max(ci%status, status)
    ci%upper = posterior_quantile(bayes_difference, shapes, max(0.0_dp, 1.0_dp - level - best_p), .true., status)
    ci%status = max(ci%status, status)
  contains
    function interval_width(p) result(width)
      real(dp), intent(in) :: p !! Lower-tail mass excluded by a candidate connected interval.
      real(dp) :: width, lower, upper
      integer :: code
      lower = posterior_quantile(bayes_difference, shapes, p, .false., code)
      ci%status = max(ci%status, code)
      upper = posterior_quantile(bayes_difference, shapes, max(0.0_dp, 1.0_dp - level - p), .true., code)
      ci%status = max(ci%status, code)
      width = upper - lower
    end function interval_width
  end function diffci_bayes_hpd

  pure subroutine prepare_posterior(x1, n1, x2, n2, a, b, c, d, level, shapes, ci)
    integer, intent(in) :: x1 !! First-sample successes, between zero and n1.
    integer, intent(in) :: n1 !! Nonnegative first-sample trials.
    integer, intent(in) :: x2 !! Second-sample successes, between zero and n2.
    integer, intent(in) :: n2 !! Nonnegative second-sample trials.
    real(dp), intent(in) :: a !! Positive first-sample prior success shape.
    real(dp), intent(in) :: b !! Positive first-sample prior failure shape.
    real(dp), intent(in) :: c !! Positive second-sample prior success shape.
    real(dp), intent(in) :: d !! Positive second-sample prior failure shape.
    real(dp), intent(in) :: level !! Posterior interval probability, strictly between zero and one.
    real(dp), intent(out) :: shapes(4) !! Posterior beta shapes [a+x1,b+n1-x1,c+x2,d+n2-x2].
    type(ci_result), intent(out) :: ci !! Initialized result with invalid status when validation fails.
    ci%conf_level = level
    ci%estimate = nan_value()
    ci%lower = nan_value()
    ci%upper = nan_value()
    ci%status = ci_invalid
    shapes = [a, b, c, d]
    if (.not. all(ieee_is_finite(shapes))) return
    if (any(shapes <= 0.0_dp) .or. .not. valid_level(level)) return
    if (n1 < 0 .or. n2 < 0 .or. x1 < 0 .or. x2 < 0 .or. x1 > n1 .or. x2 > n2) return
    shapes = shapes + real([x1, n1 - x1, x2, n2 - x2], dp)
    ci%status = ci_ok
    if (1.0_dp - level < 1.0e-8_dp) ci%status = ci_numerical
    if (minval(shapes) < 1.0e-6_dp .or. maxval(shapes) > 1.0e6_dp) ci%status = ci_numerical
  end subroutine prepare_posterior

  function bayes_cdf(distribution, threshold, shapes, status, upper_tail, tolerance) result(value)
    !! Posterior CDF or survival function, integrated on standardized beta-logit coordinates.
    integer, intent(in) :: distribution !! One of bayes_difference, bayes_risk_ratio, bayes_odds_ratio.
    real(dp), intent(in) :: threshold !! Finite difference or ratio; out-of-support thresholds are allowed.
    real(dp), intent(in) :: shapes(4) !! Positive finite posterior beta shapes, first pair then second pair.
    integer, intent(out) :: status !! ci_ok, ci_invalid, or ci_numerical if quadrature did not converge.
    logical, intent(in), optional :: upper_tail !! True returns P(parameter > threshold); false (default) returns CDF.
    real(dp), intent(in), optional :: tolerance !! Positive absolute quadrature tolerance, default 1e-10.
    real(dp) :: value, tol
    type(integration_context) :: context
    type(integration_result_t) :: integral
    logical :: upper
    status = ci_invalid
    value = nan_value()
    upper = .false.
    if (present(upper_tail)) upper = upper_tail
    tol = 1.0e-10_dp
    if (present(tolerance)) tol = tolerance
    if (.not. all(ieee_is_finite(shapes)) .or. .not. ieee_is_finite(threshold)) return
    if (.not. ieee_is_finite(tol)) return
    if (any(shapes <= 0.0_dp) .or. tol <= 0.0_dp .or. tol > 0.01_dp) return
    if (distribution < 1 .or. distribution > 3) return
    status = ci_ok
    if (distribution == bayes_difference) then
      if (threshold <= -1.0_dp) then
        value = merge(1.0_dp, 0.0_dp, upper)
        return
      else if (threshold >= 1.0_dp) then
        value = merge(0.0_dp, 1.0_dp, upper)
        return
      end if
    else if (threshold <= 0.0_dp) then
      value = merge(1.0_dp, 0.0_dp, upper)
      return
    end if
    if (minval(shapes) < 1.0e-6_dp .or. maxval(shapes) > 1.0e6_dp) then
      status = ci_numerical
      return
    end if
    context%shapes = shapes
    context%threshold = threshold
    context%distribution = distribution
    context%upper = upper
    context%location = r_digamma(shapes(3)) - r_digamma(shapes(4))
    context%scale = sqrt(r_trigamma(shapes(3)) + r_trigamma(shapes(4)))
    context%log_beta2 = r_log_beta(shapes(3), shapes(4))
    if (.not. all(ieee_is_finite([context%location, context%scale, context%log_beta2]))) then
      status = ci_numerical
      return
    end if
    call integrate_1d(posterior_integrand, -infinity(), infinity(), integral, &
      subdiv=1000, eps_abs=tol, eps_rel=1.0e-10_dp, rule=gk61, user_data=context)
    value = min(1.0_dp, max(0.0_dp, integral%value))
    if (integral%error_code /= 0 .or. .not. ieee_is_finite(integral%value)) status = ci_numerical
  end function bayes_cdf

  function posterior_integrand(coordinate, user_data) result(value)
    real(dp), intent(in) :: coordinate !! Standardized logit of the second posterior probability.
    class(*), intent(inout), optional :: user_data !! Integration context carrying posterior shapes and requested tail.
    real(dp) :: value, z, logp, logq, probability, complement, target, target_complement, weight, conditional
    value = nan_value()
    if (.not. present(user_data)) return
    select type (context => user_data)
    type is (integration_context)
      z = context%location + context%scale * coordinate
      call logit_components(z, logp, logq)
      weight = exp(context%shapes(3) * logp + context%shapes(4) * logq - &
        context%log_beta2 + log(context%scale))
      if (weight == 0.0_dp) then
        value = 0.0_dp
        return
      end if
      probability = exp(logp)
      complement = exp(logq)
      select case (context%distribution)
      case (bayes_difference)
        target = probability + context%threshold
        target_complement = complement - context%threshold
      case (bayes_risk_ratio)
        target = probability * context%threshold
        target_complement = 1.0_dp - target
      case (bayes_odds_ratio)
        call logit_components(z + log(context%threshold), logp, logq)
        target = exp(logp)
        target_complement = exp(logq)
      end select
      if (target <= 0.0_dp) then
        conditional = merge(1.0_dp, 0.0_dp, context%upper)
      else if (target_complement <= 0.0_dp) then
        conditional = merge(0.0_dp, 1.0_dp, context%upper)
      else if (context%upper) then
        conditional = r_regularized_beta(min(1.0_dp, target_complement), context%shapes(2), context%shapes(1))
      else
        conditional = r_regularized_beta(min(1.0_dp, target), context%shapes(1), context%shapes(2))
      end if
      value = weight * conditional
    end select
  end function posterior_integrand

  pure subroutine logit_components(logit, logp, logq)
    real(dp), intent(in) :: logit !! Real log odds, without exponentiation overflow.
    real(dp), intent(out) :: logp !! Log of the probability represented by the log odds.
    real(dp), intent(out) :: logq !! Log of its complement, preserved in extreme upper tails.
    real(dp) :: normalizer
    normalizer = log(1.0_dp + exp(-abs(logit)))
    if (logit >= 0.0_dp) then
      logp = -normalizer
      logq = -logit - normalizer
    else
      logp = logit - normalizer
      logq = -normalizer
    end if
  end subroutine logit_components

  function posterior_quantile(distribution, shapes, tail_probability, upper, status) result(value)
    integer, intent(in) :: distribution !! Posterior transformation: difference, risk ratio, or odds ratio.
    real(dp), intent(in) :: shapes(4) !! Positive posterior beta shapes for both samples.
    real(dp), intent(in) :: tail_probability !! Probability in the requested tail, between zero and one.
    logical, intent(in) :: upper !! True inverts survival probability; false inverts CDF probability.
    integer, intent(out) :: status !! ci_ok or ci_numerical if integration or bracketing fails.
    real(dp) :: value, lo, hi, mid, trial, probability, tol
    integer :: code, iteration
    status = ci_ok
    if (tail_probability <= 0.0_dp) then
      if (distribution == bayes_difference) then
        value = merge(1.0_dp, -1.0_dp, upper)
      else
        value = merge(infinity(), 0.0_dp, upper)
      end if
      return
    end if
    tol = min(1.0e-11_dp, max(tiny(1.0_dp), tail_probability * 1.0e-5_dp))
    if (distribution == bayes_difference) then
      lo = -1.0_dp
      hi = 1.0_dp
    else
      lo = -1.0_dp
      hi = 1.0_dp
      do iteration = 1, 12
        probability = bayes_cdf(distribution, exp(lo), shapes, code, upper, tol)
        status = max(status, code)
        if ((.not. upper .and. probability <= tail_probability) .or. &
            (upper .and. probability >= tail_probability)) exit
        lo = 2.0_dp * lo
        if (lo < -700.0_dp) exit
      end do
      do iteration = 1, 12
        probability = bayes_cdf(distribution, exp(hi), shapes, code, upper, tol)
        status = max(status, code)
        if ((.not. upper .and. probability >= tail_probability) .or. &
            (upper .and. probability <= tail_probability)) exit
        hi = 2.0_dp * hi
        if (hi > 700.0_dp) exit
      end do
      if (lo < -700.0_dp .or. hi > 700.0_dp) then
        value = nan_value()
        status = ci_numerical
        return
      end if
    end if
    do iteration = 1, 65
      mid = 0.5_dp * (lo + hi)
      trial = mid
      if (distribution /= bayes_difference) trial = exp(mid)
      probability = bayes_cdf(distribution, trial, shapes, code, upper, tol)
      status = max(status, code)
      if (.not. ieee_is_finite(probability)) then
        value = nan_value()
        status = ci_numerical
        return
      end if
      if ((.not. upper .and. probability < tail_probability) .or. &
          (upper .and. probability > tail_probability)) then
        lo = mid
      else
        hi = mid
      end if
      if (hi - lo < 2.0e-10_dp) exit
    end do
    value = 0.5_dp * (lo + hi)
    if (distribution /= bayes_difference) value = exp(value)
  end function posterior_quantile
end module propcis_bayesian
