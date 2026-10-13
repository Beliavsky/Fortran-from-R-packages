! SPDX-License-Identifier: GPL-2.0-only
module propcis_frequentist
  !! Frequentist intervals translated from PropCIs 0.3-0 by Ralph Scherer.
  use propcis_common
  use r_distributions, only: r_qnorm
  use extra_distr_math, only: beta_quantile, binom_cdf, binom_quantile
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
  implicit none
  private
  public :: acceptbin, add4ci, addz2ci, blakerci, exactci, midpci, scoreci
  public :: diffpropci_mp, diffpropci_wald_mp, scoreci_mp, oddsratioci_mp
  public :: wald2ci, diffscoreci, z2stat, orscoreci, limit, riskscoreci
contains
  pure elemental function add4ci(x,n,conf_level) result(ci)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: z, trials, phat, bound
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    z = -r_qnorm(0.5_dp * (1.0_dp - conf_level))
    trials = real(n, dp) + 4.0_dp
    phat = (real(x, dp) + 0.5_dp * (4.0_dp)) / trials
    bound = z * sqrt(phat * (1.0_dp - phat) / trials)
    ci%lower = max(0.0_dp, phat - bound)
    ci%upper = min(1.0_dp, phat + bound)
    ci%estimate = phat
  end function add4ci

  pure elemental function addz2ci(x,n,conf_level) result(ci)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: z, trials, phat, bound
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    z = -r_qnorm(0.5_dp * (1.0_dp - conf_level))
    trials = real(n, dp) + z * z
    phat = (real(x, dp) + 0.5_dp * (z * z)) / trials
    bound = z * sqrt(phat * (1.0_dp - phat) / trials)
    ci%lower = max(0.0_dp, phat - bound)
    ci%upper = min(1.0_dp, phat + bound)
    ci%estimate = phat
  end function addz2ci

  pure elemental function scoreci(x,n,conf_level) result(ci)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: z, rn, phat, bound, midpoint
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    z = -r_qnorm(0.5_dp * (1.0_dp - conf_level))
    rn = real(n, dp)
    phat = real(x, dp) / rn
    bound = z * sqrt((phat * (1.0_dp - phat) + z * z / (4.0_dp * rn)) / rn) / (1.0_dp + z * z / rn)
    midpoint = (phat + z * z / (2.0_dp * rn)) / (1.0_dp + z * z / rn)
    ci%lower = anint((midpoint - bound) * 10000.0_dp) / 10000.0_dp
    ci%upper = anint((midpoint + bound) * 10000.0_dp) / 10000.0_dp
  end function scoreci

  function exactci(x,n,conf_level) result(ci)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: alpha
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    alpha = 1.0_dp - conf_level
    if (x == 0) then
      ci%lower = 0.0_dp
      ci%upper = 1.0_dp - (alpha / 2.0_dp)**(1.0_dp / real(n, dp))
    else if (x == n) then
      ci%lower = (alpha / 2.0_dp)**(1.0_dp / real(n, dp))
      ci%upper = 1.0_dp
    else
      ci%lower = beta_quantile(alpha / 2.0_dp, real(x, dp), real(n - x + 1, dp))
      ci%upper = beta_quantile(1.0_dp - alpha / 2.0_dp, real(x + 1, dp), real(n - x, dp))
      if (alpha < 1.0e-10_dp) ci%status = ci_numerical
    end if
  end function exactci

  pure elemental real(dp) function acceptbin(x, n, p) result(value)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: p !! Null binomial probability in the closed interval zero to one.
    real(dp) :: p1, p2, a1, a2
    value = nan_value()
    if (.not. valid_counts(x, n) .or. .not. ieee_is_finite(p)) return
    if (p < 0.0_dp .or. p > 1.0_dp) return
    p1 = 1.0_dp - binom_cdf(x - 1, n, p)
    p2 = binom_cdf(x, n, p)
    a1 = p1 + binom_cdf(binom_quantile(p1, n, p) - 1, n, p)
    a2 = p2 + 1.0_dp - binom_cdf(binom_quantile(1.0_dp - p2, n, p), n, p)
    value = min(a1, a2)
  end function acceptbin

  function blakerci(x, n, conf_level, tolerance) result(ci)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    real(dp), intent(in), optional :: tolerance !! Positive grid step, default 1e-5, at most one.
    type(ci_result) :: ci
    real(dp) :: step, alpha
    integer :: k, max_steps
    step = 1.0e-5_dp
    if (present(tolerance)) step = tolerance
    ci = exactci(x, n, conf_level)
    if (ci%status /= ci_ok) return
    if (.not. ieee_is_finite(step)) then
      ci%status = ci_invalid
    else if (step < 1.0e-9_dp .or. step > 1.0_dp) then
      ci%status = ci_invalid
    end if
    if (ci%status /= ci_ok) then
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    alpha = 1.0_dp - conf_level
    if (alpha < 1.0e-10_dp) ci%status = ci_numerical
    max_steps = ceiling(1.0_dp / step)
    if (x /= 0) then
      do k = 1, max_steps
        if (ci%lower + step >= 1.0_dp) exit
        if (acceptbin(x, n, ci%lower + step) >= alpha) exit
        ci%lower = ci%lower + step
      end do
    end if
    if (x /= n) then
      do k = 1, max_steps
        if (ci%upper - step <= 0.0_dp) exit
        if (acceptbin(x, n, ci%upper - step) >= alpha) exit
        ci%upper = ci%upper - step
      end do
    end if
  end function blakerci

  pure elemental function midpci(x,n,conf_level) result(ci)
    integer, intent(in) :: x !! Number of successes, between zero and n.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: alpha, pp, a2
    integer :: k
    logical :: got_lower, got_upper
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    alpha = 1.0_dp - conf_level
    ci%lower = 0.0_dp
    ci%upper = 1.0_dp
    if (x == 0) then
      ci%upper = 1.0_dp - alpha**(1.0_dp / real(n, dp))
    else if (x == n) then
      ci%lower = alpha**(1.0_dp / real(n, dp))
    else
      got_lower = .false.
      got_upper = .false.
      do k = 0, 1999
        pp = 0.0001_dp + real(k, dp) * 0.0005_dp
        a2 = 0.5_dp * (binom_cdf(x - 1, n, pp) + binom_cdf(x, n, pp))
        if (a2 > alpha / 2.0_dp) then
          ci%upper = pp
          got_upper = .true.
        end if
        if (a2 < 1.0_dp - alpha / 2.0_dp .and. .not. got_lower) then
          ci%lower = pp
          got_lower = .true.
        end if
      end do
      if (.not. (got_lower .and. got_upper)) ci%status = ci_numerical
      if (alpha < 1.0e-10_dp) ci%status = ci_numerical
    end if
  end function midpci

  pure elemental function diffpropci_mp(b,c,n,conf_level) result(ci)
    integer, intent(in) :: b !! First discordant-pair count, nonnegative.
    integer, intent(in) :: c !! Second discordant-pair count, nonnegative.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: z, rn, difference, sd
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_paired(b, c, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    rn = real(n, dp) + 2.0_dp
    z = -r_qnorm(0.5_dp * (1.0_dp - conf_level))
    difference = real(c - b, dp) / rn
    sd = sqrt(max(0.0_dp, real(b + c, dp) + 1.0_dp - real(c - b, dp)**2 / rn)) / rn
    ci%lower = difference - z * sd
    ci%upper = difference + z * sd
    ci%estimate = difference
    ci%lower = max(-1.0_dp, ci%lower)
    ci%upper = min(1.0_dp, ci%upper)
  end function diffpropci_mp

  pure elemental function diffpropci_wald_mp(b,c,n,conf_level) result(ci)
    integer, intent(in) :: b !! First discordant-pair count, nonnegative.
    integer, intent(in) :: c !! Second discordant-pair count, nonnegative.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: z, rn, difference, sd
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_paired(b, c, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    rn = real(n, dp)
    z = -r_qnorm(0.5_dp * (1.0_dp - conf_level))
    difference = real(c - b, dp) / rn
    sd = sqrt(max(0.0_dp, real(b + c, dp) - real(c - b, dp)**2 / rn)) / rn
    ci%lower = difference - z * sd
    ci%upper = difference + z * sd
    ci%estimate = difference
  end function diffpropci_wald_mp

  pure elemental function oddsratioci_mp(b,c,conf_level) result(ci)
    integer, intent(in) :: b !! First discordant-pair count, nonnegative.
    integer, intent(in) :: c !! Second discordant-pair count, nonnegative.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: z, aa, bb, cc, low, up
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (b >= 0 .and. c >= 0 .and. real(b, dp) + real(c, dp) > 0.0_dp .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    z = r_qnorm(0.5_dp * (1.0_dp - conf_level))**2
    aa = real(b, dp) + real(c, dp) + z
    bb = 2.0_dp * real(c, dp) + z
    cc = real(c, dp)**2 / (real(b, dp) + real(c, dp))
    low = (bb - sqrt(max(0.0_dp, bb * bb - 4.0_dp * aa * cc))) / (2.0_dp * aa)
    up = (bb + sqrt(max(0.0_dp, bb * bb - 4.0_dp * aa * cc))) / (2.0_dp * aa)
    ci%lower = low / (1.0_dp - low)
    if (c == 0) ci%lower = 0.0_dp
    if (b == 0 .or. up >= 1.0_dp) then
      ci%upper = infinity()
    else
      ci%upper = up / (1.0_dp - up)
    end if
  end function oddsratioci_mp

  pure elemental function wald2ci(x1, n1, x2, n2, conf_level, adjust) result(ci)
    integer, intent(in) :: x1 !! Successes in the first independent sample, between zero and n1.
    integer, intent(in) :: n1 !! Positive trial count of the first sample.
    integer, intent(in) :: x2 !! Successes in the second independent sample, between zero and n2.
    integer, intent(in) :: n2 !! Positive trial count of the second sample.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    character(len=*), intent(in), optional :: adjust !! "AC" (default) adds one success and one failure per sample; "Wald" adds none.
    type(ci_result) :: ci
    real(dp) :: nn1, nn2, p1, p2, added, bound
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x1, n1) .and. valid_counts(x2, n2) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    added = 1.0_dp
    if (present(adjust)) then
      select case (adjust)
      case ('AC')
        added = 1.0_dp
      case ('Wald')
        added = 0.0_dp
      case default
        ci%status = ci_invalid
        ci%lower = nan_value()
        ci%upper = nan_value()
        return
      end select
    end if
    nn1 = real(n1, dp) + 2.0_dp * added
    nn2 = real(n2, dp) + 2.0_dp * added
    p1 = (real(x1, dp) + added) / nn1
    p2 = (real(x2, dp) + added) / nn2
    bound = -r_qnorm(0.5_dp * (1.0_dp - conf_level)) * &
      sqrt(p1 * (1.0_dp - p1) / nn1 + p2 * (1.0_dp - p2) / nn2)
    ci%lower = p1 - p2 - bound
    ci%upper = p1 - p2 + bound
    ci%estimate = p1 - p2
  end function wald2ci

  pure elemental function scoreci_mp(b,c,n,conf_level) result(ci)
    integer, intent(in) :: b !! First discordant-pair count, nonnegative.
    integer, intent(in) :: c !! Second discordant-pair count, nonnegative.
    integer, intent(in) :: n !! Positive total number of trials or pairs.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: rn, z, proot, step, trial, pb, pc, q21, variance, score, direction
    integer :: side, iter
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_paired(b, c, n) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    rn = real(n, dp)
    z = -r_qnorm(0.5_dp * (1.0_dp - conf_level))
    do side = 1, 2
      direction = real(2 * side - 3, dp)
      if ((side == 1 .and. b == n) .or. (side == 2 .and. c == n)) then
        trial = direction
      else
        proot = real(c - b, dp) / rn
        step = 1.0_dp - direction * proot
        do iter = 1, 50
          step = step / 2.0_dp
          trial = proot + direction * step
          pb = -real(b + c, dp) + (2.0_dp * rn - real(c - b, dp)) * trial
          pc = -real(b, dp) * trial * (1.0_dp - trial)
          q21 = (sqrt(max(0.0_dp, pb * pb - 8.0_dp * rn * pc)) - pb) / (4.0_dp * rn)
          variance = rn * (2.0_dp * q21 + trial * (1.0_dp - trial))
          if (variance <= 0.0_dp) then
            score = infinity()
          else
            score = (real(c - b, dp) - rn * trial) / sqrt(variance)
          end if
          if (abs(score) < z) proot = trial
          if (step < 1.0e-7_dp .or. abs(z - score) < 1.0e-6_dp) exit
        end do
      end if
      if (side == 1) ci%lower = trial
      if (side == 2) ci%upper = trial
    end do
  end function scoreci_mp

  pure elemental real(dp) function z2stat(p1x, nx, p1y, ny, dif) result(score)
    !! Miettinen-Nurminen score with a constrained binomial likelihood maximum.
    real(dp), intent(in) :: p1x !! Observed success proportion of the first sample, in zero to one.
    integer, intent(in) :: nx !! Positive size of the first sample.
    real(dp), intent(in) :: p1y !! Observed success proportion of the second sample, in zero to one.
    integer, intent(in) :: ny !! Positive size of the second sample.
    real(dp), intent(in) :: dif !! Null first-minus-second proportion difference, in minus one to one.
    real(dp) :: lo, hi, p1, p2, derivative, variance, total, delta
    integer :: iter
    score = nan_value()
    if (.not. all(ieee_is_finite([p1x, p1y, dif]))) return
    if (nx <= 0 .or. ny <= 0) return
    if (min(p1x, p1y) < 0.0_dp .or. max(p1x, p1y) > 1.0_dp .or. abs(dif) > 1.0_dp) return
    delta = p1x - p1y - dif
    if (delta == 0.0_dp) then
      score = 0.0_dp
      return
    end if
    lo = max(0.0_dp, dif)
    hi = min(1.0_dp, 1.0_dp + dif)
    do iter = 1, 120
      p1 = 0.5_dp * (lo + hi)
      if (p1 <= lo .or. p1 >= hi) exit
      p2 = p1 - dif
      derivative = binomial_derivative(p1x, nx, p1) + binomial_derivative(p1y, ny, p2)
      if (derivative > 0.0_dp) then
        lo = p1
      else
        hi = p1
      end if
    end do
    p1 = 0.5_dp * (lo + hi)
    p2 = min(1.0_dp, max(0.0_dp, p1 - dif))
    total = real(nx, dp) + real(ny, dp)
    variance = (p1 * (1.0_dp - p1) / real(nx, dp) + p2 * (1.0_dp - p2) / real(ny, dp)) * &
      total / (total - 1.0_dp)
    if (variance <= 0.0_dp) then
      score = infinity()
    else
      score = delta * delta / variance
    end if
  end function z2stat

  pure elemental real(dp) function binomial_derivative(observed, n, probability) result(value)
    real(dp), intent(in) :: observed !! Observed success proportion in zero to one.
    integer, intent(in) :: n !! Positive binomial trial count.
    real(dp), intent(in) :: probability !! Candidate success probability in zero to one.
    value = 0.0_dp
    if (observed > 0.0_dp) then
      if (probability <= 0.0_dp) then
        value = huge(1.0_dp) / 4.0_dp
        return
      end if
      value = real(n, dp) * observed / probability
    end if
    if (observed < 1.0_dp) then
      if (probability >= 1.0_dp) then
        value = -huge(1.0_dp) / 4.0_dp
        return
      end if
      value = value - real(n, dp) * (1.0_dp - observed) / (1.0_dp - probability)
    end if
  end function binomial_derivative

  pure elemental function diffscoreci(x1, n1, x2, n2, conf_level) result(ci)
    integer, intent(in) :: x1 !! Successes in the first sample, between zero and n1.
    integer, intent(in) :: n1 !! Positive trial count of the first sample.
    integer, intent(in) :: x2 !! Successes in the second sample, between zero and n2.
    integer, intent(in) :: n2 !! Positive trial count of the second sample.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: p1, p2, z, proot, step, trial, score, direction
    integer :: side, iter
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x1, n1) .and. valid_counts(x2, n2) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    p1 = real(x1, dp) / real(n1, dp)
    p2 = real(x2, dp) / real(n2, dp)
    z = r_qnorm(0.5_dp * (1.0_dp - conf_level))**2
    do side = 1, 2
      direction = real(2 * side - 3, dp)
      proot = p1 - p2
      step = 1.0_dp - direction * proot
      do iter = 1, 50
        step = 0.5_dp * step
        trial = proot + direction * step
        score = z2stat(p1, n1, p2, n2, trial)
        if (score < z) proot = trial
        if (step < 1.0e-7_dp .or. abs(z - score) < 1.0e-6_dp) exit
      end do
      if (side == 1) ci%lower = trial
      if (side == 2) ci%upper = trial
    end do
  end function diffscoreci

  pure elemental real(dp) function limit(x, nx, y, ny, conflev, lim, t) result(value)
    !! Original multiplicative odds-ratio grid, with a stable quadratic root.
    integer, intent(in) :: x !! First-sample success count, between zero and nx.
    integer, intent(in) :: nx !! Positive first-sample trial count.
    integer, intent(in) :: y !! Second-sample success count, between zero and ny.
    integer, intent(in) :: ny !! Positive second-sample trial count.
    real(dp), intent(in) :: conflev !! Confidence probability, strictly between zero and one.
    real(dp), intent(in) :: lim !! Positive finite initial odds-ratio candidate.
    integer, intent(in) :: t !! Search direction: zero decreases, one increases the odds ratio.
    real(dp) :: theta, score, z, aa, bb, cc, p1, p2, v1, v2, delta, total
    integer :: iter
    value = nan_value()
    if (.not. (valid_counts(x, nx) .and. valid_counts(y, ny) .and. valid_level(conflev))) return
    if (.not. ieee_is_finite(lim)) return
    if (lim <= 0.0_dp .or. t < 0 .or. t > 1) return
    if ((x == 0 .and. y == 0) .or. (x == nx .and. y == ny)) return
    z = r_qnorm(0.5_dp * (1.0_dp - conflev))**2
    theta = lim
    total = real(nx, dp) + real(ny, dp)
    do iter = 1, 2000000
      aa = real(ny, dp) * (theta - 1.0_dp)
      bb = real(nx, dp) * theta + real(ny, dp) - (real(x, dp) + real(y, dp)) * (theta - 1.0_dp)
      cc = -(real(x, dp) + real(y, dp))
      if (bb >= 0.0_dp) then
        p2 = -2.0_dp * cc / (bb + sqrt(max(0.0_dp, bb * bb - 4.0_dp * aa * cc)))
      else
        p2 = (-bb + sqrt(max(0.0_dp, bb * bb - 4.0_dp * aa * cc))) / (2.0_dp * aa)
      end if
      p1 = p2 * theta / (1.0_dp + p2 * (theta - 1.0_dp))
      v1 = real(nx, dp) * p1 * (1.0_dp - p1)
      v2 = real(ny, dp) * p2 * (1.0_dp - p2)
      delta = real(x, dp) - real(nx, dp) * p1
      if (min(v1, v2) <= 0.0_dp) then
        score = infinity()
      else
        score = delta * delta * (1.0_dp / v1 + 1.0_dp / v2) * (total - 1.0_dp) / total
      end if
      if (score >= z) then
        value = theta
        return
      end if
      if (t == 0) then
        theta = theta / 1.001_dp
      else
        theta = theta * 1.001_dp
      end if
      if (.not. ieee_is_finite(theta) .or. theta <= tiny(1.0_dp)) exit
    end do
  end function limit

  pure elemental function orscoreci(x1, n1, x2, n2, conf_level) result(ci)
    integer, intent(in) :: x1 !! Successes in the first sample, between zero and n1.
    integer, intent(in) :: n1 !! Positive first-sample trial count.
    integer, intent(in) :: x2 !! Successes in the second sample, between zero and n2.
    integer, intent(in) :: n2 !! Positive second-sample trial count.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: theta
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x1, n1) .and. valid_counts(x2, n2) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    if ((x1 == 0 .and. x2 == 0) .or. (x1 == n1 .and. x2 == n2)) then
      ci%lower = 0.0_dp
      ci%upper = infinity()
    else if (x1 == 0 .or. x2 == n2) then
      ci%lower = 0.0_dp
      ci%upper = limit(x1, n1, x2, n2, conf_level, 0.01_dp / real(n2, dp), 1)
    else if (x1 == n1 .or. x2 == 0) then
      ci%lower = limit(x1, n1, x2, n2, conf_level, 100.0_dp * real(n1, dp), 0)
      ci%upper = infinity()
    else
      theta = (real(x1, dp) / real(n1 - x1, dp)) / (real(x2, dp) / real(n2 - x2, dp))
      ci%lower = limit(x1, n1, x2, n2, conf_level, theta / 1.1_dp, 0)
      ci%upper = limit(x1, n1, x2, n2, conf_level, theta * 1.1_dp, 1)
    end if
    if (ieee_is_nan(ci%lower) .or. ieee_is_nan(ci%upper)) ci%status = ci_numerical
    if (ci%lower < 0.0_dp .or. ci%upper < ci%lower) ci%status = ci_numerical
  end function orscoreci

  pure elemental function riskscoreci(x1, n1, x2, n2, conf_level) result(ci)
    !! Inverts the same constrained-binomial score using bracketed roots.
    integer, intent(in) :: x1 !! Successes in the numerator sample, between zero and n1.
    integer, intent(in) :: n1 !! Positive trial count of the numerator sample.
    integer, intent(in) :: x2 !! Successes in the denominator sample, between zero and n2.
    integer, intent(in) :: n2 !! Positive trial count of the denominator sample.
    real(dp), intent(in) :: conf_level !! Confidence probability, strictly between zero and one.
    type(ci_result) :: ci
    real(dp) :: observed, lo, hi, trial, z
    integer :: iter
    ci%conf_level = conf_level
    ci%estimate = nan_value()
    if (.not. (valid_counts(x1, n1) .and. valid_counts(x2, n2) .and. valid_level(conf_level))) then
      ci%status = ci_invalid
      ci%lower = nan_value()
      ci%upper = nan_value()
      return
    end if
    z = r_qnorm(0.5_dp * (1.0_dp - conf_level))**2
    ci%lower = 0.0_dp
    ci%upper = infinity()
    if (x1 == 0 .and. x2 == 0) return
    observed = 1.0_dp
    if (x2 > 0) observed = (real(x1, dp) / real(n1, dp)) / (real(x2, dp) / real(n2, dp))
    if (x1 > 0) then
      lo = 0.0_dp
      hi = observed
      do iter = 1, 1020
        if (risk_score(x1, n1, x2, n2, hi) < z) exit
        hi = hi * 2.0_dp
        if (hi > sqrt(huge(1.0_dp)) / 100.0_dp) exit
      end do
      do iter = 1, 180
        trial = 0.5_dp * (lo + hi)
        if (risk_score(x1, n1, x2, n2, trial) > z) then
          lo = trial
        else
          hi = trial
        end if
        if (hi - lo <= 1.0e-12_dp * max(1.0_dp, trial)) exit
      end do
      ci%lower = 0.5_dp * (lo + hi)
    end if
    if (x2 > 0) then
      lo = max(observed, ci%lower)
      hi = max(1.0_dp, 2.0_dp * lo)
      do iter = 1, 1020
        if (risk_score(x1, n1, x2, n2, hi) > z) exit
        hi = hi * 2.0_dp
        if (hi > sqrt(huge(1.0_dp)) / 100.0_dp) exit
      end do
      do iter = 1, 180
        trial = 0.5_dp * (lo + hi)
        if (risk_score(x1, n1, x2, n2, trial) > z) then
          hi = trial
        else
          lo = trial
        end if
        if (hi - lo <= 1.0e-12_dp * max(1.0_dp, trial)) exit
      end do
      ci%upper = 0.5_dp * (lo + hi)
    end if
  end function riskscoreci

  pure elemental real(dp) function risk_score(x1, n1, x2, n2, ratio) result(score)
    integer, intent(in) :: x1 !! Successes in the numerator sample.
    integer, intent(in) :: n1 !! Positive trial count of the numerator sample.
    integer, intent(in) :: x2 !! Successes in the denominator sample.
    integer, intent(in) :: n2 !! Positive trial count of the denominator sample.
    real(dp), intent(in) :: ratio !! Strictly positive candidate risk ratio.
    real(dp) :: aa, bb, cc, p1, p2, variance, delta
    aa = (real(n1, dp) + real(n2, dp)) * ratio
    bb = -(real(x1, dp) + real(n2, dp) + ratio * (real(x2, dp) + real(n1, dp)))
    cc = real(x1, dp) + real(x2, dp)
    p2 = 2.0_dp * cc / (-bb + sqrt(max(0.0_dp, bb * bb - 4.0_dp * aa * cc)))
    p2 = min(1.0_dp, max(0.0_dp, p2))
    p1 = min(1.0_dp, max(0.0_dp, ratio * p2))
    variance = p1 * (1.0_dp - p1) / real(n1, dp) + ratio * ratio * p2 * (1.0_dp - p2) / real(n2, dp)
    delta = real(x1, dp) / real(n1, dp) - ratio * real(x2, dp) / real(n2, dp)
    if (delta == 0.0_dp) then
      score = 0.0_dp
    else if (variance <= 0.0_dp) then
      score = infinity()
    else
      score = delta * delta / variance
    end if
  end function risk_score
end module propcis_frequentist
