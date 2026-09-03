module rmpfr_probability
  use rmpfr_kinds, only: dp, i64
  use rmpfr_types
  use rmpfr_combinatorics, only: mpfr_choose
  implicit none
  private

  public :: mpfr_pnorm, mpfr_dnorm, mpfr_dt, mpfr_dpois, mpfr_dbinom, mpfr_dnbinom
  public :: mpfr_dgamma, mpfr_dchisq, mpfr_pgamma, mpfr_pbeta_integer
  public :: mpfr_log1mexp, mpfr_log1pexp

contains

  recursive function mpfr_pnorm(q, mean, sd, lower_tail, log_p) result(r)
    type(mpfr_real), intent(in) :: q !! Quantile at which the normal distribution function is evaluated.
    type(mpfr_real), intent(in), optional :: mean !! Normal mean; defaults to zero at the working precision.
    type(mpfr_real), intent(in), optional :: sd !! Positive normal standard deviation; defaults to one.
    logical, intent(in), optional :: lower_tail !! If true return P[X <= q]; otherwise return the upper tail.
    logical, intent(in), optional :: log_p !! If true return the logarithm of the requested probability.
    type(mpfr_real) :: r, mu, sigma, z, zpos, half, one, two, inv_sqrt2, epart, cutoff
    logical :: lower, give_log
    integer :: p

    p = mpfr_precision(q)
    if (present(mean)) p = max(p, mpfr_precision(mean))
    if (present(sd)) p = max(p, mpfr_precision(sd))
    mu = mpfr_zero(1, p)
    if (present(mean)) mu = mpfr_copy(mean, p)
    sigma = mpfr_from_integer(1_i64, p)
    if (present(sd)) sigma = mpfr_copy(sd, p)
    if (mpfr_sign(sigma) <= 0) then
      r = mpfr_nan(p)
      return
    end if
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    give_log = .false.
    if (present(log_p)) give_log = log_p

    z = (mpfr_copy(q, p) - mu) / sigma
    if (mpfr_sign(z) < 0) then
      zpos = -z
      r = mpfr_pnorm(zpos, lower_tail=.not. lower, log_p=give_log)
      return
    end if

    half = mpfr_from_string("0.5", p)
    one = mpfr_from_integer(1_i64, p)
    two = mpfr_from_integer(2_i64, p)
    inv_sqrt2 = mpfr_sqrt(half)
    epart = z * inv_sqrt2
    if (lower) then
      if (give_log) then
        cutoff = mpfr_from_string("0.67448975", p)
        if (z < cutoff) then
          r = mpfr_log1p(mpfr_erf(epart)) - mpfr_log(two)
        else
          r = mpfr_log1p(-(half * mpfr_erfc(epart)))
        end if
      else
        r = (one + mpfr_erf(epart)) / two
      end if
    else
      r = mpfr_erfc(epart) / two
      if (give_log) r = mpfr_log(r)
    end if
  end function mpfr_pnorm

  function mpfr_dnorm(x, mean, sd, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Point at which the normal density is evaluated.
    type(mpfr_real), intent(in), optional :: mean !! Normal mean; defaults to zero.
    type(mpfr_real), intent(in), optional :: sd !! Positive normal standard deviation; defaults to one.
    logical, intent(in), optional :: log_density !! If true return the log density.
    type(mpfr_real) :: r, mu, sigma, z, twopi, two
    logical :: give_log
    integer :: p

    p = mpfr_precision(x)
    if (present(mean)) p = max(p, mpfr_precision(mean))
    if (present(sd)) p = max(p, mpfr_precision(sd))
    mu = mpfr_zero(1, p)
    if (present(mean)) mu = mpfr_copy(mean, p)
    sigma = mpfr_from_integer(1_i64, p)
    if (present(sd)) sigma = mpfr_copy(sd, p)
    if (mpfr_sign(sigma) <= 0) then
      r = mpfr_nan(p)
      return
    end if
    give_log = .false.
    if (present(log_density)) give_log = log_density
    two = mpfr_from_integer(2_i64, p)
    twopi = two * mpfr_const_pi(p)
    z = (mpfr_copy(x, p) - mu) / sigma
    if (give_log) then
      r = -(mpfr_log(sigma) + (mpfr_log(twopi) + z * z) / two)
    else
      r = mpfr_exp(-(z * z) / two) / (sigma * mpfr_sqrt(twopi))
    end if
  end function mpfr_dnorm

  function mpfr_dt(x, df, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Point at which the central Student-t density is evaluated.
    type(mpfr_real), intent(in) :: df !! Positive degrees of freedom.
    logical, intent(in), optional :: log_density !! If true return the log density.
    type(mpfr_real) :: r, xx, nu, one, two, three, twopi, term
    logical :: give_log
    integer :: p

    p = max(mpfr_precision(x), mpfr_precision(df))
    xx = mpfr_copy(x, p)
    nu = mpfr_copy(df, p)
    if (mpfr_sign(nu) <= 0) then
      r = mpfr_nan(p)
      return
    end if
    give_log = .false.
    if (present(log_density)) give_log = log_density
    one = mpfr_from_integer(1_i64, p)
    two = mpfr_from_integer(2_i64, p)
    three = mpfr_from_integer(3_i64, p)
    twopi = two * mpfr_const_pi(p)
    term = one + xx * xx / nu
    if (give_log) then
      r = mpfr_log(nu / two) / two - mpfr_log((nu + one) / two) + &
          mpfr_lgamma((nu + three) / two) - mpfr_lgamma((nu + two) / two) - &
          (nu / two) * mpfr_log1p(xx * xx / nu) - mpfr_log(twopi * term) / two
    else
      r = mpfr_sqrt(nu / two) / ((nu + one) / two) * &
          mpfr_gamma((nu + three) / two) / mpfr_gamma((nu + two) / two) * &
          term ** (-mpfr_from_real(0.5_dp, p) * nu) / mpfr_sqrt(twopi * term)
    end if
  end function mpfr_dt

  function mpfr_dpois(x, lambda, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Nonnegative count, represented exactly when possible.
    type(mpfr_real), intent(in) :: lambda !! Nonnegative Poisson mean.
    logical, intent(in), optional :: log_density !! If true return the log probability.
    type(mpfr_real) :: r, xx, lam, one, log_r
    logical :: give_log
    integer :: p

    p = max(mpfr_precision(x), mpfr_precision(lambda))
    xx = mpfr_copy(x, p)
    lam = mpfr_copy(lambda, p)
    give_log = .false.
    if (present(log_density)) give_log = log_density
    if (mpfr_sign(xx) < 0 .or. mpfr_sign(lam) < 0 .or. .not. mpfr_is_integer(xx)) then
      if (give_log) then
        r = mpfr_inf(-1, p)
      else
        r = mpfr_zero(1, p)
      end if
      return
    end if
    if (mpfr_is_zero(lam)) then
      if (mpfr_is_zero(xx)) then
        if (give_log) then
          r = mpfr_zero(1, p)
        else
          r = mpfr_from_integer(1_i64, p)
        end if
      else
        if (give_log) then
          r = mpfr_inf(-1, p)
        else
          r = mpfr_zero(1, p)
        end if
      end if
      return
    end if
    one = mpfr_from_integer(1_i64, p)
    log_r = -lam + xx * mpfr_log(lam) - mpfr_lgamma(xx + one)
    if (give_log) then
      r = log_r
    else
      r = mpfr_exp(log_r)
    end if
  end function mpfr_dpois

  function mpfr_dbinom(x, size, prob, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Integer number of successes.
    type(mpfr_real), intent(in) :: size !! Nonnegative integer number of trials.
    type(mpfr_real), intent(in) :: prob !! Success probability in [0,1].
    logical, intent(in), optional :: log_density !! If true return the log probability.
    type(mpfr_real) :: r, xx, nn, pp, zero, one, log_r
    logical :: give_log
    integer :: p

    p = max(mpfr_precision(x), mpfr_precision(size), mpfr_precision(prob))
    xx = mpfr_copy(x, p)
    nn = mpfr_copy(size, p)
    pp = mpfr_copy(prob, p)
    zero = mpfr_zero(1, p)
    one = mpfr_from_integer(1_i64, p)
    give_log = .false.
    if (present(log_density)) give_log = log_density
    if (.not. mpfr_is_integer(xx) .or. .not. mpfr_is_integer(nn) .or. xx < zero .or. &
        nn < xx .or. pp < zero .or. pp > one) then
      if (give_log) then
        r = mpfr_inf(-1, p)
      else
        r = mpfr_zero(1, p)
      end if
      return
    end if
    if (mpfr_is_zero(pp)) then
      if (mpfr_is_zero(xx)) then
        r = merge_logical_one(give_log, p)
      else
        r = merge_logical_zero(give_log, p)
      end if
      return
    end if
    if (pp == one) then
      if (xx == nn) then
        r = merge_logical_one(give_log, p)
      else
        r = merge_logical_zero(give_log, p)
      end if
      return
    end if
    log_r = mpfr_lgamma(nn + one) - mpfr_lgamma(xx + one) - mpfr_lgamma(nn - xx + one) + &
            xx * mpfr_log(pp) + (nn - xx) * mpfr_log1p(-pp)
    if (give_log) then
      r = log_r
    else
      r = mpfr_exp(log_r)
    end if
  end function mpfr_dbinom

  function mpfr_dnbinom(x, size, prob, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Nonnegative integer failure count.
    type(mpfr_real), intent(in) :: size !! Positive negative-binomial size parameter.
    type(mpfr_real), intent(in) :: prob !! Success probability in (0,1].
    logical, intent(in), optional :: log_density !! If true return the log probability.
    type(mpfr_real) :: r, xx, ss, pp, one, log_r
    logical :: give_log
    integer :: p

    p = max(mpfr_precision(x), mpfr_precision(size), mpfr_precision(prob))
    xx = mpfr_copy(x, p)
    ss = mpfr_copy(size, p)
    pp = mpfr_copy(prob, p)
    one = mpfr_from_integer(1_i64, p)
    give_log = .false.
    if (present(log_density)) give_log = log_density
    if (.not. mpfr_is_integer(xx) .or. mpfr_sign(xx) < 0 .or. mpfr_sign(ss) <= 0 .or. &
        mpfr_sign(pp) <= 0 .or. pp > one) then
      r = mpfr_nan(p)
      return
    end if
    log_r = mpfr_lgamma(ss + xx) - mpfr_lgamma(ss) - mpfr_lgamma(xx + one) + &
            ss * mpfr_log(pp) + xx * mpfr_log1p(-pp)
    if (give_log) then
      r = log_r
    else
      r = mpfr_exp(log_r)
    end if
  end function mpfr_dnbinom

  function mpfr_dgamma(x, shape, scale, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Nonnegative point at which the gamma density is evaluated.
    type(mpfr_real), intent(in) :: shape !! Positive gamma shape parameter.
    type(mpfr_real), intent(in), optional :: scale !! Positive gamma scale; defaults to one.
    logical, intent(in), optional :: log_density !! If true return the log density.
    type(mpfr_real) :: r, xx, aa, ss, one, log_r
    logical :: give_log
    integer :: p

    p = max(mpfr_precision(x), mpfr_precision(shape))
    if (present(scale)) p = max(p, mpfr_precision(scale))
    xx = mpfr_copy(x, p)
    aa = mpfr_copy(shape, p)
    ss = mpfr_from_integer(1_i64, p)
    if (present(scale)) ss = mpfr_copy(scale, p)
    give_log = .false.
    if (present(log_density)) give_log = log_density
    if (mpfr_sign(xx) < 0 .or. mpfr_sign(aa) <= 0 .or. mpfr_sign(ss) <= 0) then
      r = mpfr_nan(p)
      return
    end if
    one = mpfr_from_integer(1_i64, p)
    if (mpfr_is_zero(xx)) then
      if (aa < one) then
        if (give_log) then
          r = mpfr_inf(1, p)
        else
          r = mpfr_inf(1, p)
        end if
        return
      else if (aa == one) then
        log_r = -mpfr_log(ss)
        if (give_log) then
          r = log_r
        else
          r = mpfr_exp(log_r)
        end if
        return
      else
        r = merge_logical_zero(give_log, p)
        return
      end if
    end if
    log_r = -aa * mpfr_log(ss) - mpfr_lgamma(aa) + (aa - one) * mpfr_log(xx) - xx / ss
    if (give_log) then
      r = log_r
    else
      r = mpfr_exp(log_r)
    end if
  end function mpfr_dgamma

  function mpfr_dchisq(x, df, log_density) result(r)
    type(mpfr_real), intent(in) :: x !! Nonnegative point at which the chi-square density is evaluated.
    type(mpfr_real), intent(in) :: df !! Positive chi-square degrees of freedom.
    logical, intent(in), optional :: log_density !! If true return the log density.
    type(mpfr_real) :: r, two
    integer :: p

    p = max(mpfr_precision(x), mpfr_precision(df))
    two = mpfr_from_integer(2_i64, p)
    r = mpfr_dgamma(mpfr_copy(x, p), mpfr_copy(df, p) / two, scale=two, log_density=log_density)
  end function mpfr_dchisq

  function mpfr_pgamma(q, shape, scale, lower_tail, log_p, rounding) result(r)
    type(mpfr_real), intent(in) :: q !! Point at which the gamma distribution function is evaluated.
    type(mpfr_real), intent(in) :: shape !! Positive gamma shape parameter.
    type(mpfr_real), intent(in), optional :: scale !! Positive gamma scale; defaults to one.
    logical, intent(in), optional :: lower_tail !! If true return the lower-tail probability.
    logical, intent(in), optional :: log_p !! If true return the logarithm of the probability.
    integer, intent(in), optional :: rounding !! MPFR rounding code used for incomplete gamma evaluation.
    type(mpfr_real) :: r, qq, aa, ss, upper_ratio, one
    logical :: lower, give_log
    integer :: p

    p = max(mpfr_precision(q), mpfr_precision(shape))
    if (present(scale)) p = max(p, mpfr_precision(scale))
    qq = mpfr_copy(q, p)
    aa = mpfr_copy(shape, p)
    ss = mpfr_from_integer(1_i64, p)
    if (present(scale)) ss = mpfr_copy(scale, p)
    one = mpfr_from_integer(1_i64, p)
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    give_log = .false.
    if (present(log_p)) give_log = log_p
    if (mpfr_sign(aa) <= 0 .or. mpfr_sign(ss) <= 0) then
      r = mpfr_nan(p)
      return
    end if
    if (mpfr_sign(qq) <= 0) then
      if (lower) then
        r = merge_logical_zero(give_log, p)
      else
        r = merge_logical_one(give_log, p)
      end if
      return
    end if
    upper_ratio = mpfr_igamma(aa, qq / ss, rounding) / mpfr_gamma(aa, rounding)
    if (lower) then
      if (give_log) then
        r = mpfr_log1p(-upper_ratio)
      else
        r = one - upper_ratio
      end if
    else
      if (give_log) then
        r = mpfr_log(upper_ratio)
      else
        r = upper_ratio
      end if
    end if
  end function mpfr_pgamma

  function mpfr_pbeta_integer(q, shape1, shape2, lower_tail, log_p) result(r)
    type(mpfr_real), intent(in) :: q !! Probability argument in the closed interval [0,1].
    integer, intent(in) :: shape1 !! Positive integer first beta shape parameter.
    integer, intent(in) :: shape2 !! Positive integer second beta shape parameter.
    logical, intent(in), optional :: lower_tail !! If true return P[X <= q]; otherwise return the upper tail.
    logical, intent(in), optional :: log_p !! If true return the logarithm of the requested probability.
    type(mpfr_real) :: r, qq, zero, one, n_mpfr, coefficient, term
    integer :: n, k, kmax, p
    logical :: lower, logarithmic

    if (shape1 <= 0 .or. shape2 <= 0) error stop "Rmpfr: beta shapes must be positive integers"
    p = mpfr_precision(q)
    qq = mpfr_copy(q, p)
    zero = mpfr_zero(1, p)
    one = mpfr_from_integer(1_i64, p)
    if (qq < zero .or. qq > one) then
      r = mpfr_nan(p)
      return
    end if
    lower = .true.
    if (present(lower_tail)) lower = lower_tail
    logarithmic = .false.
    if (present(log_p)) logarithmic = log_p
    n = shape1 + shape2 - 1
    n_mpfr = mpfr_from_integer(int(n, i64), p)
    r = mpfr_zero(1, p)
    if (lower) then
      kmax = shape2 - 1
      do k = 0, kmax
        coefficient = mpfr_choose(n_mpfr, k)
        term = coefficient * (one - qq) ** k * qq ** (n - k)
        r = r + term
      end do
    else
      kmax = shape1 - 1
      do k = 0, kmax
        coefficient = mpfr_choose(n_mpfr, k)
        term = coefficient * qq ** k * (one - qq) ** (n - k)
        r = r + term
      end do
    end if
    if (logarithmic) r = mpfr_log(r)
  end function mpfr_pbeta_integer

  function mpfr_log1mexp(a) result(r)
    type(mpfr_real), intent(in) :: a !! Nonnegative argument of log(1-exp(-a)).
    type(mpfr_real) :: r, cutoff
    integer :: p

    p = mpfr_precision(a)
    if (mpfr_sign(a) < 0) then
      r = mpfr_nan(p)
      return
    end if
    cutoff = mpfr_const_log2(p)
    if (a <= cutoff) then
      r = mpfr_log(-mpfr_expm1(-a))
    else
      r = mpfr_log1p(-mpfr_exp(-a))
    end if
  end function mpfr_log1mexp

  function mpfr_log1pexp(x) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument of the stable softplus log(1+exp(x)).
    type(mpfr_real) :: r, c0, c1, c2, one, ex
    integer :: p

    p = mpfr_precision(x)
    c0 = mpfr_from_string("-37", p)
    c1 = mpfr_from_string("18", p)
    c2 = mpfr_from_string("33.3", p)
    one = mpfr_from_integer(1_i64, p)
    if (x <= c0) then
      r = mpfr_exp(x)
    else if (x <= c1) then
      r = mpfr_log1p(mpfr_exp(x))
    else if (x <= c2) then
      ex = mpfr_exp(x)
      r = x + one / ex
    else
      r = x
    end if
  end function mpfr_log1pexp

  function merge_logical_one(log_scale, prec_bits) result(r)
    logical, intent(in) :: log_scale !! Whether a probability equal to one should be returned on log scale.
    integer, intent(in) :: prec_bits !! Binary precision for the returned value.
    type(mpfr_real) :: r

    if (log_scale) then
      r = mpfr_zero(1, prec_bits)
    else
      r = mpfr_from_integer(1_i64, prec_bits)
    end if
  end function merge_logical_one

  function merge_logical_zero(log_scale, prec_bits) result(r)
    logical, intent(in) :: log_scale !! Whether a probability equal to zero should be returned on log scale.
    integer, intent(in) :: prec_bits !! Binary precision for the returned value.
    type(mpfr_real) :: r

    if (log_scale) then
      r = mpfr_inf(-1, prec_bits)
    else
      r = mpfr_zero(1, prec_bits)
    end if
  end function merge_logical_zero

end module rmpfr_probability
