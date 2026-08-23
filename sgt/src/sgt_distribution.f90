! SPDX-License-Identifier: GPL-3.0-or-later
module sgt_distribution
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan, &
    ieee_positive_inf, ieee_negative_inf
  use sgt_kinds, only : dp
  use sgt_special, only : reg_incomplete_beta, beta_quantile, reg_lower_gamma, gamma_quantile
  implicit none
  private
  public :: dsgt, psgt, qsgt, rsgt
  public :: sgt_logpdf, sgt_pdf, sgt_cdf, sgt_quantile
  public :: sgt_mean_shift, sgt_variance_scale
contains
  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp

  pure real(dp) function pos_inf() result(x)
    x = ieee_value(0.0_dp, ieee_positive_inf)
  end function pos_inf

  pure real(dp) function neg_inf() result(x)
    x = ieee_value(0.0_dp, ieee_negative_inf)
  end function neg_inf

  pure logical function valid_params(mu, sigma, lambda, p, q, mean_cent, var_adj) result(ok)
    real(dp), intent(in) :: mu, sigma, lambda, p, q
    logical, intent(in) :: mean_cent, var_adj
    logical :: pq_mean, pq_var
    ok = ieee_is_finite(mu) .and. ieee_is_finite(sigma) .and. sigma > 0.0_dp .and. &
      lambda > -1.0_dp .and. lambda < 1.0_dp .and. p > 0.0_dp .and. q > 0.0_dp
    if (.not. ok) return
    pq_mean = (.not. mean_cent) .or. (.not. ieee_is_finite(p)) .or. &
      (.not. ieee_is_finite(q)) .or. p * q > 1.0_dp
    pq_var = (.not. var_adj) .or. (.not. ieee_is_finite(p)) .or. &
      (.not. ieee_is_finite(q)) .or. p * q > 2.0_dp
    ok = pq_mean .and. pq_var
  end function valid_params

  pure real(dp) function log_beta(a, b) result(v)
    real(dp), intent(in) :: a, b
    v = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
  end function log_beta

  pure real(dp) function beta_ratio(a1, b1, a0, b0) result(v)
    real(dp), intent(in) :: a1, b1, a0, b0
    v = exp(log_beta(a1, b1) - log_beta(a0, b0))
  end function beta_ratio

  pure real(dp) function sgt_variance_scale(lambda, p, q) result(scale)
    real(dp), intent(in) :: lambda, p, q
    real(dp) :: b31, b21, term
    if (.not. ieee_is_finite(p)) then
      scale = sqrt(3.0_dp)
      return
    end if
    if (.not. ieee_is_finite(q)) then
      term = (acos(-1.0_dp) * (1.0_dp + 3.0_dp * lambda * lambda) * gamma(3.0_dp / p) - &
        16.0_dp ** (1.0_dp / p) * lambda * lambda * gamma(0.5_dp + 1.0_dp / p) ** 2 * &
        gamma(1.0_dp / p)) / (acos(-1.0_dp) * gamma(1.0_dp / p))
      scale = sqrt(term)
      return
    end if
    b31 = beta_ratio(3.0_dp / p, q - 2.0_dp / p, 1.0_dp / p, q)
    b21 = beta_ratio(2.0_dp / p, q - 1.0_dp / p, 1.0_dp / p, q)
    term = (3.0_dp * lambda * lambda + 1.0_dp) * b31 - 4.0_dp * lambda * lambda * b21 * b21
    scale = q ** (1.0_dp / p) * sqrt(term)
  end function sgt_variance_scale

  pure real(dp) function sgt_mean_shift(sigma_internal, lambda, p, q) result(shift)
    real(dp), intent(in) :: sigma_internal, lambda, p, q
    if (.not. ieee_is_finite(p)) then
      shift = 0.0_dp
    else if (.not. ieee_is_finite(q)) then
      shift = (2.0_dp ** (2.0_dp / p) * sigma_internal * lambda * gamma(0.5_dp + 1.0_dp / p)) / &
        sqrt(acos(-1.0_dp))
    else
      shift = 2.0_dp * sigma_internal * lambda * q ** (1.0_dp / p) * &
        beta_ratio(2.0_dp / p, q - 1.0_dp / p, 1.0_dp / p, q)
    end if
  end function sgt_mean_shift

  pure real(dp) function internal_sigma(sigma, lambda, p, q, var_adj, sigma_multiplier) result(s)
    real(dp), intent(in) :: sigma, lambda, p, q
    logical, intent(in) :: var_adj
    real(dp), intent(in) :: sigma_multiplier
    s = sigma * sigma_multiplier
    if (var_adj) s = s / sgt_variance_scale(lambda, p, q)
  end function internal_sigma

  pure elemental real(dp) function sgt_logpdf(x, mu, sigma, lambda, p, q, mean_cent, var_adj, &
      sigma_multiplier) result(logf)
    real(dp), intent(in) :: x, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj
    real(dp), intent(in), optional :: sigma_multiplier
    logical :: mc, va
    real(dp) :: sm, s, xx, z, side, lden
    mc = .true.
    va = .true.
    sm = 1.0_dp
    if (present(mean_cent)) mc = mean_cent
    if (present(var_adj)) va = var_adj
    if (present(sigma_multiplier)) sm = sigma_multiplier
    if (sm <= 0.0_dp .or. .not. valid_params(mu, sigma * sm, lambda, p, q, mc, va)) then
      logf = nan_dp()
      return
    end if
    s = internal_sigma(sigma, lambda, p, q, va, sm)
    xx = x
    if (mc) xx = xx + sgt_mean_shift(s, lambda, p, q)
    z = xx - mu
    if (z >= 0.0_dp) then
      side = 1.0_dp + lambda
    else
      side = 1.0_dp - lambda
    end if
    if (.not. ieee_is_finite(p)) then
      if (abs(z) <= s) then
        logf = -log(2.0_dp * s)
      else
        logf = neg_inf()
      end if
      return
    end if
    if (.not. ieee_is_finite(q)) then
      logf = log(p) - log(2.0_dp) - log(s) - log_gamma(1.0_dp / p) - &
        (abs(z) / (s * side)) ** p
      return
    end if
    lden = log(1.0_dp + abs(z) ** p / (q * s ** p * side ** p))
    logf = log(p) - log(2.0_dp) - log(s) - log(q) / p - log_beta(1.0_dp / p, q) - &
      (1.0_dp / p + q) * lden
  end function sgt_logpdf

  pure elemental real(dp) function sgt_pdf(x, mu, sigma, lambda, p, q, mean_cent, var_adj, &
      sigma_multiplier) result(f)
    real(dp), intent(in) :: x, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj
    real(dp), intent(in), optional :: sigma_multiplier
    real(dp) :: lf
    lf = sgt_logpdf(x, mu, sigma, lambda, p, q, mean_cent, var_adj, sigma_multiplier)
    f = exp(lf)
  end function sgt_pdf

  pure elemental real(dp) function sgt_cdf(quant, mu, sigma, lambda, p, q, mean_cent, var_adj, &
      lower_tail, log_p, sigma_multiplier) result(out)
    real(dp), intent(in) :: quant, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj, lower_tail, log_p
    real(dp), intent(in), optional :: sigma_multiplier
    logical :: mc, va, lt, lp, flip
    real(dp) :: sm, s, z, lam, arg, pg
    mc = .true.
    va = .true.
    lt = .true.
    lp = .false.
    sm = 1.0_dp
    if (present(mean_cent)) mc = mean_cent
    if (present(var_adj)) va = var_adj
    if (present(lower_tail)) lt = lower_tail
    if (present(log_p)) lp = log_p
    if (present(sigma_multiplier)) sm = sigma_multiplier
    if (sm <= 0.0_dp .or. .not. valid_params(mu, sigma * sm, lambda, p, q, mc, va)) then
      out = nan_dp()
      return
    end if
    if (.not. ieee_is_finite(quant)) then
      if (quant < 0.0_dp) then
        out = 0.0_dp
      else
        out = 1.0_dp
      end if
      if (.not. lt) out = 1.0_dp - out
      if (lp) then
        if (out <= 0.0_dp) then
          out = neg_inf()
        else
          out = log(out)
        end if
      end if
      return
    end if
    s = internal_sigma(sigma, lambda, p, q, va, sm)
    z = quant
    if (mc) z = z + sgt_mean_shift(s, lambda, p, q)
    z = z - mu
    if (.not. ieee_is_finite(p)) then
      if (z <= -s) then
        out = 0.0_dp
      else if (z >= s) then
        out = 1.0_dp
      else
        out = (z + s) / (2.0_dp * s)
      end if
    else if (.not. ieee_is_finite(q)) then
      flip = z < 0.0_dp
      lam = lambda
      if (flip) then
        lam = -lam
        z = -z
      end if
      arg = (z / (s * (1.0_dp + lam))) ** p
      pg = reg_lower_gamma(1.0_dp / p, arg)
      out = (1.0_dp - lam) / 2.0_dp + (1.0_dp + lam) * pg / 2.0_dp
      if (flip) out = 1.0_dp - out
    else
      flip = z > 0.0_dp
      lam = lambda
      if (flip) then
        lam = -lam
        z = -z
      end if
      if (abs(z) <= tiny(1.0_dp)) then
        out = (1.0_dp - lam) / 2.0_dp
      else
        arg = 1.0_dp / (1.0_dp + q * (s * (1.0_dp - lam) / (-z)) ** p)
        out = (1.0_dp - lam) / 2.0_dp + (lam - 1.0_dp) * &
          reg_incomplete_beta(arg, 1.0_dp / p, q) / 2.0_dp
      end if
      if (flip) out = 1.0_dp - out
    end if
    out = max(0.0_dp, min(1.0_dp, out))
    if (.not. lt) out = 1.0_dp - out
    if (lp) then
      if (out <= 0.0_dp) then
        out = neg_inf()
      else
        out = log(out)
      end if
    end if
  end function sgt_cdf

  pure elemental real(dp) function sgt_quantile(prob, mu, sigma, lambda, p, q, mean_cent, var_adj, &
      lower_tail, log_p, sigma_multiplier) result(out)
    real(dp), intent(in) :: prob, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj, lower_tail, log_p
    real(dp), intent(in), optional :: sigma_multiplier
    logical :: mc, va, lt, lp, flip
    real(dp) :: sm, pr, s, lam, u, bq
    mc = .true.
    va = .true.
    lt = .true.
    lp = .false.
    sm = 1.0_dp
    if (present(mean_cent)) mc = mean_cent
    if (present(var_adj)) va = var_adj
    if (present(lower_tail)) lt = lower_tail
    if (present(log_p)) lp = log_p
    if (present(sigma_multiplier)) sm = sigma_multiplier
    pr = prob
    if (lp) pr = exp(pr)
    if (.not. lt) pr = 1.0_dp - pr
    if (pr < 0.0_dp .or. pr > 1.0_dp .or. sm <= 0.0_dp .or. &
        .not. valid_params(mu, sigma * sm, lambda, p, q, mc, va)) then
      out = nan_dp()
      return
    end if
    if (pr <= 0.0_dp) then
      out = neg_inf()
      if (.not. ieee_is_finite(p)) out = mu - internal_sigma(sigma, lambda, p, q, va, sm)
      if (mc .and. .not. ieee_is_finite(p)) out = out - sgt_mean_shift( &
        internal_sigma(sigma, lambda, p, q, va, sm), lambda, p, q)
      return
    end if
    if (pr >= 1.0_dp) then
      out = pos_inf()
      if (.not. ieee_is_finite(p)) out = mu + internal_sigma(sigma, lambda, p, q, va, sm)
      if (mc .and. .not. ieee_is_finite(p)) out = out - sgt_mean_shift( &
        internal_sigma(sigma, lambda, p, q, va, sm), lambda, p, q)
      return
    end if
    s = internal_sigma(sigma, lambda, p, q, va, sm)
    if (.not. ieee_is_finite(p)) then
      out = mu - s + 2.0_dp * s * pr
    else if (.not. ieee_is_finite(q)) then
      flip = pr < (1.0_dp - lambda) / 2.0_dp
      if (flip) pr = 1.0_dp - pr
      lam = lambda
      if (flip) lam = -lam
      u = 2.0_dp * pr / (1.0_dp + lam) + (lam - 1.0_dp) / (lam + 1.0_dp)
      u = max(0.0_dp, min(1.0_dp, u))
      out = s * (1.0_dp + lam) * gamma_quantile(u, 1.0_dp / p) ** (1.0_dp / p)
      if (flip) out = -out
      out = out + mu
    else
      flip = pr > (1.0_dp - lambda) / 2.0_dp
      if (flip) pr = 1.0_dp - pr
      lam = lambda
      if (flip) lam = -lam
      u = 1.0_dp - 2.0_dp * pr / (1.0_dp - lam)
      u = max(0.0_dp, min(1.0_dp, u))
      bq = beta_quantile(u, 1.0_dp / p, q)
      if (bq <= 0.0_dp) then
        out = 0.0_dp
      else
        out = s * (lam - 1.0_dp) * (1.0_dp / (q * bq) - 1.0_dp / q) ** (-1.0_dp / p)
      end if
      if (flip) out = -out
      out = out + mu
    end if
    if (mc) out = out - sgt_mean_shift(s, lambda, p, q)
  end function sgt_quantile

  pure elemental real(dp) function dsgt(x, mu, sigma, lambda, p, q, mean_cent, var_adj, log_value, &
      sigma_multiplier) result(out)
    real(dp), intent(in) :: x, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj, log_value
    real(dp), intent(in), optional :: sigma_multiplier
    logical :: lv
    lv = .false.
    if (present(log_value)) lv = log_value
    if (lv) then
      out = sgt_logpdf(x, mu, sigma, lambda, p, q, mean_cent, var_adj, sigma_multiplier)
    else
      out = sgt_pdf(x, mu, sigma, lambda, p, q, mean_cent, var_adj, sigma_multiplier)
    end if
  end function dsgt

  pure elemental real(dp) function psgt(quant, mu, sigma, lambda, p, q, mean_cent, var_adj, lower_tail, &
      log_p, sigma_multiplier) result(out)
    real(dp), intent(in) :: quant, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj, lower_tail, log_p
    real(dp), intent(in), optional :: sigma_multiplier
    out = sgt_cdf(quant, mu, sigma, lambda, p, q, mean_cent, var_adj, lower_tail, log_p, sigma_multiplier)
  end function psgt

  pure elemental real(dp) function qsgt(prob, mu, sigma, lambda, p, q, mean_cent, var_adj, lower_tail, &
      log_p, sigma_multiplier) result(out)
    real(dp), intent(in) :: prob, mu, sigma, lambda, p, q
    logical, intent(in), optional :: mean_cent, var_adj, lower_tail, log_p
    real(dp), intent(in), optional :: sigma_multiplier
    out = sgt_quantile(prob, mu, sigma, lambda, p, q, mean_cent, var_adj, lower_tail, log_p, sigma_multiplier)
  end function qsgt

  subroutine rsgt(x, mu, sigma, lambda, p, q, mean_cent, var_adj, sigma_multiplier)
    real(dp), intent(out) :: x(:)
    real(dp), intent(in), optional :: mu, sigma, lambda, p, q, sigma_multiplier
    logical, intent(in), optional :: mean_cent, var_adj
    real(dp) :: muv, sigv, lamv, pv, qv, sm, u
    integer :: i
    muv = 0.0_dp
    sigv = 1.0_dp
    lamv = 0.0_dp
    pv = 2.0_dp
    qv = pos_inf()
    sm = 1.0_dp
    if (present(mu)) muv = mu
    if (present(sigma)) sigv = sigma
    if (present(lambda)) lamv = lambda
    if (present(p)) pv = p
    if (present(q)) qv = q
    if (present(sigma_multiplier)) sm = sigma_multiplier
    do i = 1, size(x)
      call random_number(u)
      u = max(tiny(1.0_dp), min(1.0_dp - epsilon(1.0_dp), u))
      x(i) = sgt_quantile(u, muv, sigv, lamv, pv, qv, mean_cent, var_adj, sigma_multiplier=sm)
    end do
  end subroutine rsgt
end module sgt_distribution
