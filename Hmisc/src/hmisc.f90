module hmisc
  use, intrinsic :: iso_fortran_env, only : real64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
  implicit none
  private

  integer, parameter, public :: dp = real64

  public :: gini_md, trap_rule, samplesize_bin
  public :: dual_sd, bpower, bsamsize, ballocation
  public :: weibull2_fit, weibull_survival, gompertz2_fit, gompertz_survival
  public :: lognorm2_fit, lognorm_survival, improve_prob
  public :: weighted_mean, weighted_variance, weighted_quantile, weighted_rank
  public :: weighted_table, somers2
  public :: rcorr, hoeffd, rcorr_cens, rcorrp_cens
  public :: cut_gn, which_closest, which_close_pw, pseudomedian
  public :: rcspline_eval, popower, posamsize, pomodm
  public :: groupn, spearman_rho, stepfun_eval, xy_sort_no_dup_no_na
  public :: n_coincident, matxv
  public :: cpower, ciapower, wtd_ecdf, smean_sd, smean_sdl, smedian_hilow
  public :: approx_extrap, james_stein
  public :: ftupwr, ftuss, ecdf_steps, lag_numeric, nomiss, fillin, cumcategory
  public :: deff, r2_measures, binconf_wilson, binconf_asymptotic
  public :: invert_tabulated_linear, score_binary_max, score_binary_sum
  public :: mhgr, lrcum, num_denom_setup, rmultinom, qrxcenter
  public :: solvet, solvet_inverse, pc1, smean_cl_normal, which_close_k
  public :: lm_fit_qr_bare, abs_error_pred, rcorrcens_summary
  public :: logrank, spearman2_numeric, equal_bins, jitter2_numeric
  public :: hdquantile, cat_test_chisq, con_test_kw, t_test_cluster
  public :: km_quick, bystats_mean, bystats2_mean, all_digits, all_is_numeric
  public :: find_matches, seq_freq, pair_up_diff
  public :: smean_cl_boot, smearing_est_tabulated
  public :: bpower_sim_counts, spower_simulated, inverse_function_all
  public :: bezier_curve, cut2_explicit, sim_po_cuts
  public :: rcspline_function, rcspline_restate_coefficients
  public :: spearman_test, chi_square_codes, combine_levels_codes, princmp_regular
  public :: bootkm_resampled, match_cases_core, largest_empty_rexhaustive
  public :: gbayes_update, gbayes_mix_pred_density, gbayes_mix_pred_cdf
  public :: gbayes_mix_post_density, gbayes_mix_post_cdf, gbayes_mix_post_mean
  public :: gbayes1_power_np, gbayes_mix_power_np
  public :: gbayes2_tabulated, props_po_counts, props_trans_counts, ynbind_numeric
  public :: reformm_missing_order, mapply_group_mean, soprob_markov_ord
  public :: ord_group_boot_mean, clowess_smooth, curve_smooth_observed
  public :: wtd_loess_noiter, sim_markov_ord
  public :: impute_median, impute_constant, mov_stats_n_raw, summarize_mean_codes
  public :: areg_tran_numeric, areg_linear_fit
  public :: dataframe_reduce_numeric, data_rep_numeric, varclus_numeric, redun_numeric
  public :: gbayes_seq, soprob_markov_ordm_core
  public :: areg_boot_linear, fit_mult_impute_combine, impute_transcan_numeric
  public :: est_seq_sim_mean_difference, rm_boot_linear
  public :: ord_test_po_numeric, int_markov_ord_intercepts, sim_reg_ord_supplied

  interface nomiss
    module procedure nomiss_vector
    module procedure nomiss_matrix
  end interface nomiss

  interface solvet
    module procedure solvet_vector
    module procedure solvet_matrix
  end interface solvet

  interface fillin
    module procedure fillin_scalar
    module procedure fillin_vector
  end interface fillin

contains

  pure real(dp) function nan_dp() result(x)
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp


  pure real(dp) function median_sorted_copy(x) result(median)
    real(dp), intent(in) :: x(:) !! Numeric sample whose median is required; caller handles NaNs.
    real(dp), allocatable :: z(:)
    integer :: n
    n = size(x)
    if (n == 0) then
      median = nan_dp()
      return
    end if
    z = x
    call sort_real(z)
    if (mod(n,2) == 1) then
      median = z((n+1)/2)
    else
      median = 0.5_dp * (z(n/2) + z(n/2+1))
    end if
  end function median_sorted_copy

  pure real(dp) function gini_md(x) result(g)
    real(dp), intent(in) :: x(:) !! Sample values; NaNs cause a NaN result.
    real(dp), allocatable :: y(:)
    real(dp) :: xbar
    integer :: i, n
    n = size(x)
    if (n < 2 .or. any(ieee_is_nan(x))) then
      g = nan_dp()
      return
    end if
    y = x
    call sort_real(y)
    xbar = sum(y) / real(n, dp)
    g = 0.0_dp
    do i = 1, n
      g = g + 4.0_dp * (real(i, dp) - real(n - 1, dp) / 2.0_dp) / &
          real(n * (n - 1), dp) * (y(i) - xbar)
    end do
  end function gini_md

  pure real(dp) function trap_rule(x, y) result(area)
    real(dp), intent(in) :: x(:) !! Abscissae, in integration order.
    real(dp), intent(in) :: y(:) !! Ordinates corresponding one-to-one with x.
    integer :: i, n
    n = min(size(x), size(y))
    area = 0.0_dp
    do i = 1, n - 1
      area = area + (x(i + 1) - x(i)) * (y(i + 1) + y(i)) / 2.0_dp
    end do
  end function trap_rule

  pure real(dp) function inv_norm_cdf(p) result(x)
    real(dp), intent(in) :: p !! Probability strictly between zero and one.
    real(dp), parameter :: a1 = -3.969683028665376e1_dp, a2 = 2.209460984245205e2_dp
    real(dp), parameter :: a3 = -2.759285104469687e2_dp, a4 = 1.383577518672690e2_dp
    real(dp), parameter :: a5 = -3.066479806614716e1_dp, a6 = 2.506628277459239_dp
    real(dp), parameter :: b1 = -5.447609879822406e1_dp, b2 = 1.615858368580409e2_dp
    real(dp), parameter :: b3 = -1.556989798598866e2_dp, b4 = 6.680131188771972e1_dp
    real(dp), parameter :: b5 = -1.328068155288572e1_dp
    real(dp), parameter :: c1 = -7.784894002430293e-3_dp, c2 = -3.223964580411365e-1_dp
    real(dp), parameter :: c3 = -2.400758277161838_dp, c4 = -2.549732539343734_dp
    real(dp), parameter :: c5 = 4.374664141464968_dp, c6 = 2.938163982698783_dp
    real(dp), parameter :: d1 = 7.784695709041462e-3_dp, d2 = 3.224671290700398e-1_dp
    real(dp), parameter :: d3 = 2.445134137142996_dp, d4 = 3.754408661907416_dp
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r
    if (p <= 0.0_dp .or. p >= 1.0_dp) then
      x = nan_dp()
    else if (p < plow) then
      q = sqrt(-2.0_dp * log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q * q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp * log(1.0_dp - p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
  end function inv_norm_cdf

  pure real(dp) function norm_cdf(x) result(p)
    real(dp), intent(in) :: x !! Standard-normal variate.
    p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
  end function norm_cdf

  pure subroutine dual_sd(x, nmin, center, bottom, top)
    real(dp), intent(in) :: x(:) !! Complete numeric sample; NaNs yield NaN outputs.
    integer, intent(in) :: nmin !! Minimum sample size for split SDs; smaller samples use the ordinary SD twice.
    real(dp), intent(in) :: center !! Split point; observations equal to center enter both halves.
    real(dp), intent(out) :: bottom !! SD about the full-sample mean for observations at or below center.
    real(dp), intent(out) :: top !! SD about the full-sample mean for observations at or above center.
    real(dp) :: xbar, s
    integer :: n, nb, nt
    n = size(x)
    if (n < 2 .or. any(ieee_is_nan(x))) then
      bottom = nan_dp()
      top = nan_dp()
      return
    end if
    xbar = sum(x) / real(n, dp)
    if (n < nmin) then
      s = sqrt(sum((x - xbar)**2) / real(n - 1, dp))
      bottom = s
      top = s
      return
    end if
    nb = count(x <= center)
    nt = count(x >= center)
    if (nb < 2 .or. nt < 2) then
      bottom = nan_dp()
      top = nan_dp()
      return
    end if
    bottom = sqrt(sum(merge((x - xbar)**2, 0.0_dp, x <= center)) / real(nb - 1, dp))
    top = sqrt(sum(merge((x - xbar)**2, 0.0_dp, x >= center)) / real(nt - 1, dp))
  end subroutine dual_sd

  pure real(dp) function bpower(p1, p2, n1, n2, alpha) result(power)
    real(dp), intent(in) :: p1 !! Event probability in group 1, in [0,1].
    real(dp), intent(in) :: p2 !! Event probability in group 2, in [0,1].
    real(dp), intent(in) :: n1 !! Group-1 sample size, positive.
    real(dp), intent(in) :: n2 !! Group-2 sample size, positive.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability in (0,1).
    real(dp) :: z, pm, ds, ex, sd
    if (p1 < 0.0_dp .or. p1 > 1.0_dp .or. p2 < 0.0_dp .or. p2 > 1.0_dp .or. &
        n1 <= 0.0_dp .or. n2 <= 0.0_dp .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      power = nan_dp()
      return
    end if
    z = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    pm = (n1*p1 + n2*p2) / (n1 + n2)
    ds = z * sqrt((1.0_dp/n1 + 1.0_dp/n2) * pm * (1.0_dp - pm))
    ex = abs(p1 - p2)
    sd = sqrt(p1*(1.0_dp-p1)/n1 + p2*(1.0_dp-p2)/n2)
    if (sd == 0.0_dp) then
      power = merge(1.0_dp, 0.0_dp, ex > ds)
    else
      power = 1.0_dp - norm_cdf((ds - ex)/sd) + norm_cdf((-ds - ex)/sd)
    end if
  end function bpower

  pure subroutine bsamsize(p1, p2, fraction, alpha, power, n1, n2)
    real(dp), intent(in) :: p1 !! Event probability in group 1, in [0,1].
    real(dp), intent(in) :: p2 !! Event probability in group 2, in [0,1].
    real(dp), intent(in) :: fraction !! Fraction allocated to group 1, strictly between 0 and 1.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability in (0,1).
    real(dp), intent(in) :: power !! Target power in (0,1).
    real(dp), intent(out) :: n1 !! Continuous sample-size requirement for group 1.
    real(dp), intent(out) :: n2 !! Continuous sample-size requirement for group 2.
    real(dp) :: za, zb, ratio, pbar
    if (fraction <= 0.0_dp .or. fraction >= 1.0_dp .or. p1 == p2 .or. &
        p1 < 0.0_dp .or. p1 > 1.0_dp .or. p2 < 0.0_dp .or. p2 > 1.0_dp .or. &
        alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. power <= 0.0_dp .or. power >= 1.0_dp) then
      n1 = nan_dp()
      n2 = nan_dp()
      return
    end if
    za = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    zb = inv_norm_cdf(power)
    ratio = (1.0_dp - fraction) / fraction
    pbar = fraction*p1 + (1.0_dp - fraction)*p2
    n1 = (za*sqrt((ratio + 1.0_dp)*pbar*(1.0_dp-pbar)) + &
          zb*sqrt(ratio*p1*(1.0_dp-p1) + p2*(1.0_dp-p2)))**2 / &
          (ratio*(p1-p2)**2)
    n2 = ratio*n1
  end subroutine bsamsize

  pure subroutine ballocation(p1, p2, f_minvar_diff, f_minvar_ratio, f_minvar_logodds, n, f_maxpower)
    real(dp), intent(in) :: p1 !! Event probability in group 1, strictly between 0 and 1.
    real(dp), intent(in) :: p2 !! Event probability in group 2, strictly between 0 and 1.
    real(dp), intent(out) :: f_minvar_diff !! Group-1 fraction minimizing variance of the risk difference.
    real(dp), intent(out) :: f_minvar_ratio !! Group-1 fraction minimizing variance of the risk ratio.
    real(dp), intent(out) :: f_minvar_logodds !! Group-1 fraction minimizing variance of the log odds ratio.
    real(dp), intent(in), optional :: n !! Optional total sample size used to search for the allocation maximizing power.
    real(dp), intent(out), optional :: f_maxpower !! Grid-search allocation maximizing bpower; defined when n is present.
    real(dp) :: q1, q2, f, pw, bestpw, bestdist
    integer :: i
    if (p1 <= 0.0_dp .or. p1 >= 1.0_dp .or. p2 <= 0.0_dp .or. p2 >= 1.0_dp) then
      f_minvar_diff = nan_dp()
      f_minvar_ratio = nan_dp()
      f_minvar_logodds = nan_dp()
      if (present(f_maxpower)) f_maxpower = nan_dp()
      return
    end if
    q1 = 1.0_dp - p1
    q2 = 1.0_dp - p2
    f_minvar_diff = 1.0_dp / (1.0_dp + sqrt(p2*q2/(p1*q1)))
    f_minvar_ratio = 1.0_dp / (1.0_dp + sqrt(p1*q2/(p2*q1)))
    f_minvar_logodds = 1.0_dp - f_minvar_diff
    if (present(f_maxpower)) f_maxpower = nan_dp()
    if (present(n) .and. present(f_maxpower)) then
      if (n <= 0.0_dp) return
      bestpw = -1.0_dp
      bestdist = huge(1.0_dp)
      do i = 1, 1000
        f = 0.001_dp + real(i-1,dp)*(0.998_dp/999.0_dp)
        pw = bpower(p1, p2, n*f, n*(1.0_dp-f), 0.05_dp)
        if (pw > bestpw + 1.0e-15_dp) then
          bestpw = pw
          bestdist = abs(f-0.5_dp)
          f_maxpower = f
        else if (abs(pw-bestpw) <= 1.0e-15_dp .and. abs(f-0.5_dp) < bestdist) then
          bestdist = abs(f-0.5_dp)
          f_maxpower = f
        end if
      end do
    end if
  end subroutine ballocation

  pure subroutine weibull2_fit(times, surv, alpha, gamma)
    real(dp), intent(in) :: times(2) !! Two positive times defining the survival curve.
    real(dp), intent(in) :: surv(2) !! Survival probabilities in (0,1) at the two times.
    real(dp), intent(out) :: alpha !! Weibull scale coefficient in exp(-alpha*t**gamma).
    real(dp), intent(out) :: gamma !! Weibull power parameter.
    real(dp) :: z1, z2
    if (any(times <= 0.0_dp) .or. any(surv <= 0.0_dp) .or. any(surv >= 1.0_dp) .or. times(1) == times(2)) then
      alpha = nan_dp()
      gamma = nan_dp()
      return
    end if
    z1 = -log(surv(1))
    z2 = -log(surv(2))
    gamma = log(z2/z1) / log(times(2)/times(1))
    alpha = z1 / times(1)**gamma
  end subroutine weibull2_fit

  pure elemental real(dp) function weibull_survival(time, alpha, gamma) result(s)
    real(dp), intent(in) :: time !! Nonnegative evaluation time.
    real(dp), intent(in) :: alpha !! Weibull scale coefficient.
    real(dp), intent(in) :: gamma !! Weibull power parameter.
    s = exp(-alpha*time**gamma)
  end function weibull_survival

  pure subroutine gompertz2_fit(times, surv, a, b)
    real(dp), intent(in) :: times(2) !! Two distinct times defining the survival curve.
    real(dp), intent(in) :: surv(2) !! Survival probabilities in (0,1) at the two times.
    real(dp), intent(out) :: a !! Gompertz log-scale parameter.
    real(dp), intent(out) :: b !! Gompertz slope parameter.
    real(dp) :: z1, z2
    if (any(surv <= 0.0_dp) .or. any(surv >= 1.0_dp) .or. times(1) == times(2)) then
      a = nan_dp()
      b = nan_dp()
      return
    end if
    z1 = log(-log(surv(1)))
    z2 = log(-log(surv(2)))
    b = (z2-z1) / (times(2)-times(1))
    if (b <= 0.0_dp) then
      a = nan_dp()
      b = nan_dp()
      return
    end if
    a = z1 + log(b) - b*times(1)
  end subroutine gompertz2_fit

  pure elemental real(dp) function gompertz_survival(time, a, b) result(s)
    real(dp), intent(in) :: time !! Evaluation time.
    real(dp), intent(in) :: a !! Gompertz log-scale parameter.
    real(dp), intent(in) :: b !! Positive Gompertz slope parameter.
    s = exp(-exp(a + b*time)/b)
  end function gompertz_survival

  pure subroutine lognorm2_fit(times, surv, mu, sigma)
    real(dp), intent(in) :: times(2) !! Two positive times defining the survival curve.
    real(dp), intent(in) :: surv(2) !! Survival probabilities in (0,1) at the two times.
    real(dp), intent(out) :: mu !! Mean of log time for the fitted lognormal distribution.
    real(dp), intent(out) :: sigma !! Positive SD of log time for the fitted lognormal distribution.
    real(dp) :: z1, z2
    if (any(times <= 0.0_dp) .or. any(surv <= 0.0_dp) .or. any(surv >= 1.0_dp)) then
      mu = nan_dp()
      sigma = nan_dp()
      return
    end if
    z1 = inv_norm_cdf(1.0_dp - surv(1))
    z2 = inv_norm_cdf(1.0_dp - surv(2))
    if (z1 == z2) then
      mu = nan_dp()
      sigma = nan_dp()
      return
    end if
    sigma = log(times(2)/times(1)) / (z2-z1)
    mu = log(times(1)) - sigma*z1
  end subroutine lognorm2_fit

  pure elemental real(dp) function lognorm_survival(time, mu, sigma) result(s)
    real(dp), intent(in) :: time !! Positive evaluation time.
    real(dp), intent(in) :: mu !! Mean of log time.
    real(dp), intent(in) :: sigma !! Positive SD of log time.
    if (time <= 0.0_dp .or. sigma <= 0.0_dp) then
      s = nan_dp()
    else
      s = norm_cdf(-(log(time)-mu)/sigma)
    end if
  end function lognorm_survival

  pure subroutine improve_prob(x1, x2, y, nri_event, nri_nonevent, nri)
    real(dp), intent(in) :: x1(:) !! Baseline predicted event probabilities in [0,1].
    real(dp), intent(in) :: x2(:) !! Updated predicted event probabilities in [0,1].
    integer, intent(in) :: y(:) !! Binary outcomes coded 0 or 1.
    real(dp), intent(out) :: nri_event !! Event component: fraction moving up minus fraction moving down.
    real(dp), intent(out) :: nri_nonevent !! Non-event component: fraction moving down minus fraction moving up.
    real(dp), intent(out) :: nri !! Overall continuous net reclassification improvement.
    integer :: i, ne, nn, eup, edown, nup, ndown
    if (size(x1) /= size(x2) .or. size(x1) /= size(y) .or. size(y) == 0 .or. &
        any(y < 0) .or. any(y > 1) .or. any(x1 < 0.0_dp) .or. any(x1 > 1.0_dp) .or. &
        any(x2 < 0.0_dp) .or. any(x2 > 1.0_dp)) then
      nri_event = nan_dp()
      nri_nonevent = nan_dp()
      nri = nan_dp()
      return
    end if
    ne = count(y == 1)
    nn = count(y == 0)
    if (ne == 0 .or. nn == 0) then
      nri_event = nan_dp()
      nri_nonevent = nan_dp()
      nri = nan_dp()
      return
    end if
    eup = 0
    edown = 0
    nup = 0
    ndown = 0
    do i = 1, size(y)
      if (y(i) == 1) then
        if (x2(i) > x1(i)) eup = eup + 1
        if (x2(i) < x1(i)) edown = edown + 1
      else
        if (x2(i) > x1(i)) nup = nup + 1
        if (x2(i) < x1(i)) ndown = ndown + 1
      end if
    end do
    nri_event = real(eup-edown, dp) / real(ne, dp)
    nri_nonevent = real(ndown-nup, dp) / real(nn, dp)
    nri = nri_event + nri_nonevent
  end subroutine improve_prob

  pure real(dp) function samplesize_bin(alpha, power, pit, pic, rho) result(nreq)
    real(dp), intent(in) :: alpha !! One-sided type-I error probability in (0,1).
    real(dp), intent(in) :: power !! Desired power probability in (0,1).
    real(dp), intent(in) :: pit !! Treatment success probability in [0,1].
    real(dp), intent(in) :: pic !! Control success probability in [0,1].
    real(dp), intent(in) :: rho !! Fraction allocated to treatment, strictly between 0 and 1.
    real(dp) :: z
    if (rho <= 0.0_dp .or. rho >= 1.0_dp .or. pit == pic) then
      nreq = nan_dp()
      return
    end if
    z = (inv_norm_cdf(1.0_dp - alpha) + inv_norm_cdf(power)) / &
        (asin(sqrt(pit)) - asin(sqrt(pic)))
    nreq = anint(z*z / (4.0_dp*rho*(1.0_dp-rho)) + 0.5_dp)
  end function samplesize_bin

  pure real(dp) function weighted_mean(x, w) result(m)
    real(dp), intent(in) :: x(:) !! Numeric observations.
    real(dp), intent(in) :: w(:) !! Non-missing weights corresponding to x.
    real(dp) :: sw
    sw = sum(w)
    if (size(x) /= size(w) .or. sw == 0.0_dp) then
      m = nan_dp()
    else
      m = sum(w*x) / sw
    end if
  end function weighted_mean

  pure real(dp) function weighted_variance(x, w, normwt, ml) result(v)
    real(dp), intent(in) :: x(:) !! Numeric observations.
    real(dp), intent(in) :: w(:) !! Weights corresponding to x.
    logical, intent(in) :: normwt !! If true, normalize weights to sum to sample size.
    logical, intent(in) :: ml !! If true, use the maximum-likelihood denominator.
    real(dp), allocatable :: ww(:)
    real(dp) :: sw, m, denom
    if (size(x) /= size(w) .or. size(x) < 2) then
      v = nan_dp()
      return
    end if
    ww = w
    sw = sum(ww)
    if (sw <= 0.0_dp) then
      v = nan_dp()
      return
    end if
    if (normwt) ww = ww * real(size(x), dp) / sw
    sw = sum(ww)
    m = sum(ww*x) / sw
    if (.not. normwt .and. .not. ml) then
      denom = sw - 1.0_dp
    else if (ml) then
      denom = sw
    else
      denom = sw - sum(ww*ww)/sw
    end if
    if (denom <= 0.0_dp) then
      v = nan_dp()
    else
      v = sum(ww*(x-m)**2) / denom
    end if
  end function weighted_variance

  pure subroutine weighted_table(x, w, values, sums, n_unique, normwt)
    real(dp), intent(in) :: x(:) !! Numeric observations to aggregate by exact value.
    real(dp), intent(in) :: w(:) !! Weights corresponding to observations.
    real(dp), intent(out) :: values(:) !! Sorted unique values, first n_unique entries are defined.
    real(dp), intent(out) :: sums(:) !! Summed weights for values, first n_unique entries are defined.
    integer, intent(out) :: n_unique !! Number of unique values written to values and sums.
    logical, intent(in) :: normwt !! Normalize weights to sum to number of observations when true.
    real(dp), allocatable :: xx(:), ww(:)
    integer :: i, n
    n = size(x)
    n_unique = 0
    values = 0.0_dp
    sums = 0.0_dp
    if (n == 0 .or. size(w) /= n .or. size(values) < n .or. size(sums) < n) return
    xx = x
    ww = w
    if (normwt .and. sum(ww) /= 0.0_dp) ww = ww * real(n, dp) / sum(ww)
    call sort_pairs(xx, ww)
    do i = 1, n
      if (n_unique == 0) then
        n_unique = 1
        values(n_unique) = xx(i)
        sums(n_unique) = ww(i)
      else if (xx(i) /= values(n_unique)) then
        n_unique = n_unique + 1
        values(n_unique) = xx(i)
        sums(n_unique) = ww(i)
      else
        sums(n_unique) = sums(n_unique) + ww(i)
      end if
    end do
  end subroutine weighted_table

  pure real(dp) function weighted_quantile(x, w, prob, normwt) result(q)
    real(dp), intent(in) :: x(:) !! Numeric observations.
    real(dp), intent(in) :: w(:) !! Frequency-style nonnegative weights.
    real(dp), intent(in) :: prob !! Quantile probability in [0,1].
    logical, intent(in) :: normwt !! Normalize weights to sum to sample size when true.
    real(dp), allocatable :: values(:), sums(:)
    real(dp) :: order, frac, lo, hi, cum
    integer :: n, nu, i
    n = size(x)
    if (n == 0 .or. size(w) /= n .or. prob < 0.0_dp .or. prob > 1.0_dp) then
      q = nan_dp()
      return
    end if
    allocate(values(n), sums(n))
    call weighted_table(x, w, values, sums, nu, normwt)
    if (nu == 0 .or. sum(sums(:nu)) <= 0.0_dp) then
      q = nan_dp()
      return
    end if
    order = 1.0_dp + (sum(sums(:nu)) - 1.0_dp)*prob
    frac = order - floor(order)
    lo = floor(order)
    hi = min(lo + 1.0_dp, sum(sums(:nu)))
    cum = 0.0_dp
    q = values(nu)
    do i = 1, nu
      cum = cum + sums(i)
      if (cum >= lo) then
        q = (1.0_dp-frac)*values(i)
        exit
      end if
    end do
    cum = 0.0_dp
    do i = 1, nu
      cum = cum + sums(i)
      if (cum >= hi) then
        q = q + frac*values(i)
        exit
      end if
    end do
  end function weighted_quantile

  pure subroutine weighted_rank(x, w, rank_out, normwt)
    real(dp), intent(in) :: x(:) !! Numeric observations.
    real(dp), intent(in) :: w(:) !! Frequency-style weights corresponding to x.
    real(dp), intent(out) :: rank_out(:) !! Weighted midranks corresponding to x.
    logical, intent(in) :: normwt !! Normalize weights to sum to sample size when true.
    real(dp), allocatable :: values(:), sums(:), rr(:)
    integer :: n, nu, i, j
    n = size(x)
    rank_out = nan_dp()
    if (size(w) /= n .or. size(rank_out) /= n) return
    allocate(values(n), sums(n), rr(n))
    call weighted_table(x, w, values, sums, nu, normwt)
    if (nu == 0) return
    rr(1) = sums(1) - 0.5_dp*(sums(1)-1.0_dp)
    do i = 2, nu
      rr(i) = sum(sums(:i)) - 0.5_dp*(sums(i)-1.0_dp)
    end do
    do i = 1, n
      do j = 1, nu
        if (x(i) == values(j)) then
          rank_out(i) = rr(j)
          exit
        end if
      end do
    end do
  end subroutine weighted_rank

  pure subroutine somers2(x, y, w, c_index, dxy)
    real(dp), intent(in) :: x(:) !! Predictor values whose ties are retained.
    integer, intent(in) :: y(:) !! Binary outcomes coded 0 or 1.
    real(dp), intent(in) :: w(:) !! Observation weights; use all ones for unweighted analysis.
    real(dp), intent(out) :: c_index !! Concordance probability.
    real(dp), intent(out) :: dxy !! Somers' Dxy, equal to 2*(C-0.5).
    real(dp), allocatable :: r(:)
    real(dp) :: n, n1, mean_rank
    integer :: i
    if (size(x) /= size(y) .or. size(x) /= size(w) .or. any(y < 0) .or. any(y > 1)) then
      c_index = nan_dp()
      dxy = nan_dp()
      return
    end if
    allocate(r(size(x)))
    call weighted_rank(x, w, r, .false.)
    n = sum(w)
    n1 = 0.0_dp
    mean_rank = 0.0_dp
    do i = 1, size(x)
      if (y(i) == 1) then
        n1 = n1 + w(i)
        mean_rank = mean_rank + w(i)*r(i)
      end if
    end do
    if (n1 <= 0.0_dp .or. n1 >= n) then
      c_index = nan_dp()
      dxy = nan_dp()
      return
    end if
    mean_rank = mean_rank / n1
    c_index = (mean_rank - (n1 + 1.0_dp)/2.0_dp) / (n - n1)
    dxy = 2.0_dp*(c_index - 0.5_dp)
  end subroutine somers2

  pure subroutine rcorr(x, spearman, r, npair)
    real(dp), intent(in) :: x(:,:) !! Data matrix; NaNs are treated as missing pairwise.
    logical, intent(in) :: spearman !! Use Spearman midrank correlation when true, Pearson otherwise.
    real(dp), intent(out) :: r(:,:) !! Symmetric pairwise correlation matrix.
    integer, intent(out) :: npair(:,:) !! Pairwise complete observation counts.
    real(dp), allocatable :: a(:), b(:), ra(:), rb(:)
    integer :: i, j, k, m, n, p
    real(dp) :: ma, mb, va, vb
    n = size(x,1)
    p = size(x,2)
    r = nan_dp()
    npair = 0
    do i = 1, p
      npair(i,i) = count(.not. ieee_is_nan(x(:,i)))
      r(i,i) = 1.0_dp
      do j = i + 1, p
        allocate(a(n), b(n), ra(n), rb(n))
        m = 0
        do k = 1, n
          if (.not. ieee_is_nan(x(k,i)) .and. .not. ieee_is_nan(x(k,j))) then
            m = m + 1
            a(m) = x(k,i)
            b(m) = x(k,j)
          end if
        end do
        npair(i,j) = m
        npair(j,i) = m
        if (m >= 2) then
          if (spearman) then
            call midrank(a(:m), ra(:m))
            call midrank(b(:m), rb(:m))
            a(:m) = ra(:m)
            b(:m) = rb(:m)
          end if
          ma = sum(a(:m))/real(m,dp)
          mb = sum(b(:m))/real(m,dp)
          va = sum((a(:m)-ma)**2)
          vb = sum((b(:m)-mb)**2)
          if (va > 0.0_dp .and. vb > 0.0_dp) r(i,j) = sum((a(:m)-ma)*(b(:m)-mb))/sqrt(va*vb)
          r(j,i) = r(i,j)
        end if
        deallocate(a,b,ra,rb)
      end do
    end do
  end subroutine rcorr

  pure subroutine hoeffd(x, dmat, aadmat, maxadmat, npair)
    real(dp), intent(in) :: x(:,:) !! Data matrix; NaNs are omitted pairwise.
    real(dp), intent(out) :: dmat(:,:) !! Hoeffding D values on raw scale (R wrapper reports 30*D).
    real(dp), intent(out) :: aadmat(:,:) !! Mean absolute empirical dependence discrepancy.
    real(dp), intent(out) :: maxadmat(:,:) !! Maximum absolute empirical dependence discrepancy.
    integer, intent(out) :: npair(:,:) !! Pairwise complete sample sizes.
    real(dp), allocatable :: a(:), b(:)
    real(dp) :: d, aad, mad
    integer :: i, j, k, m, n, p
    n = size(x,1)
    p = size(x,2)
    dmat = nan_dp()
    aadmat = 0.0_dp
    maxadmat = 0.0_dp
    npair = 0
    do i = 1, p
      npair(i,i) = count(.not. ieee_is_nan(x(:,i)))
      dmat(i,i) = 1.0_dp/30.0_dp
      do j = i + 1, p
        allocate(a(n), b(n))
        m = 0
        do k = 1, n
          if (.not. ieee_is_nan(x(k,i)) .and. .not. ieee_is_nan(x(k,j))) then
            m = m + 1
            a(m) = x(k,i)
            b(m) = x(k,j)
          end if
        end do
        npair(i,j) = m
        npair(j,i) = m
        if (m > 4) then
          call hoeff_pair(a(:m), b(:m), d, aad, mad)
          dmat(i,j) = d
          dmat(j,i) = d
          aadmat(i,j) = aad
          aadmat(j,i) = aad
          maxadmat(i,j) = mad
          maxadmat(j,i) = mad
        end if
        deallocate(a,b)
      end do
    end do
  end subroutine hoeffd

  pure subroutine hoeff_pair(x, y, d, aad, mad)
    real(dp), intent(in) :: x(:) !! First variable for one complete pair.
    real(dp), intent(in) :: y(:) !! Second variable for one complete pair.
    real(dp), intent(out) :: d !! Hoeffding D statistic on raw scale.
    real(dp), intent(out) :: aad !! Mean absolute empirical dependence discrepancy.
    real(dp), intent(out) :: mad !! Maximum absolute empirical dependence discrepancy.
    real(dp), allocatable :: rx(:), ry(:), rj(:)
    real(dp) :: q, rr, s, z, ad
    integer :: i, n
    n = size(x)
    allocate(rx(n),ry(n),rj(n))
    call joint_rank(x,y,rx,ry,rj)
    q = 0.0_dp
    rr = 0.0_dp
    s = 0.0_dp
    aad = 0.0_dp
    mad = 0.0_dp
    z = real(n,dp)
    do i = 1,n
      ad = abs(rj(i)/z - (rx(i)/z)*(ry(i)/z))
      aad = aad + ad
      mad = max(mad,ad)
      q = q + (rx(i)-1.0_dp)*(rx(i)-2.0_dp)*(ry(i)-1.0_dp)*(ry(i)-2.0_dp)
      rr = rr + (rx(i)-2.0_dp)*(ry(i)-2.0_dp)*(rj(i)-1.0_dp)
      s = s + (rj(i)-1.0_dp)*(rj(i)-2.0_dp)
    end do
    aad = aad/z
    d = (q - 2.0_dp*(z-2.0_dp)*rr + (z-2.0_dp)*(z-3.0_dp)*s) / &
        (z*(z-1.0_dp)*(z-2.0_dp)*(z-3.0_dp)*(z-4.0_dp))
  end subroutine hoeff_pair

  pure subroutine joint_rank(x,y,rx,ry,rj)
    real(dp), intent(in) :: x(:) !! First variable.
    real(dp), intent(in) :: y(:) !! Second variable.
    real(dp), intent(out) :: rx(:) !! Marginal midranks for x.
    real(dp), intent(out) :: ry(:) !! Marginal midranks for y.
    real(dp), intent(out) :: rj(:) !! Joint ranks used by Hoeffding's statistic.
    integer :: i,j,n
    real(dp) :: cx,cy
    n=size(x)
    do i=1,n
      rx(i)=1.0_dp
      ry(i)=1.0_dp
      rj(i)=1.0_dp
      do j=1,n
        if(j==i) cycle
        if(x(j)<x(i)) then
          cx=1.0_dp
        else if(x(j)==x(i)) then
          cx=0.5_dp
        else
          cx=0.0_dp
        end if
        if(y(j)<y(i)) then
          cy=1.0_dp
        else if(y(j)==y(i)) then
          cy=0.5_dp
        else
          cy=0.0_dp
        end if
        rx(i)=rx(i)+cx
        ry(i)=ry(i)+cy
        rj(i)=rj(i)+cx*cy
      end do
    end do
  end subroutine joint_rank

  pure subroutine rcorr_cens(x,y,event,outx,nrel,nconc,nuncert,c_index,gamma,sd)
    real(dp), intent(in) :: x(:) !! Predictor values.
    real(dp), intent(in) :: y(:) !! Right-censored response times.
    logical, intent(in) :: event(:) !! True for observed events and false for right censoring.
    logical, intent(in) :: outx !! Exclude predictor ties from relevant pairs when true.
    real(dp), intent(out) :: nrel !! Number of ordered relevant comparisons.
    real(dp), intent(out) :: nconc !! Number of concordant ordered comparisons, with half credit for ties.
    real(dp), intent(out) :: nuncert !! Number of ordered uncertain comparisons.
    real(dp), intent(out) :: c_index !! Concordance index.
    real(dp), intent(out) :: gamma !! Dxy = 2*(C-0.5).
    real(dp), intent(out) :: sd !! Quade standard-error estimate used upstream.
    integer :: i,j,n
    real(dp) :: dx,dy,wi,ri,sumr,sumr2,sumw,sumw2,sumrw,tmp
    n=size(x)
    nrel=0.0_dp
 nconc=0.0_dp
 nuncert=0.0_dp
    sumr=0.0_dp
 sumr2=0.0_dp
 sumw=0.0_dp
 sumw2=0.0_dp
 sumrw=0.0_dp
    do i=1,n
      wi=0.0_dp
 ri=0.0_dp
      do j=1,n
        if(j==i) cycle
        dx=x(i)-x(j)
 dy=y(i)-y(j)
        if(dx==0.0_dp .and. outx) cycle
        if((event(i) .and. dy<0.0_dp) .or. (event(i) .and. .not.event(j) .and. dy==0.0_dp)) then
          if(dx<0.0_dp) then
 nconc=nconc+1.0_dp
 wi=wi+1.0_dp
          else if(dx==0.0_dp) then
 nconc=nconc+0.5_dp
          else
 wi=wi-1.0_dp
          end if
          nrel=nrel+1.0_dp
 ri=ri+1.0_dp
        else if((event(j) .and. dy>0.0_dp) .or. (event(j) .and. .not.event(i) .and. dy==0.0_dp)) then
          if(dx>0.0_dp) then
 nconc=nconc+1.0_dp
 wi=wi+1.0_dp
          else if(dx==0.0_dp) then
 nconc=nconc+0.5_dp
          else
 wi=wi-1.0_dp
          end if
          nrel=nrel+1.0_dp
 ri=ri+1.0_dp
        else if(.not.(event(i).and.event(j))) then
          nuncert=nuncert+1.0_dp
        end if
      end do
      sumr=sumr+ri
 sumr2=sumr2+ri*ri
 sumw=sumw+wi
 sumw2=sumw2+wi*wi
 sumrw=sumrw+ri*wi
    end do
    if(nrel==0.0_dp .or. sumr==0.0_dp) then
      c_index=nan_dp()
 gamma=nan_dp()
 sd=nan_dp()
 return
    end if
    c_index=nconc/nrel
    gamma=2.0_dp*(c_index-0.5_dp)
    tmp=sumr2*sumw**2-2.0_dp*sumr*sumw*sumrw+sumw2*sumr**2
    sd=2.0_dp*sqrt(max(tmp,0.0_dp))/(sumr*sumr)
  end subroutine rcorr_cens

  pure subroutine rcorrp_cens(x1,x2,y,event,method,outx,gamma,sd,c12,c21,c1,c2,gamma1,gamma2,nrel,nuncert)
    real(dp), intent(in) :: x1(:) !! First set of predictions.
    real(dp), intent(in) :: x2(:) !! Second set of predictions.
    real(dp), intent(in) :: y(:) !! Right-censored response times.
    logical, intent(in) :: event(:) !! True for observed events.
    integer, intent(in) :: method !! 1 compares prediction changes, 2 compares concordance reversals.
    logical, intent(in) :: outx !! Exclude comparisons tied on both predictors when true.
    real(dp), intent(out) :: gamma !! Paired concordance contrast.
    real(dp), intent(out) :: sd !! Quade standard-error estimate for gamma.
    real(dp), intent(out) :: c12 !! Fraction favoring x1 under the chosen method.
    real(dp), intent(out) :: c21 !! Fraction favoring x2 under the chosen method.
    real(dp), intent(out) :: c1 !! Concordance index for x1.
    real(dp), intent(out) :: c2 !! Concordance index for x2.
    real(dp), intent(out) :: gamma1 !! Dxy for x1.
    real(dp), intent(out) :: gamma2 !! Dxy for x2.
    real(dp), intent(out) :: nrel !! Number of ordered relevant comparisons.
    real(dp), intent(out) :: nuncert !! Number of ordered uncertain comparisons.
    integer :: i,j,n
    real(dp) :: dx1,dx2,dy,wi,ri,nconc1,nconc2,sumr,sumr2,sumw,sumw2,sumrw,sumc,tmp
    n=size(x1)
 nconc1=0.0_dp
 nconc2=0.0_dp
 nrel=0.0_dp
 nuncert=0.0_dp
    sumr=0.0_dp
 sumr2=0.0_dp
 sumw=0.0_dp
 sumw2=0.0_dp
 sumrw=0.0_dp
 sumc=0.0_dp
    do i=1,n
      wi=0.0_dp
 ri=0.0_dp
      do j=1,n
        if(i==j) cycle
        dx1=x1(i)-x1(j)
 dx2=x2(i)-x2(j)
        if(outx .and. dx1==0.0_dp .and. dx2==0.0_dp) cycle
        dy=y(i)-y(j)
        if((event(i).and.dy<0.0_dp) .or. (event(i).and..not.event(j).and.dy==0.0_dp)) then
          nrel=nrel+1.0_dp
 ri=ri+1.0_dp
          nconc1=nconc1+bool_real(dx1<0.0_dp)+0.5_dp*bool_real(dx1==0.0_dp)
          nconc2=nconc2+bool_real(dx2<0.0_dp)+0.5_dp*bool_real(dx2==0.0_dp)
          if(method==1) then
            wi=wi+bool_real(dx1<dx2)-bool_real(dx1>dx2)
 sumc=sumc+bool_real(dx1<dx2)
          else
            wi=wi+bool_real(dx1<0.0_dp.and.dx2>=0.0_dp)-bool_real(dx1>0.0_dp.and.dx2<=0.0_dp)
            sumc=sumc+bool_real(dx1<0.0_dp.and.dx2>=0.0_dp)
          end if
        else if((event(j).and.dy>0.0_dp) .or. (event(j).and..not.event(i).and.dy==0.0_dp)) then
          nrel=nrel+1.0_dp
 ri=ri+1.0_dp
          nconc1=nconc1+bool_real(dx1>0.0_dp)+0.5_dp*bool_real(dx1==0.0_dp)
          nconc2=nconc2+bool_real(dx2>0.0_dp)+0.5_dp*bool_real(dx2==0.0_dp)
          if(method==1) then
            wi=wi+bool_real(dx1>dx2)-bool_real(dx1<dx2)
 sumc=sumc+bool_real(dx1>dx2)
          else
            wi=wi+bool_real(dx1>0.0_dp.and.dx2<=0.0_dp)-bool_real(dx1<0.0_dp.and.dx2>=0.0_dp)
            sumc=sumc+bool_real(dx1>0.0_dp.and.dx2<=0.0_dp)
          end if
        else if(.not.(event(i).and.event(j))) then
          nuncert=nuncert+1.0_dp
        end if
      end do
      sumr=sumr+ri
 sumr2=sumr2+ri*ri
 sumw=sumw+wi
 sumw2=sumw2+wi*wi
 sumrw=sumrw+ri*wi
    end do
    if(nrel==0.0_dp .or. sumr==0.0_dp) then
      gamma=nan_dp()
 sd=nan_dp()
 c12=nan_dp()
 c21=nan_dp()
 c1=nan_dp()
 c2=nan_dp()
 gamma1=nan_dp()
 gamma2=nan_dp()
 return
    end if
    c1=nconc1/nrel
 c2=nconc2/nrel
 gamma1=2.0_dp*(c1-0.5_dp)
 gamma2=2.0_dp*(c2-0.5_dp)
    gamma=sumw/sumr
 c12=sumc/sumr
 c21=c12-gamma
    tmp=sumr2*sumw**2-2.0_dp*sumr*sumw*sumrw+sumw2*sumr**2
    sd=2.0_dp*sqrt(max(tmp,0.0_dp))/(sumr*sumr)
  end subroutine rcorrp_cens

  pure real(dp) function bool_real(a) result(z)
    logical, intent(in) :: a !! Logical condition converted to zero or one.
    if(a) then
      z=1.0_dp
    else
      z=0.0_dp
    end if
  end function bool_real

  pure subroutine cut_gn(x,m,g)
    real(dp), intent(in) :: x(:) !! Sorted numeric values to divide into groups.
    integer, intent(in) :: m !! Target minimum observations per group; must be positive.
    integer, intent(out) :: g(:) !! Group number for each sorted observation.
    integer :: is,ie,ng,k,i,lte,n
    real(dp) :: lastval
    n=size(x)
 ie=0
 g=0
 ng=0
    if(m<=0 .or. size(g)/=n) return
    do
      is=ie+1
 lte=is+m-1
      if(lte>n) then
 if(is<=n) g(is:n)=ng
 return
 end if
      ng=ng+1
      if(lte==n) then
 g(is:n)=ng
 return
 end if
      lastval=x(lte)
      if(x(lte+1)==lastval) then
        k=1
        do i=lte+2,n
          if(x(i)/=lastval) exit
          k=k+1
        end do
        ie=lte+k
      else
        ie=lte
      end if
      g(is:ie)=ng
      if(ie==n) return
    end do
  end subroutine cut_gn

  pure subroutine which_closest(x,w,j)
    real(dp), intent(in) :: x(:) !! Reference values.
    real(dp), intent(in) :: w(:) !! Query values for nearest-reference lookup.
    integer, intent(out) :: j(:) !! One-based nearest indices in x, ties resolved to the first.
    integer :: i,k,m
    real(dp) :: d,dmin
    do i=1,size(w)
      dmin=huge(1.0_dp)
 m=0
      do k=1,size(x)
        d=abs(x(k)-w(i))
        if(d<dmin) then
 dmin=d
 m=k
 end if
      end do
      j(i)=m
    end do
  end subroutine which_closest

  pure subroutine which_close_pw(x,w,r,f,j)
    real(dp), intent(in) :: x(:) !! Reference values.
    real(dp), intent(in) :: w(:) !! Query values.
    real(dp), intent(in) :: r(:) !! Uniform draws in [0,1) supplied explicitly for reproducibility.
    real(dp), intent(in) :: f !! Distance-bandwidth multiplier used by the tricube probabilities.
    integer, intent(out) :: j(:) !! One-based sampled indices in x.
    real(dp), allocatable :: wt(:)
    real(dp) :: dmean,sump,prob,z
    integer :: i,k,m
    allocate(wt(size(x)))
    do i=1,size(w)
      dmean=sum(abs(x-w(i)))*f/real(size(x),dp)
      if(dmean<=0.0_dp) then
        call which_closest(x,w(i:i),j(i:i))
        cycle
      end if
      do k=1,size(x)
        z=min(abs(x(k)-w(i))/dmean,1.0_dp)
        wt(k)=(1.0_dp-z**3)**3
      end do
      sump=sum(wt)
 prob=0.0_dp
 m=1
      if(sump<=0.0_dp) then
 call which_closest(x,w(i:i),j(i:i))
 cycle
 end if
      do k=1,size(x)
        prob=prob+wt(k)/sump
        if(r(i)>prob) m=m+1
      end do
      j(i)=min(m,size(x))
    end do
  end subroutine which_close_pw

  pure subroutine rcspline_eval(x, knots, inclx, integral, norm, basis)
    real(dp), intent(in) :: x(:) !! Predictor values at which to evaluate the restricted cubic spline basis.
    real(dp), intent(in) :: knots(:) !! Strictly increasing unique knot locations; at least three are required.
    logical, intent(in) :: inclx !! Include the original predictor as the first column for an ordinary spline basis.
    logical, intent(in) :: integral !! Return the integrated spline basis, including x and x squared over two.
    integer, intent(in) :: norm !! Basis normalization mode: 0 none, 1 cubic spacing, or 2 squared outer-knot spacing.
    real(dp), allocatable, intent(out) :: basis(:,:) !! Evaluated basis matrix; zero columns signal invalid knots or norm.
    real(dp) :: kd, knot1, knotnk, knotnk1, z1, z2, z3
    integer :: i, j, nk, nc, power, offset

    nk = size(knots)
    if (nk < 3 .or. any(knots(2:nk) <= knots(1:nk-1)) .or. norm < 0 .or. norm > 2) then
      allocate(basis(size(x), 0))
      return
    end if
    knot1 = knots(1)
    knotnk = knots(nk)
    knotnk1 = knots(nk - 1)
    if (knotnk == knotnk1 .or. knotnk == knot1) then
      allocate(basis(size(x), 0))
      return
    end if
    select case (norm)
    case (0)
      kd = 1.0_dp
    case (1)
      kd = knotnk - knotnk1
    case default
      kd = (knotnk - knot1)**(2.0_dp / 3.0_dp)
    end select
    power = merge(4, 3, integral)
    if (integral) then
      nc = nk
      offset = 2
    else if (inclx) then
      nc = nk - 1
      offset = 1
    else
      nc = nk - 2
      offset = 0
    end if
    allocate(basis(size(x), nc))
    basis = 0.0_dp
    if (integral) then
      basis(:, 1) = x
      basis(:, 2) = x*x / 2.0_dp
    else if (inclx) then
      basis(:, 1) = x
    end if
    do j = 1, nk - 2
      do i = 1, size(x)
        if (ieee_is_nan(x(i))) then
          basis(i, j + offset) = nan_dp()
        else
          z1 = max((x(i) - knots(j)) / kd, 0.0_dp)**power
          z2 = max((x(i) - knotnk) / kd, 0.0_dp)**power
          z3 = max((x(i) - knotnk1) / kd, 0.0_dp)**power
          basis(i, j + offset) = z1 + ((knotnk1 - knots(j))*z2 - &
            (knotnk - knots(j))*z3) / (knotnk - knotnk1)
          if (integral) basis(i, j + offset) = basis(i, j + offset) * kd / 4.0_dp
        end if
      end do
    end do
  end subroutine rcspline_eval

  pure subroutine popower(p, odds_ratio, n1, n2, alpha, power, efficiency, approx_se)
    real(dp), intent(in) :: p(:) !! Ordinal-category probabilities, nonnegative and summing to one.
    real(dp), intent(in) :: odds_ratio !! Proportional odds ratio, strictly positive.
    real(dp), intent(in) :: n1 !! Sample size in treatment group one, positive.
    real(dp), intent(in) :: n2 !! Sample size in treatment group two, positive.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability in (0,1).
    real(dp), intent(out) :: power !! Approximate power of the proportional-odds test.
    real(dp), intent(out) :: efficiency !! Efficiency relative to a continuous response.
    real(dp), intent(out) :: approx_se !! Approximate standard error of the log odds ratio.
    real(dp) :: n, ps, v, z

    if (size(p) == 0 .or. any(p < 0.0_dp) .or. abs(sum(p) - 1.0_dp) > 1.0e-5_dp .or. &
        odds_ratio <= 0.0_dp .or. n1 <= 0.0_dp .or. n2 <= 0.0_dp .or. &
        alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      power = nan_dp()
      efficiency = nan_dp()
      approx_se = nan_dp()
      return
    end if
    n = n1 + n2
    ps = 1.0_dp - sum(p**3)
    v = n1*n2*n*ps / (3.0_dp*(n + 1.0_dp)**2)
    if (v <= 0.0_dp) then
      power = nan_dp()
      efficiency = nan_dp()
      approx_se = nan_dp()
      return
    end if
    z = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    power = norm_cdf(abs(log(odds_ratio))*sqrt(v) - z)
    efficiency = ps / (1.0_dp - 1.0_dp/(n*n))
    approx_se = 1.0_dp / sqrt(v)
  end subroutine popower

  pure subroutine posamsize(p, odds_ratio, fraction, alpha, target_power, n, efficiency)
    real(dp), intent(in) :: p(:) !! Ordinal-category probabilities, nonnegative and summing to one.
    real(dp), intent(in) :: odds_ratio !! Proportional odds ratio, positive and different from one.
    real(dp), intent(in) :: fraction !! Fraction allocated to group one, strictly between zero and one.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability in (0,1).
    real(dp), intent(in) :: target_power !! Desired power in (0,1).
    real(dp), intent(out) :: n !! Approximate required total sample size.
    real(dp), intent(out) :: efficiency !! Efficiency relative to a continuous response at the computed n.
    real(dp) :: a, log_or, ps, za, zb

    if (size(p) == 0 .or. any(p < 0.0_dp) .or. abs(sum(p) - 1.0_dp) > 1.0e-5_dp .or. &
        odds_ratio <= 0.0_dp .or. odds_ratio == 1.0_dp .or. fraction <= 0.0_dp .or. fraction >= 1.0_dp .or. &
        alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. target_power <= 0.0_dp .or. target_power >= 1.0_dp) then
      n = nan_dp()
      efficiency = nan_dp()
      return
    end if
    a = (1.0_dp - fraction) / fraction
    log_or = log(odds_ratio)
    za = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    zb = inv_norm_cdf(target_power)
    ps = 1.0_dp - sum(p**3)
    if (ps <= 0.0_dp) then
      n = nan_dp()
      efficiency = nan_dp()
      return
    end if
    n = 3.0_dp*(a + 1.0_dp)**2*(za + zb)**2 / (a*log_or**2*ps)
    efficiency = ps / (1.0_dp - 1.0_dp/(n*n))
  end subroutine posamsize

  pure subroutine pomodm(p, odds_ratio, p_out)
    real(dp), intent(in) :: p(:) !! Baseline individual category probabilities, nonnegative and summing to one.
    real(dp), intent(in) :: odds_ratio !! Odds ratio applied to all exceedance probabilities; must be positive.
    real(dp), intent(out) :: p_out(:) !! Modified category probabilities; same length as p.
    real(dp), allocatable :: ep(:), ep2(:)
    integer :: i, k

    k = size(p)
    if (size(p_out) /= k) return
    if (k == 0 .or. any(p < 0.0_dp) .or. abs(sum(p) - 1.0_dp) > 1.0e-5_dp .or. odds_ratio <= 0.0_dp) then
      p_out = nan_dp()
      return
    end if
    allocate(ep(k), ep2(k))
    ep(1) = 1.0_dp
    do i = 2, k
      ep(i) = 1.0_dp - sum(p(1:i-1))
    end do
    do i = 1, k
      ep2(i) = odds_ratio*ep(i) / (1.0_dp - ep(i) + odds_ratio*ep(i))
    end do
    do i = 1, k - 1
      p_out(i) = ep2(i) - ep2(i + 1)
    end do
    p_out(k) = ep2(k)
  end subroutine pomodm

  pure subroutine groupn(x, y, m, mean_x, mean_y, status)
    real(dp), intent(in) :: x(:) !! First numeric variable; pairs with missing x or y are discarded.
    real(dp), intent(in) :: y(:) !! Second numeric variable paired elementwise with x.
    integer, intent(in) :: m !! Number of sorted observations per group; must not exceed complete-pair count.
    real(dp), allocatable, intent(out) :: mean_x(:) !! Successive group means, retaining the upstream final overlap behavior.
    real(dp), allocatable, intent(out) :: mean_y(:) !! Means of y in the same groups as mean_x.
    integer, intent(out), optional :: status !! Zero on success; nonzero for bad lengths, invalid m, or too few complete pairs.
    real(dp), allocatable :: xx(:), yy(:)
    integer, allocatable :: ord(:)
    integer :: i, j, n, ng, start_i, end_i, ncomplete

    if (present(status)) status = 0
    if (size(x) /= size(y) .or. m <= 0) then
      allocate(mean_x(0), mean_y(0))
      if (present(status)) status = 1
      return
    end if
    ncomplete = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (ncomplete < m) then
      allocate(mean_x(0), mean_y(0))
      if (present(status)) status = 2
      return
    end if
    allocate(xx(ncomplete), yy(ncomplete), ord(ncomplete))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        j = j + 1
        xx(j) = x(i)
        yy(j) = y(i)
      end if
    end do
    call sort_index_real(xx, ord)
    xx = xx(ord)
    yy = yy(ord)
    n = ncomplete
    ng = n / m + 1
    allocate(mean_x(ng), mean_y(ng))
    j = 0
    start_i = 1
    end_i = m
    do while (end_i <= n)
      j = j + 1
      mean_x(j) = sum(xx(start_i:end_i)) / real(m, dp)
      mean_y(j) = sum(yy(start_i:end_i)) / real(m, dp)
      start_i = start_i + m
      end_i = end_i + m
    end do
    j = j + 1
    mean_x(j) = sum(xx(n-m+1:n)) / real(m, dp)
    mean_y(j) = sum(yy(n-m+1:n)) / real(m, dp)
  end subroutine groupn

  pure real(dp) function spearman_rho(x, y) result(rho)
    real(dp), intent(in) :: x(:) !! First numeric vector; pairs containing NaN are discarded.
    real(dp), intent(in) :: y(:) !! Second numeric vector paired elementwise with x.
    real(dp), allocatable :: xx(:), yy(:), rx(:), ry(:)
    real(dp) :: mx, my, sx, sy
    integer :: i, j, n

    if (size(x) /= size(y)) then
      rho = nan_dp()
      return
    end if
    n = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (n < 3) then
      rho = nan_dp()
      return
    end if
    allocate(xx(n), yy(n), rx(n), ry(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        j = j + 1
        xx(j) = x(i)
        yy(j) = y(i)
      end if
    end do
    call midrank(xx, rx)
    call midrank(yy, ry)
    mx = sum(rx) / real(n, dp)
    my = sum(ry) / real(n, dp)
    sx = sqrt(sum((rx-mx)**2))
    sy = sqrt(sum((ry-my)**2))
    if (sx == 0.0_dp .or. sy == 0.0_dp) then
      rho = nan_dp()
    else
      rho = sum((rx-mx)*(ry-my)) / (sx*sy)
    end if
  end function spearman_rho

  pure subroutine stepfun_eval(x, y, xout, right, yout)
    real(dp), intent(in) :: x(:) !! Step-function abscissae; NaN pairs are discarded and remaining x must be increasing.
    real(dp), intent(in) :: y(:) !! Step-function ordinates corresponding to x.
    real(dp), intent(in) :: xout(:) !! Abscissae at which to evaluate the constant interpolation.
    logical, intent(in) :: right !! Use upstream type='right' convention when true; otherwise type='left'.
    real(dp), intent(out) :: yout(:) !! Interpolated ordinates; points outside the observed x range are NaN.
    real(dp), allocatable :: xx(:), yy(:)
    integer :: i, j, k, n

    if (size(yout) /= size(xout)) return
    yout = nan_dp()
    if (size(x) /= size(y)) return
    n = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (n == 0) return
    allocate(xx(n), yy(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        j = j + 1
        xx(j) = x(i)
        yy(j) = y(i)
      end if
    end do
    if (n > 1) then
      if (any(xx(2:n) <= xx(1:n-1))) return
    end if
    do k = 1, size(xout)
      if (ieee_is_nan(xout(k)) .or. xout(k) < xx(1) .or. xout(k) > xx(n)) cycle
      if (xout(k) == xx(n)) then
        yout(k) = yy(n)
        cycle
      end if
      do i = 1, n - 1
        if (xout(k) == xx(i)) then
          yout(k) = yy(i)
          exit
        else if (xout(k) > xx(i) .and. xout(k) < xx(i+1)) then
          if (right) then
            yout(k) = yy(i+1)
          else
            yout(k) = yy(i)
          end if
          exit
        end if
      end do
    end do
  end subroutine stepfun_eval

  pure subroutine xy_sort_no_dup_no_na(x, y, x_out, y_out)
    real(dp), intent(in) :: x(:) !! First numeric vector; NaN pairs are removed before sorting.
    real(dp), intent(in) :: y(:) !! Second numeric vector paired elementwise with x.
    real(dp), allocatable, intent(out) :: x_out(:) !! Sorted unique x values, retaining the first y after sorting for duplicate x.
    real(dp), allocatable, intent(out) :: y_out(:) !! y values paired with x_out.
    real(dp), allocatable :: xx(:), yy(:)
    integer, allocatable :: ord(:)
    integer :: i, j, n, nu

    if (size(x) /= size(y)) then
      allocate(x_out(0), y_out(0))
      return
    end if
    n = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (n == 0) then
      allocate(x_out(0), y_out(0))
      return
    end if
    allocate(xx(n), yy(n), ord(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        j = j + 1
        xx(j) = x(i)
        yy(j) = y(i)
      end if
    end do
    call sort_index_real(xx, ord)
    xx = xx(ord)
    yy = yy(ord)
    nu = 1
    do i = 2, n
      if (xx(i) /= xx(i-1)) nu = nu + 1
    end do
    allocate(x_out(nu), y_out(nu))
    j = 1
    x_out(1) = xx(1)
    y_out(1) = yy(1)
    do i = 2, n
      if (xx(i) /= xx(i-1)) then
        j = j + 1
        x_out(j) = xx(i)
        y_out(j) = yy(i)
      end if
    end do
  end subroutine xy_sort_no_dup_no_na

  pure integer function n_coincident(x, y, bins) result(nc)
    real(dp), intent(in) :: x(:) !! Scatterplot x coordinates; pairs containing NaN are discarded.
    real(dp), intent(in) :: y(:) !! Scatterplot y coordinates paired elementwise with x.
    integer, intent(in) :: bins !! Compatibility argument; upstream hard-codes 300-by-300 scaling and ignores this value.
    integer, allocatable :: bx(:), by(:)
    real(dp) :: xmin, xmax, ymin, ymax
    integer :: i, j, n, nu

    nc = 0
    if (size(x) /= size(y) .or. bins <= 0) return
    n = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (n <= 1) return
    xmin = minval(x, mask=.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    xmax = maxval(x, mask=.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    ymin = minval(y, mask=.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    ymax = maxval(y, mask=.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (xmax == xmin .or. ymax == ymin) return
    allocate(bx(n), by(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        j = j + 1
        bx(j) = nint((x(i)-xmin)/(xmax-xmin)*300.0_dp)
        by(j) = nint((y(i)-ymin)/(ymax-ymin)*300.0_dp)
      end if
    end do
    nu = 0
    do i = 1, n
      if (i == 1) then
        nu = 1
      else
        do j = 1, i - 1
          if (bx(j) == bx(i) .and. by(j) == by(i)) exit
        end do
        if (j == i) nu = nu + 1
      end if
    end do
    nc = n - nu
  end function n_coincident

  pure subroutine matxv(a, b, result, kint)
    real(dp), intent(in) :: a(:,:) !! Numeric design matrix to multiply by the trailing coefficients of b.
    real(dp), intent(in) :: b(:) !! Coefficient vector; may contain leading intercept coefficients beyond the columns of a.
    real(dp), intent(out) :: result(:) !! Matrix-vector product plus selected intercept; length must equal rows of a.
    integer, intent(in), optional :: kint !! Optional one-based intercept index when b is longer than the column count of a.
    integer :: nc, lb, excess, ki
    real(dp) :: intercept

    if (size(result) /= size(a,1)) return
    nc = size(a,2)
    lb = size(b)
    if (lb < nc) then
      result = nan_dp()
      return
    end if
    intercept = 0.0_dp
    if (lb > nc .and. present(kint)) then
      excess = lb - nc
      ki = kint
      if (ki < 1 .or. ki > excess) then
        result = nan_dp()
        return
      end if
      intercept = b(ki)
    end if
    result = intercept + matmul(a, b(lb-nc+1:lb))
  end subroutine matxv


  pure real(dp) function cpower(tref, n, mc, reduction, accrual, tmin, noncomp_c, noncomp_i, alpha, nc, ni) result(power)
    real(dp), intent(in) :: tref !! Reference time at which control mortality is specified; must be positive.
    real(dp), intent(in) :: n !! Total sample size used when nc and ni are not supplied; must be positive.
    real(dp), intent(in) :: mc !! Control mortality probability at tref, strictly between zero and one.
    real(dp), intent(in) :: reduction !! Percent reduction in mortality for the intervention arm.
    real(dp), intent(in) :: accrual !! Accrual duration in the same time units as tref; nonnegative.
    real(dp), intent(in) :: tmin !! Minimum follow-up time in the same time units as tref; positive.
    real(dp), intent(in), optional :: noncomp_c !! Control-arm noncompliance percent; defaults to zero.
    real(dp), intent(in), optional :: noncomp_i !! Intervention-arm noncompliance percent; defaults to zero.
    real(dp), intent(in), optional :: alpha !! Two-sided type-I error probability; defaults to 0.05.
    real(dp), intent(in), optional :: nc !! Optional control-arm sample size; paired with ni when supplied.
    real(dp), intent(in), optional :: ni !! Optional intervention-arm sample size; paired with nc when supplied.
    real(dp) :: mi, lamc, lami, pc, pi, ec, ei, delta, sd, z, aa, ncc, nii, nccmp, nicmp, tmax

    power = nan_dp()
    if (tref <= 0.0_dp .or. n <= 0.0_dp .or. mc <= 0.0_dp .or. mc >= 1.0_dp .or. tmin <= 0.0_dp .or. accrual < 0.0_dp) return
    aa = 0.05_dp
    if (present(alpha)) aa = alpha
    if (aa <= 0.0_dp .or. aa >= 1.0_dp) return
    nccmp = 0.0_dp
    nicmp = 0.0_dp
    if (present(noncomp_c)) nccmp = noncomp_c
    if (present(noncomp_i)) nicmp = noncomp_i
    if (present(nc) .neqv. present(ni)) return
    if (present(nc)) then
      ncc = nc
      nii = ni
    else
      ncc = n / 2.0_dp
      nii = n / 2.0_dp
    end if
    if (ncc <= 0.0_dp .or. nii <= 0.0_dp) return
    mi = (1.0_dp - reduction / 100.0_dp) * mc
    if (mi <= 0.0_dp .or. mi >= 1.0_dp) return
    lamc = -log(1.0_dp - mc) / tref
    lami = -log(1.0_dp - mi) / tref
    tmax = tmin + accrual
    if (accrual == 0.0_dp) then
      pc = 1.0_dp - exp(-lamc * tmin)
      pi = 1.0_dp - exp(-lami * tmin)
    else
      pc = 1.0_dp - (exp(-tmin*lamc) - exp(-tmax*lamc)) / (accrual*lamc)
      pi = 1.0_dp - (exp(-tmin*lami) - exp(-tmax*lami)) / (accrual*lami)
    end if
    ec = pc * ncc
    ei = pi * nii
    if (ec <= 0.0_dp .or. ei <= 0.0_dp) return
    delta = log(lami / lamc)
    delta = delta * (1.0_dp - (nccmp + nicmp) / 100.0_dp)
    sd = sqrt(1.0_dp / ec + 1.0_dp / ei)
    z = inv_norm_cdf(1.0_dp - aa / 2.0_dp)
    power = 1.0_dp - (norm_cdf(z - abs(delta)/sd) - norm_cdf(-z - abs(delta)/sd))
  end function cpower

  pure real(dp) function ciapower(tref, n1, n2, m1c, m2c, r1, r2, accrual, tmin, alpha) result(power)
    real(dp), intent(in) :: tref !! Reference time at which stratum-specific control mortalities are specified; positive.
    real(dp), intent(in) :: n1 !! Total sample size in stratum 1, split equally between treatment arms.
    real(dp), intent(in) :: n2 !! Total sample size in stratum 2, split equally between treatment arms.
    real(dp), intent(in) :: m1c !! Stratum-1 control mortality probability at tref, in (0,1).
    real(dp), intent(in) :: m2c !! Stratum-2 control mortality probability at tref, in (0,1).
    real(dp), intent(in) :: r1 !! Percent mortality reduction in stratum 1 intervention subjects.
    real(dp), intent(in) :: r2 !! Percent mortality reduction in stratum 2 intervention subjects.
    real(dp), intent(in) :: accrual !! Accrual duration; positive for the upstream closed-form calculation.
    real(dp), intent(in) :: tmin !! Minimum follow-up duration; positive.
    real(dp), intent(in), optional :: alpha !! Two-sided type-I error probability; defaults to 0.05.
    real(dp) :: m1i, m2i, l1c, l2c, l1i, l2i, p1c, p2c, p1i, p2i
    real(dp) :: e1c, e2c, e1i, e2i, delta, sd, z, aa, tmax

    power = nan_dp()
    if (tref <= 0.0_dp .or. n1 <= 0.0_dp .or. n2 <= 0.0_dp .or. accrual <= 0.0_dp .or. tmin <= 0.0_dp) return
    if (m1c <= 0.0_dp .or. m1c >= 1.0_dp .or. m2c <= 0.0_dp .or. m2c >= 1.0_dp) return
    aa = 0.05_dp
    if (present(alpha)) aa = alpha
    if (aa <= 0.0_dp .or. aa >= 1.0_dp) return
    m1i = (1.0_dp - r1/100.0_dp) * m1c
    m2i = (1.0_dp - r2/100.0_dp) * m2c
    if (m1i <= 0.0_dp .or. m1i >= 1.0_dp .or. m2i <= 0.0_dp .or. m2i >= 1.0_dp) return
    l1c = -log(1.0_dp-m1c)/tref
    l2c = -log(1.0_dp-m2c)/tref
    l1i = -log(1.0_dp-m1i)/tref
    l2i = -log(1.0_dp-m2i)/tref
    tmax = tmin + accrual
    p1c = 1.0_dp - (exp(-tmin*l1c)-exp(-tmax*l1c))/(accrual*l1c)
    p2c = 1.0_dp - (exp(-tmin*l2c)-exp(-tmax*l2c))/(accrual*l2c)
    p1i = 1.0_dp - (exp(-tmin*l1i)-exp(-tmax*l1i))/(accrual*l1i)
    p2i = 1.0_dp - (exp(-tmin*l2i)-exp(-tmax*l2i))/(accrual*l2i)
    e1c = p1c*n1/2.0_dp
    e2c = p2c*n2/2.0_dp
    e1i = p1i*n1/2.0_dp
    e2i = p2i*n2/2.0_dp
    if (min(e1c,e2c,e1i,e2i) <= 0.0_dp) return
    delta = log((l1i/l1c)/(l2i/l2c))
    sd = sqrt(1.0_dp/e1c + 1.0_dp/e2c + 1.0_dp/e1i + 1.0_dp/e2i)
    z = inv_norm_cdf(1.0_dp-aa/2.0_dp)
    power = 1.0_dp - (norm_cdf(z-abs(delta)/sd) - norm_cdf(-z-abs(delta)/sd))
  end function ciapower

  pure subroutine wtd_ecdf(x, weights, ecdf_type, xout, fout)
    real(dp), intent(in) :: x(:) !! Numeric observations; NaNs are omitted together with their weights.
    real(dp), intent(in), optional :: weights(:) !! Optional nonnegative observation weights; unit weights are used when absent.
    integer, intent(in), optional :: ecdf_type !! 1 for i/n, 2 for (i-1)/(n-1), or 3 for i/(n+1); defaults to 1.
    real(dp), allocatable, intent(out) :: xout(:) !! Sorted unique support, repeating the first point when a leading zero is needed.
    real(dp), allocatable, intent(out) :: fout(:) !! Weighted empirical CDF values corresponding to xout.
    real(dp), allocatable :: xx(:), ww(:), ux(:), uw(:)
    real(dp) :: a, b, total, c
    integer :: i, j, n, nu, typ, lead

    if (present(weights)) then
      if (size(weights) /= size(x)) then
        allocate(xout(0), fout(0))
        return
      end if
      n = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(weights) .and. weights /= 0.0_dp)
    else
      n = count(.not. ieee_is_nan(x))
    end if
    if (n == 0) then
      allocate(xout(0), fout(0))
      return
    end if
    allocate(xx(n), ww(n))
    j = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i))) cycle
      if (present(weights)) then
        if (ieee_is_nan(weights(i)) .or. weights(i) == 0.0_dp) cycle
        j = j + 1
        xx(j) = x(i)
        ww(j) = weights(i)
      else
        j = j + 1
        xx(j) = x(i)
        ww(j) = 1.0_dp
      end if
    end do
    call sort_pairs(xx, ww)
    nu = 1
    do i = 2, n
      if (xx(i) /= xx(i-1)) nu = nu + 1
    end do
    allocate(ux(nu), uw(nu))
    j = 1
    ux(1) = xx(1)
    uw(1) = ww(1)
    do i = 2, n
      if (xx(i) == xx(i-1)) then
        uw(j) = uw(j) + ww(i)
      else
        j = j + 1
        ux(j) = xx(i)
        uw(j) = ww(i)
      end if
    end do
    typ = 1
    if (present(ecdf_type)) typ = ecdf_type
    select case (typ)
    case (1)
      a = 0.0_dp
      b = 0.0_dp
    case (2)
      a = -1.0_dp
      b = -1.0_dp
    case (3)
      a = 0.0_dp
      b = 1.0_dp
    case default
      allocate(xout(0), fout(0))
      return
    end select
    total = sum(uw)
    if (total + b == 0.0_dp) then
      allocate(xout(0), fout(0))
      return
    end if
    c = (uw(1)+a)/(total+b)
    lead = merge(1, 0, c > 0.0_dp)
    allocate(xout(nu+lead), fout(nu+lead))
    if (lead == 1) then
      xout(1) = ux(1)
      fout(1) = 0.0_dp
    end if
    c = 0.0_dp
    do i = 1, nu
      c = c + uw(i)
      xout(i+lead) = ux(i)
      fout(i+lead) = (c+a)/(total+b)
    end do
  end subroutine wtd_ecdf

  pure subroutine smean_sd(x, mean_value, sd_value)
    real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted as with upstream na.rm=TRUE.
    real(dp), intent(out) :: mean_value !! Arithmetic mean of complete observations, or NaN when none exist.
    real(dp), intent(out) :: sd_value !! Sample standard deviation, or NaN when fewer than two complete observations exist.
    real(dp) :: s, ss
    integer :: i, n
    n = count(.not. ieee_is_nan(x))
    if (n == 0) then
      mean_value = nan_dp()
      sd_value = nan_dp()
      return
    end if
    s = 0.0_dp
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) s = s + x(i)
    end do
    mean_value = s / real(n,dp)
    if (n < 2) then
      sd_value = nan_dp()
      return
    end if
    ss = 0.0_dp
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) ss = ss + (x(i)-mean_value)**2
    end do
    sd_value = sqrt(ss/real(n-1,dp))
  end subroutine smean_sd

  pure subroutine smean_sdl(x, mult, mean_value, lower, upper)
    real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted as with upstream na.rm=TRUE.
    real(dp), intent(in), optional :: mult !! Number of sample SDs defining each limit; defaults to 2.
    real(dp), intent(out) :: mean_value !! Arithmetic mean of complete observations.
    real(dp), intent(out) :: lower !! Mean minus mult times sample SD.
    real(dp), intent(out) :: upper !! Mean plus mult times sample SD.
    real(dp) :: sdv, m
    m = 2.0_dp
    if (present(mult)) m = mult
    call smean_sd(x, mean_value, sdv)
    if (ieee_is_nan(mean_value) .or. ieee_is_nan(sdv)) then
      lower = nan_dp()
      upper = nan_dp()
    else
      lower = mean_value - m*sdv
      upper = mean_value + m*sdv
    end if
  end subroutine smean_sdl

  pure subroutine smedian_hilow(x, conf_int, median, lower, upper)
    real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted before quantile calculation.
    real(dp), intent(in), optional :: conf_int !! Central probability interval width in (0,1); defaults to 0.95.
    real(dp), intent(out) :: median !! Type-7 sample median.
    real(dp), intent(out) :: lower !! Type-7 lower central quantile.
    real(dp), intent(out) :: upper !! Type-7 upper central quantile.
    real(dp), allocatable :: xx(:)
    real(dp) :: ci
    integer :: i, j, n
    ci = 0.95_dp
    if (present(conf_int)) ci = conf_int
    n = count(.not. ieee_is_nan(x))
    if (n == 0 .or. ci <= 0.0_dp .or. ci >= 1.0_dp) then
      median = nan_dp()
      lower = nan_dp()
      upper = nan_dp()
      return
    end if
    allocate(xx(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) then
        j = j + 1
        xx(j) = x(i)
      end if
    end do
    call sort_real(xx)
    median = quantile_type7_sorted(xx, 0.5_dp)
    lower = quantile_type7_sorted(xx, (1.0_dp-ci)/2.0_dp)
    upper = quantile_type7_sorted(xx, (1.0_dp+ci)/2.0_dp)
  end subroutine smedian_hilow

  pure subroutine approx_extrap(x, y, xout, yout)
    real(dp), intent(in) :: x(:) !! Interpolation abscissae; NaN pairs are omitted and duplicate x values keep the first pair.
    real(dp), intent(in) :: y(:) !! Numeric ordinates paired elementwise with x.
    real(dp), intent(in) :: xout(:) !! Query points; NaNs return NaN.
    real(dp), intent(out) :: yout(:) !! Linear interpolation within the range and linear extrapolation beyond it.
    real(dp), allocatable :: xx(:), yy(:)
    integer :: i, k, n
    yout = nan_dp()
    if (size(x) /= size(y) .or. size(yout) /= size(xout)) return
    call xy_sort_no_dup_no_na(x, y, xx, yy)
    n = size(xx)
    if (n < 2) return
    do k = 1, size(xout)
      if (ieee_is_nan(xout(k))) cycle
      if (xout(k) <= xx(1)) then
        yout(k) = yy(1) + (yy(2)-yy(1))/(xx(2)-xx(1))*(xout(k)-xx(1))
      else if (xout(k) >= xx(n)) then
        yout(k) = yy(n-1) + (yy(n)-yy(n-1))/(xx(n)-xx(n-1))*(xout(k)-xx(n-1))
      else
        do i = 1, n-1
          if (xout(k) >= xx(i) .and. xout(k) <= xx(i+1)) then
            yout(k) = yy(i) + (yy(i+1)-yy(i))/(xx(i+1)-xx(i))*(xout(k)-xx(i))
            exit
          end if
        end do
      end if
    end do
  end subroutine approx_extrap

  pure subroutine james_stein(y, group, group_ids, n_group, mean_group, shrunk_mean, shrink)
    real(dp), intent(in) :: y(:) !! Numeric responses; NaN responses are omitted.
    integer, intent(in) :: group(:) !! Integer group labels paired with y.
    integer, allocatable, intent(out) :: group_ids(:) !! Sorted distinct group labels represented in complete observations.
    integer, allocatable, intent(out) :: n_group(:) !! Complete-observation count for each returned group.
    real(dp), allocatable, intent(out) :: mean_group(:) !! Raw mean for each returned group.
    real(dp), allocatable, intent(out) :: shrunk_mean(:) !! James-Stein shrunk mean for each returned group.
    real(dp), allocatable, intent(out) :: shrink(:) !! Shrinkage multiplier constrained to [0,1].
    integer, allocatable :: ids(:)
    real(dp), allocatable :: ss(:), vmean(:)
    real(dp) :: grand, ssb
    integer :: i, j, k, ng, n

    if (size(y) /= size(group)) then
      allocate(group_ids(0), n_group(0), mean_group(0), shrunk_mean(0), shrink(0))
      return
    end if
    n = count(.not. ieee_is_nan(y))
    if (n == 0) then
      allocate(group_ids(0), n_group(0), mean_group(0), shrunk_mean(0), shrink(0))
      return
    end if
    allocate(ids(n))
    j = 0
    do i = 1, size(y)
      if (.not. ieee_is_nan(y(i))) then
        j = j + 1
        ids(j) = group(i)
      end if
    end do
    call sort_int(ids)
    ng = 1
    do i = 2, n
      if (ids(i) /= ids(i-1)) ng = ng + 1
    end do
    allocate(group_ids(ng), n_group(ng), mean_group(ng), shrunk_mean(ng), shrink(ng), ss(ng), vmean(ng))
    j = 1
    group_ids(1) = ids(1)
    do i = 2, n
      if (ids(i) /= ids(i-1)) then
        j = j + 1
        group_ids(j) = ids(i)
      end if
    end do
    n_group = 0
    mean_group = 0.0_dp
    do i = 1, size(y)
      if (ieee_is_nan(y(i))) cycle
      do k = 1, ng
        if (group(i) == group_ids(k)) then
          n_group(k) = n_group(k) + 1
          mean_group(k) = mean_group(k) + y(i)
          exit
        end if
      end do
    end do
    do k = 1, ng
      mean_group(k) = mean_group(k)/real(n_group(k),dp)
    end do
    grand = sum(pack(y, .not. ieee_is_nan(y)))/real(n,dp)
    if (ng < 3) then
      shrink = nan_dp()
      shrunk_mean = nan_dp()
      return
    end if
    ss = 0.0_dp
    do i = 1, size(y)
      if (ieee_is_nan(y(i))) cycle
      do k = 1, ng
        if (group(i) == group_ids(k)) then
          ss(k) = ss(k) + (y(i)-mean_group(k))**2
          exit
        end if
      end do
    end do
    do k = 1, ng
      if (n_group(k) > 1) then
        vmean(k) = ss(k)/(real(n_group(k),dp)*real(n_group(k)-1,dp))
      else
        vmean(k) = 0.0_dp
      end if
    end do
    ssb = sum((mean_group-sum(mean_group)/real(ng,dp))**2)
    if (ssb == 0.0_dp) then
      shrink = 0.0_dp
    else
      shrink = 1.0_dp - real(ng-3,dp)*vmean/ssb
      do k = 1, ng
        if (n_group(k) == 1) shrink(k) = 0.0_dp
        shrink(k) = min(1.0_dp, max(0.0_dp, shrink(k)))
      end do
    end if
    shrunk_mean = grand*(1.0_dp-shrink) + shrink*mean_group
  end subroutine james_stein

  pure real(dp) function quantile_type7_sorted(x, p) result(q)
    real(dp), intent(in) :: x(:) !! Sorted numeric sample with at least one element.
    real(dp), intent(in) :: p !! Quantile probability in [0,1].
    real(dp) :: h, g
    integer :: j, n
    n = size(x)
    if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      q = nan_dp()
      return
    end if
    if (n == 1) then
      q = x(1)
      return
    end if
    h = 1.0_dp + real(n-1,dp)*p
    j = int(floor(h))
    g = h-real(j,dp)
    if (j >= n) then
      q = x(n)
    else
      q = (1.0_dp-g)*x(j) + g*x(j+1)
    end if
  end function quantile_type7_sorted

  pure subroutine sort_int(x)
    integer, intent(inout) :: x(:) !! Integer vector sorted in ascending order in place.
    integer :: i, j, key
    do i = 2, size(x)
      key = x(i)
      j = i-1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j-1
      end do
      x(j+1) = key
    end do
  end subroutine sort_int

  pure subroutine sort_index_real(x, ord)
    real(dp), intent(in) :: x(:) !! Numeric values whose ascending stable ordering is requested.
    integer, intent(out) :: ord(:) !! One-based permutation indices placing x in ascending order.
    integer :: i, j, key
    if (size(ord) /= size(x)) return
    do i = 1, size(x)
      ord(i) = i
    end do
    do i = 2, size(x)
      key = ord(i)
      j = i - 1
      do while (j >= 1)
        if (x(ord(j)) <= x(key)) exit
        ord(j+1) = ord(j)
        j = j - 1
      end do
      ord(j+1) = key
    end do
  end subroutine sort_index_real

  pure real(dp) function pseudomedian(x) result(pm)
    real(dp), intent(in) :: x(:) !! Numeric observations; NaNs are not permitted.
    real(dp), allocatable :: mids(:)
    integer :: i,j,k,n,np
    n=size(x)
    if(n==0 .or. any(ieee_is_nan(x))) then
 pm=nan_dp()
 return
 end if
    np=n*(n+1)/2
    allocate(mids(np))
 k=0
    do i=1,n
      do j=i,n
        k=k+1
 mids(k)=(x(i)+x(j))/2.0_dp
      end do
    end do
    call sort_real(mids)
    if(mod(np,2)==1) then
      pm=mids((np+1)/2)
    else
      pm=(mids(np/2)+mids(np/2+1))/2.0_dp
    end if
  end function pseudomedian


  pure real(dp) function ftupwr(p1, p2, bign, r, alpha) result(power)
    real(dp), intent(in) :: p1 !! Event probability in group 1; intended range is [0,1].
    real(dp), intent(in) :: p2 !! Event probability in group 2; intended range is [0,1].
    real(dp), intent(in) :: bign !! Total sample size across both groups; must be positive.
    real(dp), intent(in) :: r !! Group-2 to group-1 sample-size ratio; upstream approximation is recommended for roughly 0.33 to 3.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability; must lie strictly between zero and one.
    real(dp) :: mstar, delta, rp1, zalp, pbar, qbar, term, num, den, zbet

    if (p1 < 0.0_dp .or. p1 > 1.0_dp .or. p2 < 0.0_dp .or. p2 > 1.0_dp .or. &
        bign <= 0.0_dp .or. r <= 0.0_dp .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      power = nan_dp()
      return
    end if
    mstar = bign / (r + 1.0_dp)
    delta = abs(p2 - p1)
    rp1 = r + 1.0_dp
    zalp = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    pbar = (p1 + r * p2) / rp1
    qbar = 1.0_dp - pbar
    term = r * delta**2 * mstar - rp1 * delta
    den = sqrt(r * p1 * (1.0_dp - p1) + p2 * (1.0_dp - p2))
    if (term < 0.0_dp .or. den <= 0.0_dp) then
      power = nan_dp()
      return
    end if
    num = sqrt(term) - zalp * sqrt(rp1 * pbar * qbar)
    zbet = num / den
    power = norm_cdf(zbet)
  end function ftupwr

  pure subroutine ftuss(p1, p2, r, alpha, beta, n1, n2)
    real(dp), intent(in) :: p1 !! Event probability in group 1; intended range is [0,1].
    real(dp), intent(in) :: p2 !! Event probability in group 2; intended range is [0,1].
    real(dp), intent(in) :: r !! Planned group-2 to group-1 sample-size ratio; must be positive.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability; must lie strictly between zero and one.
    real(dp), intent(in) :: beta !! Type-II error probability; target power is one minus beta and beta must be in (0,1).
    integer, intent(out) :: n1 !! Integer sample-size recommendation for group 1 using the upstream floor(m+1) convention.
    integer, intent(out) :: n2 !! Integer sample-size recommendation for group 2 using the upstream floor(r*m+1) convention.
    real(dp) :: zalp, zbet, rp1, pbar, qbar, q1, q2, delta, num, den, mp, m

    n1 = -1
    n2 = -1
    if (p1 < 0.0_dp .or. p1 > 1.0_dp .or. p2 < 0.0_dp .or. p2 > 1.0_dp .or. &
        p1 == p2 .or. r <= 0.0_dp .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. &
        beta <= 0.0_dp .or. beta >= 1.0_dp) return
    zalp = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    zbet = inv_norm_cdf(1.0_dp - beta)
    rp1 = r + 1.0_dp
    pbar = (p1 + r * p2) / rp1
    qbar = 1.0_dp - pbar
    q1 = 1.0_dp - p1
    q2 = 1.0_dp - p2
    delta = abs(p2 - p1)
    num = (zalp * sqrt(rp1 * pbar * qbar) + zbet * sqrt(r * p1 * q1 + p2 * q2))**2
    den = r * delta**2
    mp = num / den
    m = 0.25_dp * mp * (1.0_dp + sqrt(1.0_dp + 2.0_dp * rp1 / (r * mp * delta)))**2
    n1 = floor(m + 1.0_dp)
    n2 = floor(m * r + 1.0_dp)
  end subroutine ftuss

  pure subroutine ecdf_steps(x, xout, yout, extend, do_extend)
    real(dp), intent(in) :: x(:) !! Sample values; NaNs are ignored as in the R implementation.
    real(dp), allocatable, intent(out) :: xout(:) !! Sorted unique sample values, optionally bracketed by extension points.
    real(dp), allocatable, intent(out) :: yout(:) !! Empirical CDF evaluated at each returned x coordinate.
    real(dp), optional, intent(in) :: extend(:) !! Optional lower/upper coordinates overriding default range extension.
    logical, optional, intent(in) :: do_extend !! Whether to extend the range; defaults to true as in the R interface.
    real(dp), allocatable :: clean(:), uniq(:), work(:)
    real(dp) :: eps, lo, hi
    logical :: ext
    integer :: i, j, n, nu

    n = count(.not. ieee_is_nan(x))
    if (n == 0) then
      allocate(xout(0), yout(0))
      return
    end if
    allocate(clean(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) then
        j = j + 1
        clean(j) = x(i)
      end if
    end do
    call sort_real(clean)
    allocate(work(n))
    nu = 0
    nu = 1
    work(1) = clean(1)
    do i = 2, n
      if (clean(i) /= clean(i - 1)) then
        nu = nu + 1
        work(nu) = clean(i)
      end if
    end do
    allocate(uniq(nu))
    uniq = work(1:nu)
    ext = .true.
    if (present(do_extend)) ext = do_extend
    if (.not. ext) then
      allocate(xout(nu), yout(nu))
      xout = uniq
    else
      if (present(extend)) then
        if (size(extend) /= 2) then
          allocate(xout(0), yout(0))
          return
        end if
        lo = extend(1)
        hi = extend(2)
      else
        lo = uniq(1)
        hi = uniq(nu)
        eps = (hi - lo) / 20.0_dp
        lo = lo - eps
        hi = hi + eps
      end if
      allocate(xout(nu + 2), yout(nu + 2))
      xout(1) = lo
      xout(2:nu + 1) = uniq
      xout(nu + 2) = hi
    end if
    do i = 1, size(xout)
      yout(i) = real(count(clean <= xout(i)), dp) / real(n, dp)
    end do
  end subroutine ecdf_steps

  pure subroutine lag_numeric(x, shift, out)
    real(dp), intent(in) :: x(:) !! Numeric vector to lag; NaNs in x are preserved when shifted.
    integer, intent(in) :: shift !! Number of observations to shift; positive values lag and negative values lead.
    real(dp), intent(out) :: out(:) !! Shifted vector, same size as x, padded with quiet NaNs.
    integer :: n

    n = size(x)
    if (size(out) /= n) return
    out = nan_dp()
    if (shift == 0) then
      out = x
    else if (abs(shift) < n) then
      if (shift > 0) then
        out(shift + 1:n) = x(1:n - shift)
      else
        out(1:n + shift) = x(1 - shift:n)
      end if
    end if
  end subroutine lag_numeric

  pure subroutine nomiss_vector(x, out)
    real(dp), intent(in) :: x(:) !! Numeric vector from which quiet-NaN elements are removed.
    real(dp), allocatable, intent(out) :: out(:) !! Vector containing only non-NaN elements in original order.
    integer :: i, j, n

    n = count(.not. ieee_is_nan(x))
    allocate(out(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) then
        j = j + 1
        out(j) = x(i)
      end if
    end do
  end subroutine nomiss_vector

  pure subroutine nomiss_matrix(x, out)
    real(dp), intent(in) :: x(:,:) !! Numeric matrix; any row containing a quiet NaN is omitted.
    real(dp), allocatable, intent(out) :: out(:,:) !! Matrix of complete rows, preserving column count and row order.
    integer :: i, j, nr

    nr = 0
    do i = 1, size(x, 1)
      if (.not. any(ieee_is_nan(x(i, :)))) nr = nr + 1
    end do
    allocate(out(nr, size(x, 2)))
    j = 0
    do i = 1, size(x, 1)
      if (.not. any(ieee_is_nan(x(i, :)))) then
        j = j + 1
        out(j, :) = x(i, :)
      end if
    end do
  end subroutine nomiss_matrix

  pure subroutine fillin_scalar(v, p, out)
    real(dp), intent(in) :: v(:) !! Numeric vector whose quiet-NaN entries are to be replaced.
    real(dp), intent(in) :: p !! Scalar replacement value used for every NaN entry.
    real(dp), intent(out) :: out(:) !! Filled vector with the same shape as v.
    integer :: i

    if (size(out) /= size(v)) return
    out = v
    do i = 1, size(v)
      if (ieee_is_nan(v(i))) out(i) = p
    end do
  end subroutine fillin_scalar

  pure subroutine fillin_vector(v, p, out)
    real(dp), intent(in) :: v(:) !! Numeric vector whose quiet-NaN entries are to be replaced.
    real(dp), intent(in) :: p(:) !! Replacement values recycled by position using R-style vector recycling.
    real(dp), intent(out) :: out(:) !! Filled vector with the same shape as v.
    integer :: i

    if (size(out) /= size(v) .or. size(p) == 0) return
    out = v
    do i = 1, size(v)
      if (ieee_is_nan(v(i))) out(i) = p(mod(i - 1, size(p)) + 1)
    end do
  end subroutine fillin_vector

  pure subroutine cumcategory(y, nlevels, out)
    integer, intent(in) :: y(:) !! Integer category codes corresponding to ordered R factor levels, normally 1 through nlevels.
    integer, intent(in) :: nlevels !! Number of ordered levels; output has nlevels minus one columns.
    integer, allocatable, intent(out) :: out(:,:) !! Indicator matrix whose column j is one when y is at least level j+1.
    integer :: i, j

    if (nlevels < 2) then
      allocate(out(size(y), 0))
      return
    end if
    allocate(out(size(y), nlevels - 1))
    do j = 2, nlevels
      do i = 1, size(y)
        out(i, j - 1) = merge(1, 0, y(i) >= j)
      end do
    end do
  end subroutine cumcategory

  pure subroutine midrank(x,r)
    real(dp), intent(in) :: x(:) !! Values to rank with average ranks for ties.
    real(dp), intent(out) :: r(:) !! Midranks corresponding to original positions.
    integer, allocatable :: idx(:)
    integer :: i,j,n,k,t
    real(dp) :: rr
    n=size(x)
 allocate(idx(n))
 idx=[(i,i=1,n)]
    do i=2,n
      t=idx(i)
 j=i-1
      do while (j >= 1)
        if (x(idx(j)) <= x(t)) exit
        idx(j+1) = idx(j)
        j = j - 1
      end do
      idx(j+1)=t
    end do
    i=1
    do while(i<=n)
      j=i
      do while (j < n)
        if (x(idx(j+1)) /= x(idx(i))) exit
        j = j + 1
      end do
      rr=(real(i,dp)+real(j,dp))/2.0_dp
      do k=i,j
 r(idx(k))=rr
 end do
      i=j+1
    end do
  end subroutine midrank

  pure subroutine sort_real(x)
    real(dp), intent(inout) :: x(:) !! Real vector sorted in ascending order in place.
    integer :: i,j
    real(dp) :: key
    do i=2,size(x)
      key=x(i)
 j=i-1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j - 1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  pure subroutine sort_pairs(x,w)
    real(dp), intent(inout) :: x(:) !! Values sorted ascending in place.
    real(dp), intent(inout) :: w(:) !! Companion weights permuted with x.
    integer :: i,j
    real(dp) :: kx,kw
    do i=2,size(x)
      kx=x(i)
 kw=w(i)
 j=i-1
      do while (j >= 1)
        if (x(j) <= kx) exit
        x(j+1) = x(j)
        w(j+1) = w(j)
        j = j - 1
      end do
      x(j+1)=kx
 w(j+1)=kw
    end do
  end subroutine sort_pairs


  pure subroutine deff(y, cluster, nobs, nclusters, rho, design_effect)
    real(dp), intent(in) :: y(:) !! Numeric response values; NaNs are omitted together with their cluster labels.
    integer, intent(in) :: cluster(:) !! Integer cluster labels corresponding one-to-one with y; labels need not be consecutive.
    integer, intent(out) :: nobs !! Number of complete observations used in the calculation.
    integer, intent(out) :: nclusters !! Number of distinct cluster labels among complete observations.
    real(dp), intent(out) :: rho !! Estimated intracluster correlation implied by the one-way ANOVA decomposition.
    real(dp), intent(out) :: design_effect !! Estimated variance inflation factor 1 + (B - 1) * rho.
    integer, allocatable :: ids(:), counts(:)
    real(dp), allocatable :: sums(:), means(:)
    real(dp) :: overall, sst, sse, r2, fstat, g, b
    integer :: i, j, k

    nobs = 0
    nclusters = 0
    rho = nan_dp()
    design_effect = nan_dp()
    if (size(cluster) /= size(y)) return
    if (size(y) == 0) return

    allocate(ids(size(y)), counts(size(y)), sums(size(y)), means(size(y)))
    ids = 0
    counts = 0
    sums = 0.0_dp
    means = 0.0_dp

    do i = 1, size(y)
      if (ieee_is_nan(y(i))) cycle
      nobs = nobs + 1
      j = 0
      do k = 1, nclusters
        if (ids(k) == cluster(i)) then
          j = k
          exit
        end if
      end do
      if (j == 0) then
        nclusters = nclusters + 1
        j = nclusters
        ids(j) = cluster(i)
      end if
      counts(j) = counts(j) + 1
      sums(j) = sums(j) + y(i)
    end do
    if (nobs < 2 .or. nclusters < 1 .or. nobs <= nclusters) return

    overall = 0.0_dp
    do j = 1, nclusters
      means(j) = sums(j) / real(counts(j), dp)
      overall = overall + sums(j)
    end do
    overall = overall / real(nobs, dp)

    sst = 0.0_dp
    sse = 0.0_dp
    do i = 1, size(y)
      if (ieee_is_nan(y(i))) cycle
      sst = sst + (y(i) - overall)**2
      do j = 1, nclusters
        if (ids(j) == cluster(i)) then
          sse = sse + (y(i) - means(j))**2
          exit
        end if
      end do
    end do
    if (sst <= 0.0_dp) return

    r2 = 1.0_dp - sse / sst
    if (r2 >= 1.0_dp) then
      rho = 1.0_dp
    else
      fstat = r2 * real(nobs - nclusters, dp) / ((1.0_dp - r2) * real(nclusters, dp))
      g = (fstat - 1.0_dp) * real(nclusters, dp) / real(nobs, dp)
      rho = g / (1.0_dp + g)
    end if
    b = sum(real(counts(1:nclusters), dp)**2) / real(nobs, dp)
    design_effect = 1.0_dp + (b - 1.0_dp) * rho
  end subroutine deff

  pure subroutine r2_measures(lr, p, n, padj, r2, r2adj, ess, frequencies, r2_ess, r2adj_ess)
    real(dp), intent(in) :: lr !! Likelihood-ratio chi-square statistic, normally nonnegative.
    integer, intent(in) :: p !! Number of fitted non-intercept regression parameters.
    real(dp), intent(in) :: n !! Raw sample size, strictly positive.
    integer, intent(in) :: padj !! Adjustment selector: 1 subtracts p from lr; 2 uses the classical adjusted-R2 multiplier.
    real(dp), intent(out) :: r2 !! Maddala-Cox-Snell R-squared based on n.
    real(dp), intent(out) :: r2adj !! Parameter-adjusted R-squared based on n.
    real(dp), intent(in), optional :: ess !! Optional directly supplied effective sample size.
    real(dp), intent(in), optional :: frequencies(:) !! Outcome frequencies used to derive effective sample size when ess is absent.
    real(dp), intent(out), optional :: r2_ess !! R-squared using effective sample size; NaN when it equals the raw size.
    real(dp), intent(out), optional :: r2adj_ess !! Adjusted R-squared using effective sample size; NaN when sizes are equal.
    real(dp) :: neff, totalf

    r2 = nan_dp()
    r2adj = nan_dp()
    if (present(r2_ess)) r2_ess = nan_dp()
    if (present(r2adj_ess)) r2adj_ess = nan_dp()
    if (n <= 0.0_dp .or. p < 0) return

    r2 = 1.0_dp - exp(-lr / n)
    if (padj == 1) then
      r2adj = 1.0_dp - exp(-max(lr - real(p, dp), 0.0_dp) / n)
    else if (padj == 2) then
      if (n <= real(p + 1, dp)) return
      r2adj = 1.0_dp - exp(-lr / n) * (n - 1.0_dp) / (n - real(p, dp) - 1.0_dp)
    else
      r2adj = nan_dp()
      return
    end if

    if (.not. present(r2_ess) .and. .not. present(r2adj_ess)) return
    if (present(ess)) then
      neff = ess
    else if (present(frequencies)) then
      if (size(frequencies) == 0 .or. any(frequencies < 0.0_dp)) return
      totalf = sum(frequencies)
      if (totalf <= 0.0_dp) return
      neff = n * (1.0_dp - sum((frequencies / totalf)**3))
    else
      return
    end if
    if (neff <= 0.0_dp) return
    if (abs(neff - n) <= epsilon(n) * max(1.0_dp, abs(n))) return

    if (present(r2_ess)) r2_ess = 1.0_dp - exp(-lr / neff)
    if (present(r2adj_ess)) then
      if (padj == 1) then
        r2adj_ess = 1.0_dp - exp(-max(lr - real(p, dp), 0.0_dp) / neff)
      else if (neff > real(p + 1, dp)) then
        r2adj_ess = 1.0_dp - exp(-lr / neff) * (neff - 1.0_dp) / (neff - real(p, dp) - 1.0_dp)
      end if
    end if
  end subroutine r2_measures

  pure subroutine binconf_wilson(x, n, alpha, estimate, lower, upper)
    real(dp), intent(in) :: x !! Number of successes, between zero and n inclusive.
    real(dp), intent(in) :: n !! Number of Bernoulli trials, strictly positive.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability in (0,1).
    real(dp), intent(out) :: estimate !! Observed success proportion x/n.
    real(dp), intent(out) :: lower !! Wilson lower confidence limit, including Hmisc's x=1 small-count correction.
    real(dp), intent(out) :: upper !! Wilson upper confidence limit, including Hmisc's x=n-1 small-count correction.
    real(dp) :: zcrit, z2, p0, half

    estimate = nan_dp()
    lower = nan_dp()
    upper = nan_dp()
    if (n <= 0.0_dp .or. x < 0.0_dp .or. x > n .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) return
    zcrit = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    z2 = zcrit * zcrit
    p0 = x / n
    estimate = p0
    half = zcrit * sqrt((p0 * (1.0_dp - p0) + z2 / (4.0_dp * n)) / n)
    lower = (p0 + z2 / (2.0_dp * n) - half) / (1.0_dp + z2 / n)
    upper = (p0 + z2 / (2.0_dp * n) + half) / (1.0_dp + z2 / n)
    if (x == 1.0_dp) lower = -log(1.0_dp - alpha) / n
    if (x == n - 1.0_dp) upper = 1.0_dp + log(1.0_dp - alpha) / n
  end subroutine binconf_wilson

  pure subroutine binconf_asymptotic(x, n, alpha, estimate, lower, upper)
    real(dp), intent(in) :: x !! Number of successes, between zero and n inclusive.
    real(dp), intent(in) :: n !! Number of Bernoulli trials, strictly positive.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability in (0,1).
    real(dp), intent(out) :: estimate !! Observed success proportion x/n.
    real(dp), intent(out) :: lower !! Unclipped normal-approximation lower confidence limit used by Hmisc.
    real(dp), intent(out) :: upper !! Unclipped normal-approximation upper confidence limit used by Hmisc.
    real(dp) :: zcrit, se

    estimate = nan_dp()
    lower = nan_dp()
    upper = nan_dp()
    if (n <= 0.0_dp .or. x < 0.0_dp .or. x > n .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) return
    estimate = x / n
    zcrit = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)
    se = sqrt(estimate * (1.0_dp - estimate) / n)
    lower = estimate - zcrit * se
    upper = estimate + zcrit * se
  end subroutine binconf_asymptotic

  pure subroutine invert_tabulated_linear(x, y, aty, xinv)
    real(dp), intent(in) :: x(:) !! Tabulated x values corresponding one-to-one with y.
    real(dp), intent(in) :: y(:) !! Tabulated transformed values to invert; duplicate y values are averaged in x.
    real(dp), intent(in) :: aty(:) !! Target y values at which the inverse mapping is requested.
    real(dp), allocatable, intent(out) :: xinv(:) !! Interpolated inverse x values; endpoints are clamped as in approx(rule=2).
    real(dp), allocatable :: ys(:), xs(:), uy(:), ux(:)
    real(dp) :: sy, sx
    integer :: i, j, m, nu

    allocate(xinv(size(aty)))
    xinv = nan_dp()
    if (size(x) /= size(y) .or. size(x) == 0) return
    m = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (m == 0) return
    allocate(xs(m), ys(m))
    j = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
      j = j + 1
      xs(j) = x(i)
      ys(j) = y(i)
    end do
    call sort_pairs(ys, xs)
    allocate(uy(m), ux(m))
    nu = 0
    i = 1
    do while (i <= m)
      sy = ys(i)
      sx = 0.0_dp
      j = 0
      do while (i <= m)
        if (ys(i) /= sy) exit
        sx = sx + xs(i)
        j = j + 1
        i = i + 1
      end do
      nu = nu + 1
      uy(nu) = sy
      ux(nu) = sx / real(j, dp)
    end do
    if (nu == 1) then
      where (.not. ieee_is_nan(aty)) xinv = ux(1)
      return
    end if
    do i = 1, size(aty)
      if (ieee_is_nan(aty(i))) cycle
      if (aty(i) <= uy(1)) then
        xinv(i) = ux(1)
      else if (aty(i) >= uy(nu)) then
        xinv(i) = ux(nu)
      else
        do j = 1, nu - 1
          if (aty(i) >= uy(j) .and. aty(i) <= uy(j + 1)) then
            xinv(i) = ux(j) + (ux(j + 1) - ux(j)) * (aty(i) - uy(j)) / (uy(j + 1) - uy(j))
            exit
          end if
        end do
      end if
    end do
  end subroutine invert_tabulated_linear

  pure subroutine score_binary_max(x, points, scores)
    real(dp), intent(in) :: x(:,:) !! Binary/numeric item matrix with observations in rows and items in columns; NaNs are ignored.
    real(dp), intent(in) :: points(:) !! Per-item numeric scores, one for each column of x.
    real(dp), intent(out) :: scores(:) !! Row-wise maximum of x(:,j)*points(j), corresponding to score.binary(fun=max, na.rm=TRUE).
    integer :: i, j
    real(dp) :: v

    scores = nan_dp()
    if (size(points) /= size(x,2) .or. size(scores) /= size(x,1)) return
    do i = 1, size(x,1)
      do j = 1, size(x,2)
        if (ieee_is_nan(x(i,j))) cycle
        v = x(i,j) * points(j)
        if (ieee_is_nan(scores(i))) then
          scores(i) = v
        else
          scores(i) = max(scores(i), v)
        end if
      end do
    end do
  end subroutine score_binary_max

  pure subroutine score_binary_sum(x, points, na_rm, scores)
    real(dp), intent(in) :: x(:,:) !! Binary/numeric item matrix with observations in rows and items in columns.
    real(dp), intent(in) :: points(:) !! Per-item numeric scores, one for each column of x.
    logical, intent(in) :: na_rm !! If true, omit NaNs; if false, any NaN in a row yields a NaN row score.
    real(dp), intent(out) :: scores(:) !! Row-wise sum of x(:,j)*points(j), with missing-value behavior controlled by na_rm.
    integer :: i, j
    logical :: any_value

    scores = nan_dp()
    if (size(points) /= size(x,2) .or. size(scores) /= size(x,1)) return
    do i = 1, size(x,1)
      scores(i) = 0.0_dp
      any_value = .false.
      do j = 1, size(x,2)
        if (ieee_is_nan(x(i,j))) then
          if (.not. na_rm) then
            scores(i) = nan_dp()
            exit
          end if
          cycle
        end if
        scores(i) = scores(i) + x(i,j) * points(j)
        any_value = .true.
      end do
      if (na_rm .and. .not. any_value) scores(i) = 0.0_dp
    end do
  end subroutine score_binary_sum


  pure subroutine mhgr(y, group, strata, conf_int, rr, lower, upper, group_counts)
    real(dp), intent(in) :: y(:) !! Binary or nonnegative event values; NaNs are omitted before stratified aggregation.
    integer, intent(in) :: group(:) !! Treatment-group codes, restricted to 1 or 2 for included observations.
    integer, intent(in) :: strata(:) !! Positive integer stratum identifiers; nonpositive values are omitted.
    real(dp), intent(in) :: conf_int !! Two-sided confidence level strictly between zero and one.
    real(dp), intent(out) :: rr !! Mantel-Haenszel common risk ratio comparing group 1 with group 2.
    real(dp), intent(out) :: lower !! Greenland-Robins lower confidence limit for the common risk ratio.
    real(dp), intent(out) :: upper !! Greenland-Robins upper confidence limit for the common risk ratio.
    integer, intent(out) :: group_counts(2) !! Numbers of included observations in groups 1 and 2 after omission.
    integer, allocatable :: n1(:), n2(:)
    real(dp), allocatable :: x1(:), x2(:)
    integer :: i, ns, st
    real(dp) :: nn, tk, rsum, ssum, dsum, sigma, zcrit, dterm

    rr = nan_dp()
    lower = nan_dp()
    upper = nan_dp()
    group_counts = 0
    if (size(group) /= size(y) .or. size(strata) /= size(y)) return
    if (conf_int <= 0.0_dp .or. conf_int >= 1.0_dp) return
    if (size(y) == 0) return
    ns = maxval(strata, mask=strata > 0)
    if (ns <= 0) return
    allocate(n1(ns), n2(ns), x1(ns), x2(ns))
    n1 = 0
    n2 = 0
    x1 = 0.0_dp
    x2 = 0.0_dp
    do i = 1, size(y)
      if (ieee_is_nan(y(i))) cycle
      if (strata(i) <= 0) cycle
      if (group(i) == 1) then
        n1(strata(i)) = n1(strata(i)) + 1
        x1(strata(i)) = x1(strata(i)) + y(i)
        group_counts(1) = group_counts(1) + 1
      else if (group(i) == 2) then
        n2(strata(i)) = n2(strata(i)) + 1
        x2(strata(i)) = x2(strata(i)) + y(i)
        group_counts(2) = group_counts(2) + 1
      else
        return
      end if
    end do
    if (any(group_counts == 0)) return
    rsum = 0.0_dp
    ssum = 0.0_dp
    dsum = 0.0_dp
    do st = 1, ns
      nn = real(n1(st) + n2(st), dp)
      if (nn <= 0.0_dp) cycle
      tk = x1(st) + x2(st)
      rsum = rsum + x1(st) * real(n2(st), dp) / nn
      ssum = ssum + x2(st) * real(n1(st), dp) / nn
      dterm = (real(n1(st) * n2(st), dp) * tk - x1(st) * x2(st) * nn) / (nn * nn)
      dsum = dsum + dterm
    end do
    if (rsum <= 0.0_dp .or. ssum <= 0.0_dp .or. dsum < 0.0_dp) return
    rr = rsum / ssum
    sigma = sqrt(dsum / (rsum * ssum))
    zcrit = inv_norm_cdf((1.0_dp + conf_int) / 2.0_dp)
    lower = rr * exp(-zcrit * sigma)
    upper = rr * exp(zcrit * sigma)
  end subroutine mhgr

  pure subroutine lrcum(a, b, c, d, conf_int, lrpos, lower_pos, upper_pos, lrneg, lower_neg, upper_neg, &
                        lrpos_cum, lower_pos_cum, upper_pos_cum, lrneg_cum, lower_neg_cum, upper_neg_cum)
    real(dp), intent(in) :: a(:) !! True-positive counts for ordered diagnostic thresholds.
    real(dp), intent(in) :: b(:) !! False-positive counts for ordered diagnostic thresholds.
    real(dp), intent(in) :: c(:) !! False-negative counts for ordered diagnostic thresholds.
    real(dp), intent(in) :: d(:) !! True-negative counts for ordered diagnostic thresholds.
    real(dp), intent(in) :: conf_int !! Two-sided confidence level strictly between zero and one.
    real(dp), intent(out) :: lrpos(:) !! Positive likelihood ratio at each threshold.
    real(dp), intent(out) :: lower_pos(:) !! Lower confidence limits for threshold-specific positive likelihood ratios.
    real(dp), intent(out) :: upper_pos(:) !! Upper confidence limits for threshold-specific positive likelihood ratios.
    real(dp), intent(out) :: lrneg(:) !! Negative likelihood ratio at each threshold.
    real(dp), intent(out) :: lower_neg(:) !! Lower confidence limits for threshold-specific negative likelihood ratios.
    real(dp), intent(out) :: upper_neg(:) !! Upper confidence limits for threshold-specific negative likelihood ratios.
    real(dp), intent(out) :: lrpos_cum(:) !! Cumulative products of positive likelihood ratios.
    real(dp), intent(out) :: lower_pos_cum(:) !! Lower confidence limits for cumulative positive likelihood ratios.
    real(dp), intent(out) :: upper_pos_cum(:) !! Upper confidence limits for cumulative positive likelihood ratios.
    real(dp), intent(out) :: lrneg_cum(:) !! Cumulative products of negative likelihood ratios.
    real(dp), intent(out) :: lower_neg_cum(:) !! Lower confidence limits for cumulative negative likelihood ratios.
    real(dp), intent(out) :: upper_neg_cum(:) !! Upper confidence limits for cumulative negative likelihood ratios.
    real(dp), allocatable :: aa(:), bb(:), cc(:), dd(:), vp(:), vn(:)
    real(dp) :: zcrit, cvp, cvn
    integer :: i, n

    n = size(a)
    lrpos = nan_dp()
    lower_pos = nan_dp()
    upper_pos = nan_dp()
    lrneg = nan_dp()
    lower_neg = nan_dp()
    upper_neg = nan_dp()
    lrpos_cum = nan_dp()
    lower_pos_cum = nan_dp()
    upper_pos_cum = nan_dp()
    lrneg_cum = nan_dp()
    lower_neg_cum = nan_dp()
    upper_neg_cum = nan_dp()
    if (size(b) /= n .or. size(c) /= n .or. size(d) /= n) return
    if (size(lrpos) /= n .or. size(lower_pos) /= n .or. size(upper_pos) /= n) return
    if (size(lrneg) /= n .or. size(lower_neg) /= n .or. size(upper_neg) /= n) return
    if (size(lrpos_cum) /= n .or. size(lower_pos_cum) /= n .or. size(upper_pos_cum) /= n) return
    if (size(lrneg_cum) /= n .or. size(lower_neg_cum) /= n .or. size(upper_neg_cum) /= n) return
    if (n == 0 .or. conf_int <= 0.0_dp .or. conf_int >= 1.0_dp) return
    if (any(ieee_is_nan(a)) .or. any(ieee_is_nan(b)) .or. any(ieee_is_nan(c)) .or. any(ieee_is_nan(d))) return
    if (any(a < 0.0_dp) .or. any(b < 0.0_dp) .or. any(c < 0.0_dp) .or. any(d < 0.0_dp)) return
    aa = a
    bb = b
    cc = c
    dd = d
    if (any(aa == 0.0_dp) .or. any(bb == 0.0_dp) .or. any(cc == 0.0_dp) .or. any(dd == 0.0_dp)) then
      aa = aa + 0.5_dp
      bb = bb + 0.5_dp
      cc = cc + 0.5_dp
      dd = dd + 0.5_dp
    end if
    allocate(vp(n), vn(n))
    lrpos = (aa / (aa + cc)) / (bb / (bb + dd))
    lrneg = (cc / (aa + cc)) / (dd / (bb + dd))
    vp = 1.0_dp / aa - 1.0_dp / (aa + cc) + 1.0_dp / bb - 1.0_dp / (bb + dd)
    vn = 1.0_dp / dd - 1.0_dp / (bb + dd) + 1.0_dp / cc - 1.0_dp / (aa + cc)
    zcrit = inv_norm_cdf((1.0_dp + conf_int) / 2.0_dp)
    lower_pos = exp(log(lrpos) - zcrit * sqrt(vp))
    upper_pos = exp(log(lrpos) + zcrit * sqrt(vp))
    lower_neg = exp(log(lrneg) - zcrit * sqrt(vn))
    upper_neg = exp(log(lrneg) + zcrit * sqrt(vn))
    lrpos_cum(1) = lrpos(1)
    lrneg_cum(1) = lrneg(1)
    cvp = vp(1)
    cvn = vn(1)
    lower_pos_cum(1) = exp(log(lrpos_cum(1)) - zcrit * sqrt(cvp))
    upper_pos_cum(1) = exp(log(lrpos_cum(1)) + zcrit * sqrt(cvp))
    lower_neg_cum(1) = exp(log(lrneg_cum(1)) - zcrit * sqrt(cvn))
    upper_neg_cum(1) = exp(log(lrneg_cum(1)) + zcrit * sqrt(cvn))
    do i = 2, n
      lrpos_cum(i) = lrpos_cum(i - 1) * lrpos(i)
      lrneg_cum(i) = lrneg_cum(i - 1) * lrneg(i)
      cvp = cvp + vp(i)
      cvn = cvn + vn(i)
      lower_pos_cum(i) = exp(log(lrpos_cum(i)) - zcrit * sqrt(cvp))
      upper_pos_cum(i) = exp(log(lrpos_cum(i)) + zcrit * sqrt(cvp))
      lower_neg_cum(i) = exp(log(lrneg_cum(i)) - zcrit * sqrt(cvn))
      upper_neg_cum(i) = exp(log(lrneg_cum(i)) + zcrit * sqrt(cvn))
    end do
  end subroutine lrcum

  pure subroutine num_denom_setup(num, denom, subs, weights, y)
    real(dp), intent(in) :: num(:) !! Numerator counts; NaNs and entries paired with zero denominators are omitted.
    real(dp), intent(in) :: denom(:) !! Denominator counts corresponding one-to-one with num.
    integer, allocatable, intent(out) :: subs(:) !! One-based source indices for numerator and complement contributions.
    real(dp), allocatable, intent(out) :: weights(:) !! Numerator or denominator-minus-numerator weights for expanded records.
    integer, allocatable, intent(out) :: y(:) !! Expanded binary outcomes: one for numerator contributions and zero for complements.
    integer, allocatable :: valid(:)
    real(dp) :: other
    integer :: i, nv, m, k

    if (size(num) /= size(denom)) then
      allocate(subs(0), weights(0), y(0))
      return
    end if
    nv = count(.not. ieee_is_nan(num) .and. .not. ieee_is_nan(denom) .and. denom /= 0.0_dp)
    allocate(valid(nv))
    k = 0
    do i = 1, size(num)
      if (ieee_is_nan(num(i)) .or. ieee_is_nan(denom(i)) .or. denom(i) == 0.0_dp) cycle
      k = k + 1
      valid(k) = i
    end do
    m = 0
    do k = 1, nv
      i = valid(k)
      if (num(i) > 0.0_dp) m = m + 1
      other = denom(i) - num(i)
      if (other > 0.0_dp) m = m + 1
    end do
    allocate(subs(m), weights(m), y(m))
    k = 0
    do i = 1, nv
      if (num(valid(i)) > 0.0_dp) then
        k = k + 1
        subs(k) = valid(i)
        weights(k) = num(valid(i))
        y(k) = 1
      end if
    end do
    do i = 1, nv
      other = denom(valid(i)) - num(valid(i))
      if (other > 0.0_dp) then
        k = k + 1
        subs(k) = valid(i)
        weights(k) = other
        y(k) = 0
      end if
    end do
  end subroutine num_denom_setup

  pure subroutine rmultinom(probs, uniforms, draws, status)
    real(dp), intent(in) :: probs(:,:) !! Row-wise multinomial probabilities; each row must sum to one within numerical tolerance.
    real(dp), intent(in) :: uniforms(:,:) !! Uniform variates in [0,1), rows matching probs and columns repeated draws.
    integer, intent(out) :: draws(:,:) !! Simulated one-based category indices, with shape matching uniforms.
    integer, intent(out) :: status !! Zero on success; nonzero for dimension, probability, or uniform-input errors.
    real(dp) :: total, cumulative
    integer :: i, j, k, ncat

    draws = 0
    status = 0
    if (size(uniforms,1) /= size(probs,1) .or. any(shape(draws) /= shape(uniforms))) then
      status = 1
      return
    end if
    ncat = size(probs,2)
    if (ncat < 1) then
      status = 1
      return
    end if
    if (any(probs < 0.0_dp) .or. any(ieee_is_nan(probs))) then
      status = 2
      return
    end if
    if (any(uniforms < 0.0_dp) .or. any(uniforms >= 1.0_dp) .or. any(ieee_is_nan(uniforms))) then
      status = 3
      return
    end if
    do i = 1, size(probs,1)
      total = sum(probs(i,:))
      if (abs(total - 1.0_dp) > 1.0e-5_dp) then
        status = 2
        return
      end if
      do j = 1, size(uniforms,2)
        cumulative = 0.0_dp
        draws(i,j) = ncat
        do k = 1, ncat - 1
          cumulative = cumulative + probs(i,k)
          if (uniforms(i,j) <= cumulative) then
            draws(i,j) = k
            exit
          end if
        end do
      end do
    end do
  end subroutine rmultinom

  pure subroutine qrxcenter(x, x_transformed, r_transform, ri_transform, xbar, status)
    real(dp), intent(in) :: x(:,:) !! Observations by variables; centered columns must have full rank.
    real(dp), intent(out) :: x_transformed(:,:) !! Mean-centered orthogonalized columns scaled to sample standard deviation one.
    real(dp), intent(out) :: r_transform(:,:) !! Matrix mapping centered raw data to x_transformed by matrix multiplication.
    real(dp), intent(out) :: ri_transform(:,:) !! Inverse transform mapping x_transformed back to centered raw data.
    real(dp), intent(out) :: xbar(:) !! Arithmetic means of the original columns.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or rank deficiency.
    real(dp), allocatable :: a(:,:), qmat(:,:), rr(:,:), v(:), invrr(:,:)
    real(dp) :: rij, normv, sn, tol
    integer :: i, j, n, p

    x_transformed = nan_dp()
    r_transform = nan_dp()
    ri_transform = nan_dp()
    xbar = nan_dp()
    status = 0
    n = size(x,1)
    p = size(x,2)
    if (n <= 1 .or. p < 1 .or. n < p) then
      status = 1
      return
    end if
    if (size(x_transformed,1) /= n .or. size(x_transformed,2) /= p) then
      status = 1
      return
    end if
    if (size(r_transform,1) /= p .or. size(r_transform,2) /= p) then
      status = 1
      return
    end if
    if (size(ri_transform,1) /= p .or. size(ri_transform,2) /= p .or. size(xbar) /= p) then
      status = 1
      return
    end if
    if (any(ieee_is_nan(x))) then
      status = 1
      return
    end if
    allocate(a(n,p), qmat(n,p), rr(p,p), v(n), invrr(p,p))
    do j = 1, p
      xbar(j) = sum(x(:,j)) / real(n, dp)
      a(:,j) = x(:,j) - xbar(j)
    end do
    qmat = 0.0_dp
    rr = 0.0_dp
    tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(a)))
    do j = 1, p
      v = a(:,j)
      do i = 1, j - 1
        rij = dot_product(qmat(:,i), v)
        rr(i,j) = rij
        v = v - rij * qmat(:,i)
      end do
      normv = sqrt(dot_product(v, v))
      if (normv <= tol) then
        status = 2
        return
      end if
      rr(j,j) = normv
      qmat(:,j) = v / normv
    end do
    sn = sqrt(real(n - 1, dp))
    x_transformed = qmat * sn
    ri_transform = rr / sn
    call invert_upper_triangular(ri_transform, invrr, status)
    if (status /= 0) then
      x_transformed = nan_dp()
      r_transform = nan_dp()
      ri_transform = nan_dp()
      xbar = nan_dp()
      return
    end if
    r_transform = invrr
  end subroutine qrxcenter

  pure subroutine invert_upper_triangular(a, ainv, status)
    real(dp), intent(in) :: a(:,:) !! Nonsingular upper-triangular square matrix.
    real(dp), intent(out) :: ainv(:,:) !! Inverse of a when status is zero.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shape or a zero diagonal.
    integer :: i, j, k, n
    real(dp) :: s

    ainv = 0.0_dp
    status = 0
    n = size(a,1)
    if (size(a,2) /= n .or. any(shape(ainv) /= shape(a))) then
      status = 1
      return
    end if
    do i = n, 1, -1
      if (a(i,i) == 0.0_dp) then
        status = 2
        return
      end if
      ainv(i,i) = 1.0_dp / a(i,i)
      do j = i + 1, n
        s = 0.0_dp
        do k = i + 1, j
          s = s + a(i,k) * ainv(k,j)
        end do
        ainv(i,j) = -s / a(i,i)
      end do
    end do
  end subroutine invert_upper_triangular


  pure subroutine solvet_vector(a, b, x, tol, status)
    real(dp), intent(in) :: a(:,:) !! Design/coefficient matrix with observations in rows and unknowns in columns.
    real(dp), intent(in) :: b(:) !! Right-hand-side vector with one value for each row of a.
    real(dp), intent(out) :: x(:) !! Least-squares coefficient vector; length must equal the number of columns of a.
    real(dp), intent(in), optional :: tol !! Relative rank tolerance; defaults to 1e-9, matching the R function default.
    integer, intent(out) :: status !! Zero on success; 1 for shape errors and 2 for an apparently rank-deficient matrix.
    real(dp), allocatable :: q(:,:), r(:,:), v(:)
    real(dp) :: eps, scale, rii
    integer :: i, j, m, n

    status = 0
    x = nan_dp()
    m = size(a, 1)
    n = size(a, 2)
    if (size(b) /= m .or. size(x) /= n .or. m < n .or. n < 1) then
      status = 1
      return
    end if
    eps = 1.0e-9_dp
    if (present(tol)) eps = tol
    allocate(q(m,n), r(n,n), v(m))
    q = 0.0_dp
    r = 0.0_dp
    scale = max(1.0_dp, maxval(abs(a)))
    do j = 1, n
      v = a(:,j)
      do i = 1, j - 1
        r(i,j) = dot_product(q(:,i), v)
        v = v - r(i,j) * q(:,i)
      end do
      rii = sqrt(max(0.0_dp, dot_product(v, v)))
      r(j,j) = rii
      if (rii <= eps * scale) then
        status = 2
        return
      end if
      q(:,j) = v / rii
    end do
    x = matmul(transpose(q), b)
    do i = n, 1, -1
      if (i < n) x(i) = x(i) - dot_product(r(i,i+1:n), x(i+1:n))
      x(i) = x(i) / r(i,i)
    end do
  end subroutine solvet_vector

  pure subroutine solvet_matrix(a, b, x, tol, status)
    real(dp), intent(in) :: a(:,:) !! Design/coefficient matrix with observations in rows and unknowns in columns.
    real(dp), intent(in) :: b(:,:) !! Right-hand-side matrix; its rows must match a and each column is solved independently.
    real(dp), intent(out) :: x(:,:) !! Least-squares coefficient matrix with one row per column of a.
    real(dp), intent(in), optional :: tol !! Relative rank tolerance; defaults to 1e-9, matching the R function default.
    integer, intent(out) :: status !! Zero on success; 1 for shape errors and 2 for an apparently rank-deficient matrix.
    real(dp), allocatable :: q(:,:), r(:,:), y(:,:)
    real(dp) :: eps, scale, rii
    integer :: i, j, k, m, n, nrhs

    status = 0
    x = nan_dp()
    m = size(a, 1)
    n = size(a, 2)
    nrhs = size(b, 2)
    if (size(b,1) /= m .or. size(x,1) /= n .or. size(x,2) /= nrhs .or. m < n .or. n < 1) then
      status = 1
      return
    end if
    eps = 1.0e-9_dp
    if (present(tol)) eps = tol
    allocate(q(m,n), r(n,n), y(n,nrhs))
    q = 0.0_dp
    r = 0.0_dp
    scale = max(1.0_dp, maxval(abs(a)))
    do j = 1, n
      q(:,j) = a(:,j)
      do i = 1, j - 1
        r(i,j) = dot_product(q(:,i), q(:,j))
        q(:,j) = q(:,j) - r(i,j) * q(:,i)
      end do
      rii = sqrt(max(0.0_dp, dot_product(q(:,j), q(:,j))))
      r(j,j) = rii
      if (rii <= eps * scale) then
        status = 2
        return
      end if
      q(:,j) = q(:,j) / rii
    end do
    y = matmul(transpose(q), b)
    x = y
    do k = 1, nrhs
      do i = n, 1, -1
        if (i < n) x(i,k) = x(i,k) - dot_product(r(i,i+1:n), x(i+1:n,k))
        x(i,k) = x(i,k) / r(i,i)
      end do
    end do
  end subroutine solvet_matrix

  pure subroutine solvet_inverse(a, ainv, tol, status)
    real(dp), intent(in) :: a(:,:) !! Square nonsingular matrix to invert through the solvet QR path.
    real(dp), intent(out) :: ainv(:,:) !! Matrix inverse when status is zero.
    real(dp), intent(in), optional :: tol !! Relative rank tolerance; defaults to 1e-9.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shape or rank deficiency.
    real(dp), allocatable :: ident(:,:)
    integer :: i, n

    ainv = nan_dp()
    n = size(a,1)
    if (size(a,2) /= n .or. any(shape(ainv) /= [n,n]) .or. n < 1) then
      status = 1
      return
    end if
    allocate(ident(n,n))
    ident = 0.0_dp
    do i = 1, n
      ident(i,i) = 1.0_dp
    end do
    call solvet_matrix(a, ident, ainv, tol, status)
  end subroutine solvet_inverse

  pure subroutine smean_cl_normal(x, mult, stats)
    real(dp), intent(in) :: x(:) !! Sample values; NaNs are removed as in the R default na.rm=TRUE behavior.
    real(dp), intent(in) :: mult !! Confidence multiplier, e.g. a Student-t quantile computed by the caller.
    real(dp), intent(out) :: stats(3) !! Mean, lower confidence limit, and upper confidence limit.
    real(dp) :: xbar, se
    integer :: i, n

    n = count(.not. ieee_is_nan(x))
    if (n == 0) then
      stats = nan_dp()
      return
    end if
    xbar = 0.0_dp
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) xbar = xbar + x(i)
    end do
    xbar = xbar / real(n, dp)
    if (n < 2) then
      stats = [xbar, nan_dp(), nan_dp()]
      return
    end if
    se = 0.0_dp
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) se = se + (x(i) - xbar)**2
    end do
    se = sqrt(se / real(n * (n - 1), dp))
    stats = [xbar, xbar - mult * se, xbar + mult * se]
  end subroutine smean_cl_normal

  pure subroutine which_close_k(x, w, k, jitter_u, pick_u, indices, status)
    real(dp), intent(in) :: x(:) !! Reference values; values must be finite and no missing values are allowed.
    real(dp), intent(in) :: w(:) !! Query values for which one of the k closest references is selected.
    integer, intent(in) :: k !! Number of nearest candidates; must be between one and size(x).
    real(dp), intent(in) :: jitter_u(:) !! Uniform draws in [0,1] used to break ties in x, one per reference value.
    real(dp), intent(in) :: pick_u(:) !! Uniform draws in [0,1] selecting among k nearest candidates, one per query.
    integer, intent(out) :: indices(:) !! One-based indices into x selected for each query value.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shapes, k, uniforms, or non-finite inputs.
    real(dp), allocatable :: xs(:), diffs(:), work(:)
    integer, allocatable :: order_idx(:)
    real(dp) :: mindif, temp
    integer :: i, j, n, m, p, itmp

    status = 0
    indices = 0
    n = size(x)
    m = size(w)
    if (n < 1 .or. k < 1 .or. k > n .or. size(jitter_u) /= n .or. size(pick_u) /= m .or. size(indices) /= m) then
      status = 1
      return
    end if
    if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(w)) .or. any(jitter_u < 0.0_dp) .or. &
        any(jitter_u > 1.0_dp) .or. any(pick_u < 0.0_dp) .or. any(pick_u > 1.0_dp)) then
      status = 2
      return
    end if
    allocate(xs(n), diffs(n), work(n), order_idx(n))
    work = x
    call sort_real(work)
    mindif = huge(1.0_dp)
    do i = 1, n - 1
      if (work(i+1) > work(i)) mindif = min(mindif, work(i+1) - work(i))
    end do
    if (mindif == huge(1.0_dp)) mindif = 1.0_dp
    xs = x + (2.0_dp * jitter_u - 1.0_dp) * mindif / 100.0_dp
    do i = 1, m
      diffs = abs(w(i) - xs)
      order_idx = [(j, j=1,n)]
      do j = 1, n - 1
        do p = j + 1, n
          if (diffs(p) < diffs(j)) then
            temp = diffs(j)
            diffs(j) = diffs(p)
            diffs(p) = temp
            itmp = order_idx(j)
            order_idx(j) = order_idx(p)
            order_idx(p) = itmp
          end if
        end do
      end do
      if (k == 1) then
        indices(i) = order_idx(1)
      else
        p = min(k, int(pick_u(i) * real(k,dp)) + 1)
        indices(i) = order_idx(p)
      end if
    end do
  end subroutine which_close_k


  pure subroutine largest_eigen_symmetric(a, eigenvalue, eigenvector, status)
    real(dp), intent(in) :: a(:,:) !! Real symmetric matrix whose largest algebraic eigenpair is required.
    real(dp), intent(out) :: eigenvalue !! Largest algebraic eigenvalue on successful convergence.
    real(dp), allocatable, intent(out) :: eigenvector(:) !! Unit eigenvector corresponding to eigenvalue.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shape or failure to converge.
    real(dp), allocatable :: d(:,:), v(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, sn, maxoff, dpi, dqi, vip, viq
    integer :: i, j, p, q, n, sweep, imax

    n = size(a,1)
    status = 0
    eigenvalue = nan_dp()
    allocate(eigenvector(max(0,n)))
    if (n < 1 .or. size(a,2) /= n) then
      status = 1
      return
    end if
    allocate(d(n,n), v(n,n))
    d = 0.5_dp * (a + transpose(a))
    v = 0.0_dp
    do i = 1, n
      v(i,i) = 1.0_dp
    end do
    do sweep = 1, max(50, 20*n*n)
      maxoff = 0.0_dp
      p = 1
      q = min(2,n)
      do j = 2, n
        do i = 1, j - 1
          if (abs(d(i,j)) > maxoff) then
            maxoff = abs(d(i,j))
            p = i
            q = j
          end if
        end do
      end do
      if (maxoff <= 1.0e-13_dp * max(1.0_dp, maxval(abs(d)))) exit
      app = d(p,p)
      aqq = d(q,q)
      apq = d(p,q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t)
      sn = t * c
      do i = 1, n
        if (i /= p .and. i /= q) then
          dpi = d(i,p)
          dqi = d(i,q)
          d(i,p) = c*dpi - sn*dqi
          d(p,i) = d(i,p)
          d(i,q) = sn*dpi + c*dqi
          d(q,i) = d(i,q)
        end if
      end do
      d(p,p) = c*c*app - 2.0_dp*c*sn*apq + sn*sn*aqq
      d(q,q) = sn*sn*app + 2.0_dp*c*sn*apq + c*c*aqq
      d(p,q) = 0.0_dp
      d(q,p) = 0.0_dp
      do i = 1, n
        vip = v(i,p)
        viq = v(i,q)
        v(i,p) = c*vip - sn*viq
        v(i,q) = sn*vip + c*viq
      end do
    end do
    if (sweep > max(50,20*n*n)) then
      status = 2
      eigenvector = nan_dp()
      return
    end if
    imax = 1
    do i = 2, n
      if (d(i,i) > d(imax,imax)) imax = i
    end do
    eigenvalue = d(imax,imax)
    eigenvector = v(:,imax)
  end subroutine largest_eigen_symmetric

  pure subroutine pc1(x, scores, coefficients, fraction_variance, status, hi)
    real(dp), intent(in) :: x(:,:) !! Input observations by variables; rows containing any NaN are omitted.
    real(dp), allocatable, intent(out) :: scores(:) !! PC1 scores for complete rows, optionally rescaled to [0,hi].
    real(dp), allocatable, intent(out) :: coefficients(:) !! Intercept and slopes reproducing scores from complete-case variables.
    real(dp), intent(out) :: fraction_variance !! Fraction of standardized total variance explained by the first component.
    integer, intent(out) :: status !! Zero on success; nonzero for too little data, constant columns, or rank deficiency.
    real(dp), intent(in), optional :: hi !! Positive upper bound for oriented score rescaling to [0,hi].
    real(dp), allocatable :: xo(:,:), z(:,:), means(:), sds(:), v(:), corr(:,:), design(:,:)
    real(dp) :: eig, lo, hi_score
    integer :: i, j, m, n, nc, st, neg

    status = 0
    fraction_variance = nan_dp()
    n = size(x,1)
    m = size(x,2)
    nc = 0
    do i = 1, n
      if (.not. any(ieee_is_nan(x(i,:)))) nc = nc + 1
    end do
    if (m < 1 .or. nc < 2) then
      status = 1
      allocate(scores(0), coefficients(0))
      return
    end if
    allocate(xo(nc,m), z(nc,m), means(m), sds(m), v(m), design(nc,m+1), coefficients(m+1), scores(nc))
    nc = 0
    do i = 1, n
      if (.not. any(ieee_is_nan(x(i,:)))) then
        nc = nc + 1
        xo(nc,:) = x(i,:)
      end if
    end do
    do j = 1, m
      means(j) = sum(xo(:,j)) / real(nc,dp)
      sds(j) = sqrt(sum((xo(:,j) - means(j))**2) / real(nc - 1,dp))
      if (sds(j) <= tiny(1.0_dp)) then
        status = 2
        scores = nan_dp()
        coefficients = nan_dp()
        return
      end if
      z(:,j) = (xo(:,j) - means(j)) / sds(j)
    end do
    allocate(corr(m,m))
    corr = matmul(transpose(z), z) / real(nc - 1,dp)
    call largest_eigen_symmetric(corr, eig, v, st)
    if (st /= 0) then
      status = 3
      scores = nan_dp()
      coefficients = nan_dp()
      return
    end if
    fraction_variance = eig / real(m,dp)
    scores = matmul(z,v)
    design(:,1) = 1.0_dp
    design(:,2:m+1) = xo
    call solvet_vector(design, scores, coefficients, 1.0e-9_dp, st)
    if (st /= 0) then
      status = 4
      scores = nan_dp()
      coefficients = nan_dp()
      return
    end if
    if (present(hi)) then
      if (hi > 0.0_dp) then
        neg = count(coefficients(2:) < 0.0_dp)
        if (real(neg,dp) >= real(m,dp) / 2.0_dp) scores = -scores
        lo = minval(scores)
        hi_score = maxval(scores)
        if (hi_score > lo) scores = hi * (scores - lo) / (hi_score - lo)
        call solvet_vector(design, scores, coefficients, 1.0e-9_dp, st)
        if (st /= 0) status = 4
      end if
    end if
  end subroutine pc1


  pure subroutine lm_fit_qr_bare(x, y, intercept, coefficients, residuals, &
                                 fitted_values, rsquared, status)
    real(dp), intent(in) :: x(:,:) !! Predictor matrix, observations by predictors; must contain no NaNs.
    real(dp), intent(in) :: y(:) !! Response vector with one value per row of x; must contain no NaNs.
    logical, intent(in) :: intercept !! If true, prepend an intercept column as in Hmisc lm.fit.qr.bare.
    real(dp), allocatable, intent(out) :: coefficients(:) !! Least-squares coefficients, including intercept when requested.
    real(dp), allocatable, intent(out) :: residuals(:) !! Response minus fitted values for each observation.
    real(dp), allocatable, intent(out) :: fitted_values(:) !! Least-squares fitted response for each observation.
    real(dp), intent(out) :: rsquared !! 1 - SSE/SST, using response deviations from their sample mean.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, NaNs, or singular design.
    real(dp), allocatable :: a(:,:)
    real(dp) :: sse, sst, ybar
    integer :: n, p, st

    n = size(x,1)
    p = size(x,2)
    status = 0
    rsquared = nan_dp()
    if (size(y) /= n .or. n < 1 .or. any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) then
      status = 1
      allocate(coefficients(0), residuals(0), fitted_values(0))
      return
    end if
    if (intercept) then
      allocate(a(n,p+1), coefficients(p+1))
      a(:,1) = 1.0_dp
      if (p > 0) a(:,2:p+1) = x
    else
      allocate(a(n,p), coefficients(p))
      if (p > 0) a = x
    end if
    allocate(residuals(n), fitted_values(n))
    call solvet_vector(a, y, coefficients, 1.0e-7_dp, st)
    if (st /= 0) then
      status = 2
      coefficients = nan_dp()
      residuals = nan_dp()
      fitted_values = nan_dp()
      return
    end if
    fitted_values = matmul(a, coefficients)
    residuals = y - fitted_values
    sse = sum(residuals**2)
    ybar = sum(y) / real(n,dp)
    sst = sum((y-ybar)**2)
    if (sst <= tiny(1.0_dp)) then
      rsquared = nan_dp()
    else
      rsquared = 1.0_dp - sse / sst
    end if
  end subroutine lm_fit_qr_bare

  pure subroutine abs_error_pred(y, lp, differences, ratios, status)
    real(dp), intent(in) :: y(:) !! Observed outcomes; pairs containing NaN are omitted.
    real(dp), intent(in) :: lp(:) !! Predicted values or linear predictors paired one-to-one with y.
    real(dp), intent(out) :: differences(3,2) !! Absolute-difference summaries: rows are spread/error/outcome, columns mean/median.
    real(dp), intent(out) :: ratios(2,2) !! Spread/error divided by outcome spread; columns are mean and median ratios.
    integer, intent(out) :: status !! Zero on success; nonzero for mismatched lengths or no complete pairs.
    real(dp), allocatable :: yy(:), pp(:), ay(:), ap(:), ae(:)
    real(dp) :: my, mp, meant, meanr, meane, medt, medr, mede
    integer :: i, n, nc

    differences = nan_dp()
    ratios = nan_dp()
    status = 0
    n = size(y)
    if (size(lp) /= n) then
      status = 1
      return
    end if
    nc = 0
    do i = 1, n
      if (.not. ieee_is_nan(y(i)) .and. .not. ieee_is_nan(lp(i))) nc = nc + 1
    end do
    if (nc == 0) then
      status = 2
      return
    end if
    allocate(yy(nc), pp(nc), ay(nc), ap(nc), ae(nc))
    nc = 0
    do i = 1, n
      if (.not. ieee_is_nan(y(i)) .and. .not. ieee_is_nan(lp(i))) then
        nc = nc + 1
        yy(nc) = y(i)
        pp(nc) = lp(i)
      end if
    end do
    my = median_sorted_copy(yy)
    mp = median_sorted_copy(pp)
    ay = abs(yy-my)
    ap = abs(pp-mp)
    ae = abs(pp-yy)
    meant = sum(ay) / real(size(ay),dp)
    meanr = sum(ap) / real(size(ap),dp)
    meane = sum(ae) / real(size(ae),dp)
    medt = median_sorted_copy(ay)
    medr = median_sorted_copy(ap)
    mede = median_sorted_copy(ae)
    differences(:,1) = [meanr, meane, meant]
    differences(:,2) = [medr, mede, medt]
    if (meant > 0.0_dp) ratios(:,1) = [meanr/meant, meane/meant]
    if (medt > 0.0_dp) ratios(:,2) = [medr/medt, mede/medt]
  end subroutine abs_error_pred

  pure subroutine rcorrcens_summary(x, y, event, outx, c_index, dxy, abs_dxy, &
                                    sd, z, p_value, status)
    real(dp), intent(in) :: x(:) !! Numeric predictor values paired with censored outcomes.
    real(dp), intent(in) :: y(:) !! Observed event or censoring times.
    logical, intent(in) :: event(:) !! True for observed events and false for censored observations.
    logical, intent(in) :: outx !! If true, predictor ties are excluded from relevant pairs.
    real(dp), intent(out) :: c_index !! Harrell concordance index returned by rcorr.cens.
    real(dp), intent(out) :: dxy !! Somers Dxy = 2*(C - 0.5).
    real(dp), intent(out) :: abs_dxy !! Absolute value of Dxy.
    real(dp), intent(out) :: sd !! Quade standard-error estimate returned by rcorr.cens.
    real(dp), intent(out) :: z !! Absolute Dxy divided by its standard error.
    real(dp), intent(out) :: p_value !! Two-sided normal-approximation P value.
    integer, intent(out) :: status !! Zero on success; nonzero when concordance or standard error is undefined.
    real(dp) :: gamma
    real(dp) :: nrel, nconc, nuncert

    call rcorr_cens(x, y, event, outx, nrel, nconc, nuncert, c_index, gamma, sd)
    dxy = gamma
    abs_dxy = abs(dxy)
    if (ieee_is_nan(c_index) .or. ieee_is_nan(sd) .or. sd <= 0.0_dp) then
      z = nan_dp()
      p_value = nan_dp()
      status = 1
      return
    end if
    z = abs_dxy / sd
    p_value = 2.0_dp * norm_cdf(-z)
    status = 0
  end subroutine rcorrcens_summary


  pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a !! First beta-shape parameter; must be positive.
    real(dp), intent(in) :: b !! Second beta-shape parameter; must be positive.
    real(dp), intent(in) :: x !! Evaluation point in the closed unit interval.
    integer, parameter :: maxit = 200
    real(dp), parameter :: eps = 3.0e-14_dp
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
    real(dp) :: qab, qap, qam, c, d, h, aa, del
    integer :: m, m2

    qab = a + b
    qap = a + 1.0_dp
    qam = a - 1.0_dp
    c = 1.0_dp
    d = 1.0_dp - qab * x / qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp / d
    h = d
    do m = 1, maxit
      m2 = 2 * m
      aa = real(m, dp) * (b - real(m, dp)) * x / &
           ((qam + real(m2, dp)) * (a + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      h = h * d * c
      aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
           ((a + real(m2, dp)) * (qap + real(m2, dp)))
      d = 1.0_dp + aa * d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp + aa / c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
    cf = h
  end function beta_continued_fraction

  pure real(dp) function regularized_beta(x, a, b) result(p)
    real(dp), intent(in) :: x !! Evaluation point in the closed unit interval.
    real(dp), intent(in) :: a !! First positive beta-shape parameter.
    real(dp), intent(in) :: b !! Second positive beta-shape parameter.
    real(dp) :: bt

    if (x <= 0.0_dp) then
      p = 0.0_dp
      return
    end if
    if (x >= 1.0_dp) then
      p = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
             a * log(x) + b * log(1.0_dp - x))
    if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
      p = bt * beta_continued_fraction(a, b, x) / a
    else
      p = 1.0_dp - bt * beta_continued_fraction(b, a, 1.0_dp - x) / b
    end if
  end function regularized_beta

  pure real(dp) function f_upper_tail(fstat, df1, df2) result(p)
    real(dp), intent(in) :: fstat !! Nonnegative F statistic.
    real(dp), intent(in) :: df1 !! Positive numerator degrees of freedom.
    real(dp), intent(in) :: df2 !! Positive denominator degrees of freedom.
    real(dp) :: z

    if (fstat < 0.0_dp .or. df1 <= 0.0_dp .or. df2 <= 0.0_dp) then
      p = nan_dp()
      return
    end if
    z = df2 / (df2 + df1 * fstat)
    p = regularized_beta(z, 0.5_dp * df2, 0.5_dp * df1)
  end function f_upper_tail

  pure subroutine logrank(time, event, group, chisq, hazard_ratio, status)
    real(dp), intent(in) :: time(:) !! Event or censoring times; NaNs are omitted with their paired observations.
    logical, intent(in) :: event(:) !! True for observed events and false for censored observations.
    integer, intent(in) :: group(:) !! Two-group labels; exactly two distinct values are required after NaN removal.
    real(dp), intent(out) :: chisq !! One-degree-of-freedom log-rank chi-square statistic.
    real(dp), intent(out) :: hazard_ratio !! Observed/expected hazard-ratio estimate for second versus first group.
    integer, intent(out) :: status !! Zero on success; nonzero for shape, grouping, or degenerate-risk errors.
    real(dp), allocatable :: t(:)
    integer, allocatable :: g(:), ord(:)
    logical, allocatable :: e(:)
    integer :: n, nc, i, j, k, g1, g2, nr1, nr2, d1, d2, rd, rs, nt
    real(dp) :: oe, vv, e1, e2, o1, o2, den, tt

    chisq = nan_dp()
    hazard_ratio = nan_dp()
    status = 0
    n = size(time)
    if (size(event) /= n .or. size(group) /= n) then
      status = 1
      return
    end if
    nc = count(.not. ieee_is_nan(time))
    if (nc < 2) then
      status = 2
      return
    end if
    allocate(t(nc), e(nc), g(nc), ord(nc))
    k = 0
    do i = 1, n
      if (.not. ieee_is_nan(time(i))) then
        k = k + 1
        t(k) = time(i)
        e(k) = event(i)
        g(k) = group(i)
        ord(k) = k
      end if
    end do
    g1 = g(1)
    g2 = g1
    do i = 2, nc
      if (g(i) /= g1) then
        g2 = g(i)
        exit
      end if
    end do
    if (g2 == g1) then
      status = 3
      return
    end if
    do i = 1, nc
      if (g(i) /= g1 .and. g(i) /= g2) then
        status = 3
        return
      end if
    end do
    do i = 2, nc
      j = i
      do while (j > 1)
        if (t(ord(j - 1)) >= t(ord(j))) exit
        k = ord(j - 1)
        ord(j - 1) = ord(j)
        ord(j) = k
        j = j - 1
      end do
    end do
    nr1 = 0
    nr2 = 0
    oe = 0.0_dp
    vv = 0.0_dp
    e1 = 0.0_dp
    e2 = 0.0_dp
    o1 = 0.0_dp
    o2 = 0.0_dp
    i = 1
    do while (i <= nc)
      tt = t(ord(i))
      j = i
      do while (j <= nc)
        if (t(ord(j)) /= tt) exit
        if (g(ord(j)) == g1) then
          nr1 = nr1 + 1
        else
          nr2 = nr2 + 1
        end if
        j = j + 1
      end do
      d1 = 0
      d2 = 0
      do k = i, j - 1
        if (e(ord(k))) then
          if (g(ord(k)) == g1) then
            d1 = d1 + 1
          else
            d2 = d2 + 1
          end if
        end if
      end do
      rd = d1 + d2
      nt = nr1 + nr2
      if (rd > 0 .and. nt > 0) then
        oe = oe + real(d1, dp) - real(rd, dp) * real(nr1, dp) / real(nt, dp)
        e1 = e1 + real(nr1, dp) * real(rd, dp) / real(nt, dp)
        e2 = e2 + real(nr2, dp) * real(rd, dp) / real(nt, dp)
        o1 = o1 + real(d1, dp)
        o2 = o2 + real(d2, dp)
        rs = nt - rd
        if (nt > 1) then
          vv = vv + real(rd, dp) * real(rs, dp) * real(nr1, dp) * real(nr2, dp) / &
               (real(nt, dp) * real(nt, dp) * real(nt - 1, dp))
        end if
      end if
      i = j
    end do
    if (vv <= 0.0_dp .or. e1 <= 0.0_dp .or. e2 <= 0.0_dp .or. o1 <= 0.0_dp) then
      status = 4
      return
    end if
    chisq = oe * oe / vv
    den = o1 / e1
    if (den <= 0.0_dp) then
      status = 4
      chisq = nan_dp()
      return
    end if
    hazard_ratio = (o2 / e2) / den
  end subroutine logrank

  pure subroutine spearman2_numeric(x, y, power, rho2, fstat, df1, df2, &
                                    p_value, adjusted_rho2, n_complete, status)
    real(dp), intent(in) :: x(:) !! Numeric predictor values; NaN pairs are omitted.
    real(dp), intent(in) :: y(:) !! Numeric response values; NaN pairs are omitted before ranking.
    integer, intent(in) :: power !! Predictor rank-polynomial degree, either 1 or 2 as in Hmisc spearman2.default.
    real(dp), intent(out) :: rho2 !! Squared rank-association measure.
    real(dp), intent(out) :: fstat !! F statistic for testing the rank association.
    real(dp), intent(out) :: df1 !! Numerator degrees of freedom, equal to power after numeric simplification.
    real(dp), intent(out) :: df2 !! Denominator degrees of freedom, n - power - 1.
    real(dp), intent(out) :: p_value !! Upper-tail F-test probability.
    real(dp), intent(out) :: adjusted_rho2 !! Degrees-of-freedom adjusted squared association.
    integer, intent(out) :: n_complete !! Number of complete x/y pairs used.
    integer, intent(out) :: status !! Zero on success; nonzero for bad shapes, power, or insufficient/rank-deficient data.
    real(dp), allocatable :: xx(:), yy(:), rx(:), ry(:), design(:,:), coef(:), res(:), fit(:)
    real(dp) :: rxy, rsq
    integer :: i, k, lmstat, u

    rho2 = nan_dp()
    fstat = nan_dp()
    df1 = 0.0_dp
    df2 = 0.0_dp
    p_value = nan_dp()
    adjusted_rho2 = nan_dp()
    n_complete = 0
    status = 0
    if (size(x) /= size(y)) then
      status = 1
      return
    end if
    if (power < 1 .or. power > 2) then
      status = 2
      return
    end if
    n_complete = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (n_complete < 3) then
      status = 3
      df2 = real(n_complete, dp)
      return
    end if
    allocate(xx(n_complete), yy(n_complete), rx(n_complete), ry(n_complete))
    k = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        k = k + 1
        xx(k) = x(i)
        yy(k) = y(i)
      end if
    end do
    call midrank(xx, rx)
    call midrank(yy, ry)
    u = 1
    do i = 2, n_complete
      if (all(xx(i) /= xx(1:i-1))) u = u + 1
    end do
    if (u < 3) then
      df1 = 1.0_dp
      rxy = sum((rx - sum(rx)/real(n_complete,dp)) * &
                (ry - sum(ry)/real(n_complete,dp))) / &
            sqrt(sum((rx - sum(rx)/real(n_complete,dp))**2) * &
                 sum((ry - sum(ry)/real(n_complete,dp))**2))
      rsq = rxy * rxy
    else if (power == 1) then
      df1 = 1.0_dp
      rxy = sum((rx - sum(rx)/real(n_complete,dp)) * &
                (ry - sum(ry)/real(n_complete,dp))) / &
            sqrt(sum((rx - sum(rx)/real(n_complete,dp))**2) * &
                 sum((ry - sum(ry)/real(n_complete,dp))**2))
      rsq = rxy * rxy
    else
      df1 = 2.0_dp
      allocate(design(n_complete,2), coef(3), res(n_complete), fit(n_complete))
      design(:,1) = rx
      design(:,2) = rx * rx
      call lm_fit_qr_bare(design, ry, .true., coef, res, fit, rsq, lmstat)
      if (lmstat /= 0) then
        status = 4
        return
      end if
    end if
    df2 = real(n_complete, dp) - df1 - 1.0_dp
    if (df2 <= 0.0_dp .or. rsq >= 1.0_dp) then
      if (rsq >= 1.0_dp .and. df2 > 0.0_dp) then
        rho2 = rsq
        fstat = huge(1.0_dp)
        p_value = 0.0_dp
        adjusted_rho2 = 1.0_dp
        return
      end if
      status = 5
      return
    end if
    rho2 = rsq
    fstat = (rho2 / df1) / ((1.0_dp - rho2) / df2)
    p_value = f_upper_tail(fstat, df1, df2)
    adjusted_rho2 = 1.0_dp - (1.0_dp - rho2) * &
                    real(n_complete - 1, dp) / df2
  end subroutine spearman2_numeric

  pure subroutine equal_bins(widths, subwidths, counts, result, status)
    integer, intent(in) :: widths(:) !! Total display widths for each outer group, including column spacers.
    integer, intent(in) :: subwidths(:,:) !! Existing widths for subcolumns; only the first counts(j) entries of row j are used.
    integer, intent(in) :: counts(:) !! Number of valid subcolumns in each row of subwidths; each value must be positive.
    integer, allocatable, intent(out) :: result(:) !! Concatenated adjusted subcolumn widths, matching Hmisc equalBins ordering.
    integer, intent(out) :: status !! Zero on success; nonzero for incompatible dimensions or invalid counts/widths.
    integer :: j, i, k, total, adjusted, swsum, divv, modv

    status = 0
    if (size(widths) /= size(counts) .or. size(subwidths,1) /= size(widths)) then
      status = 1
      allocate(result(0))
      return
    end if
    total = sum(counts)
    if (total < 0 .or. any(counts <= 0) .or. any(counts > size(subwidths,2))) then
      status = 2
      allocate(result(0))
      return
    end if
    allocate(result(total))
    k = 0
    do j = 1, size(widths)
      adjusted = widths(j) - counts(j) + 1
      if (adjusted < 0) then
        status = 3
        result = 0
        return
      end if
      swsum = sum(subwidths(j,1:counts(j)))
      if (swsum < adjusted) then
        divv = adjusted / counts(j)
        modv = modulo(adjusted, counts(j))
        do i = 1, counts(j)
          k = k + 1
          result(k) = divv
          if (i <= modv) result(k) = result(k) + 1
        end do
      else
        do i = 1, counts(j)
          k = k + 1
          result(k) = subwidths(j,i)
        end do
      end if
    end do
  end subroutine equal_bins

  pure subroutine jitter2_numeric(x, fill, limit, eps, uniforms, out, status)
    real(dp), intent(in) :: x(:) !! Numeric values to spread when exact or eps-rounded ties occur; NaNs are preserved.
    real(dp), intent(in) :: fill !! Fraction of neighboring spacing occupied by each tied run, matching jitter2.default.
    real(dp), intent(in) :: limit !! Positive maximum half-width; nonpositive values disable the explicit cap.
    real(dp), intent(in) :: eps !! If positive, ties are defined after rounding x/eps to nearest integer multiples of eps.
    real(dp), intent(in) :: uniforms(:) !! Caller-supplied tie-order uniforms; at least one is required per tied observation.
    real(dp), intent(out) :: out(:) !! Jittered values in original observation order with NaNs unchanged.
    integer, intent(out) :: status !! Zero on success; nonzero for output shape or insufficient supplied uniforms.
    real(dp), allocatable :: xs(:), keys(:), values(:), dist(:), half(:), vals(:), us(:)
    integer, allocatable :: orig(:), lengths(:), runstart(:), perm(:)
    integer :: n, nc, i, j, k, nr, ntied, r, pos, m, a, b, tmpi
    real(dp) :: tmp, dleft, dright, h, step

    status = 0
    n = size(x)
    if (size(out) /= n) then
      status = 1
      return
    end if
    out = x
    nc = count(.not. ieee_is_nan(x))
    if (nc < 2) return
    allocate(xs(nc), keys(nc), orig(nc))
    k = 0
    do i = 1, n
      if (.not. ieee_is_nan(x(i))) then
        k = k + 1
        xs(k) = x(i)
        orig(k) = i
      end if
    end do
    do i = 2, nc
      j = i
      do while (j > 1)
        if (xs(j-1) <= xs(j)) exit
        tmp = xs(j-1)
        xs(j-1) = xs(j)
        xs(j) = tmp
        tmpi = orig(j-1)
        orig(j-1) = orig(j)
        orig(j) = tmpi
        j = j - 1
      end do
    end do
    if (eps > 0.0_dp) then
      keys = anint(xs / eps) * eps
    else
      keys = xs
    end if
    nr = 1
    do i = 2, nc
      if (keys(i) /= keys(i-1)) nr = nr + 1
    end do
    allocate(values(nr), lengths(nr), runstart(nr), dist(nr), half(nr))
    r = 1
    runstart(1) = 1
    values(1) = keys(1)
    lengths = 0
    lengths(1) = 1
    do i = 2, nc
      if (keys(i) == keys(i-1)) then
        lengths(r) = lengths(r) + 1
      else
        r = r + 1
        runstart(r) = i
        values(r) = keys(i)
        lengths(r) = 1
      end if
    end do
    if (maxval(lengths) < 2 .or. nr < 2) return
    do r = 1, nr
      if (r == 1) then
        dleft = abs(values(2) - values(1))
      else
        dleft = abs(values(r) - values(r-1))
      end if
      if (r == nr) then
        dright = abs(values(nr) - values(nr-1))
      else
        dright = abs(values(r+1) - values(r))
      end if
      dist(r) = min(dleft, dright)
      half(r) = 0.5_dp * fill * dist(r)
      if (limit > 0.0_dp) half(r) = min(half(r), limit)
    end do
    ntied = 0
    do r = 1, nr
      if (lengths(r) > 1) ntied = ntied + lengths(r)
    end do
    if (size(uniforms) < ntied) then
      status = 2
      return
    end if
    allocate(vals(ntied), us(ntied), perm(ntied))
    pos = 0
    k = 0
    do r = 1, nr
      if (lengths(r) <= 1) cycle
      m = lengths(r)
      h = half(r)
      if (m > 1) then
        step = 2.0_dp * h / real(m - 1, dp)
      else
        step = 0.0_dp
      end if
      a = pos + 1
      b = pos + m
      do i = 1, m
        pos = pos + 1
        vals(pos) = values(r) - h + real(i - 1, dp) * step
        k = k + 1
        us(pos) = uniforms(k)
        perm(pos) = i
      end do
      do i = a + 1, b
        j = i
        do while (j > a)
          if (us(j-1) <= us(j)) exit
          tmp = us(j-1)
          us(j-1) = us(j)
          us(j) = tmp
          tmpi = perm(j-1)
          perm(j-1) = perm(j)
          perm(j) = tmpi
          j = j - 1
        end do
      end do
      do i = 1, m
        out(orig(runstart(r) + i - 1)) = vals(a + perm(a + i - 1) - 1)
      end do
    end do
  end subroutine jitter2_numeric


  pure subroutine hdquantile(x, probs, estimates, status)
    real(dp), intent(in) :: x(:) !! Numeric sample; NaNs are omitted before estimating Harrell-Davis quantiles.
    real(dp), intent(in) :: probs(:) !! Quantile probabilities in the closed unit interval.
    real(dp), intent(out) :: estimates(:) !! Harrell-Davis quantile estimates, with one element per probability.
    integer, intent(out) :: status !! Zero on success; nonzero for shape, probability, or insufficient-data errors.
    real(dp), allocatable :: xs(:)
    real(dp) :: a, b, lo, hi, weight, m, temp
    integer :: i, j, k, n, nc

    estimates = nan_dp()
    status = 0
    if (size(estimates) /= size(probs)) then
      status = 1
      return
    end if
    if (any(probs < 0.0_dp) .or. any(probs > 1.0_dp) .or. any(ieee_is_nan(probs))) then
      status = 2
      return
    end if
    nc = count(.not. ieee_is_nan(x))
    if (nc < 2) then
      status = 3
      return
    end if
    allocate(xs(nc))
    j = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i))) cycle
      j = j + 1
      xs(j) = x(i)
    end do
    do i = 2, nc
      temp = xs(i)
      k = i - 1
      do while (k >= 1)
        if (xs(k) <= temp) exit
        xs(k + 1) = xs(k)
        k = k - 1
      end do
      xs(k + 1) = temp
    end do
    n = size(xs)
    m = real(n + 1, dp)
    do j = 1, size(probs)
      if (probs(j) == 0.0_dp) then
        estimates(j) = xs(1)
      else if (probs(j) == 1.0_dp) then
        estimates(j) = xs(n)
      else
        a = probs(j) * m
        b = (1.0_dp - probs(j)) * m
        estimates(j) = 0.0_dp
        do i = 1, n
          lo = real(i - 1, dp) / real(n, dp)
          hi = real(i, dp) / real(n, dp)
          weight = regularized_beta(hi, a, b) - regularized_beta(lo, a, b)
          estimates(j) = estimates(j) + xs(i) * weight
        end do
      end if
    end do
  end subroutine hdquantile

  pure real(dp) function gamma_q(a, x) result(q)
    real(dp), intent(in) :: a !! Positive gamma shape parameter.
    real(dp), intent(in) :: x !! Nonnegative evaluation point.
    integer, parameter :: maxit = 300
    real(dp), parameter :: eps = 3.0e-14_dp
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
    real(dp) :: ap, del, sumv, b, c, d, h, an, gln
    integer :: i

    if (a <= 0.0_dp .or. x < 0.0_dp .or. ieee_is_nan(a) .or. ieee_is_nan(x)) then
      q = nan_dp()
      return
    end if
    if (x == 0.0_dp) then
      q = 1.0_dp
      return
    end if
    gln = log_gamma(a)
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp / a
      del = sumv
      do i = 1, maxit
        ap = ap + 1.0_dp
        del = del * x / ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv) * eps) exit
      end do
      q = 1.0_dp - sumv * exp(-x + a * log(x) - gln)
      q = max(0.0_dp, min(1.0_dp, q))
      return
    end if
    b = x + 1.0_dp - a
    c = 1.0_dp / fpmin
    d = 1.0_dp / max(abs(b), fpmin)
    if (b < 0.0_dp) d = -d
    h = d
    do i = 1, maxit
      an = -real(i, dp) * (real(i, dp) - a)
      b = b + 2.0_dp
      d = an * d + b
      if (abs(d) < fpmin) d = sign(fpmin, d)
      c = b + an / c
      if (abs(c) < fpmin) c = sign(fpmin, c)
      d = 1.0_dp / d
      del = d * c
      h = h * del
      if (abs(del - 1.0_dp) <= eps) exit
    end do
    q = exp(-x + a * log(x) - gln) * h
    q = max(0.0_dp, min(1.0_dp, q))
  end function gamma_q

  pure subroutine cat_test_chisq(tab, statistic, df, p_value, status)
    real(dp), intent(in) :: tab(:,:) !! Nonnegative contingency-table counts; zero-total rows are discarded.
    real(dp), intent(out) :: statistic !! Pearson chi-square statistic without continuity correction.
    real(dp), intent(out) :: df !! Pearson chi-square degrees of freedom after removing zero-total rows.
    real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid or degenerate tables.
    real(dp), allocatable :: rows(:), cols(:)
    real(dp) :: total, expected
    integer :: i, j, nr

    statistic = nan_dp()
    df = nan_dp()
    p_value = nan_dp()
    status = 0
    if (size(tab, 1) < 2 .or. size(tab, 2) < 2 .or. any(tab < 0.0_dp) .or. any(ieee_is_nan(tab))) then
      status = 1
      return
    end if
    allocate(rows(size(tab,1)), cols(size(tab,2)))
    rows = sum(tab, dim=2)
    nr = count(rows > 0.0_dp)
    if (nr < 2) then
      status = 2
      return
    end if
    cols = sum(tab, dim=1)
    if (count(cols > 0.0_dp) < 2) then
      status = 2
      return
    end if
    total = sum(rows)
    statistic = 0.0_dp
    do i = 1, size(tab,1)
      if (rows(i) <= 0.0_dp) cycle
      do j = 1, size(tab,2)
        if (cols(j) <= 0.0_dp) cycle
        expected = rows(i) * cols(j) / total
        statistic = statistic + (tab(i,j) - expected)**2 / expected
      end do
    end do
    df = real((nr - 1) * (count(cols > 0.0_dp) - 1), dp)
    p_value = gamma_q(0.5_dp * df, 0.5_dp * statistic)
  end subroutine cat_test_chisq

  pure subroutine con_test_kw(group, x, statistic, df1, df2, p_value, status)
    real(dp), intent(in) :: group(:) !! Numeric group scores passed to Hmisc spearman2 with its default power one.
    real(dp), intent(in) :: x(:) !! Numeric response values paired with group scores.
    real(dp), intent(out) :: statistic !! Spearman2 F statistic used by conTestkw.
    real(dp), intent(out) :: df1 !! Numerator degrees of freedom for the F statistic.
    real(dp), intent(out) :: df2 !! Denominator degrees of freedom for the F statistic.
    real(dp), intent(out) :: p_value !! Upper-tail F probability.
    integer, intent(out) :: status !! Zero on success; propagated nonzero status from the numeric spearman2 path.
    real(dp) :: rho2, adjusted
    integer :: n

    call spearman2_numeric(group, x, 1, rho2, statistic, df1, df2, p_value, adjusted, n, status)
  end subroutine con_test_kw

  pure subroutine t_test_cluster(y, cluster, group, conf_level, group_n, cluster_n, group_mean, &
                                 rho, correction, variance_effect, variance_unadjusted, &
                                 design_effect, effect, se, conf_limits, z, p_value, status)
    real(dp), intent(in) :: y(:) !! Numeric outcome values; NaNs are omitted with paired cluster and group labels.
    integer, intent(in) :: cluster(:) !! Cluster identifiers paired with outcomes; identifiers may be arbitrary integers.
    integer, intent(in) :: group(:) !! Treatment-group identifiers; exactly two distinct values are required.
    real(dp), intent(in) :: conf_level !! Two-sided normal confidence level strictly between zero and one.
    integer, intent(out) :: group_n(2) !! Numbers of retained observations in the two sorted treatment groups.
    integer, intent(out) :: cluster_n(2) !! Numbers of distinct clusters in the two treatment groups.
    real(dp), intent(out) :: group_mean(2) !! Outcome means for the two treatment groups.
    real(dp), intent(out) :: rho !! Estimated intracluster correlation before truncation at zero for variance correction.
    real(dp), intent(out) :: correction(2) !! Variance correction factors C1 and C2 for the two treatment groups.
    real(dp), intent(out) :: variance_effect !! Cluster-adjusted variance of the difference in group means.
    real(dp), intent(out) :: variance_unadjusted !! Variance of the difference ignoring cluster correlation.
    real(dp), intent(out) :: design_effect !! Ratio of adjusted to unadjusted effect variance.
    real(dp), intent(out) :: effect !! Mean of group two minus mean of group one.
    real(dp), intent(out) :: se !! Standard error of the group-mean difference.
    real(dp), intent(out) :: conf_limits(2) !! Lower and upper normal confidence limits for the effect.
    real(dp), intent(out) :: z !! Normal Z statistic for the group-mean difference.
    real(dp), intent(out) :: p_value !! Two-sided normal-approximation P value.
    integer, intent(out) :: status !! Zero on success; nonzero for shape, grouping, cluster, or variance degeneracy errors.
    real(dp), allocatable :: yy(:)
    integer, allocatable :: cc(:), gg(:), labels(:), clabels(:), counts(:)
    real(dp), allocatable :: cmeans(:), cvars(:)
    integer :: i, j, k, nc, ng, gidx, ccount, df_msc, df_mse
    real(dp) :: ssc(2), ssw(2), msc, mse, na, ruse, csum, ss, zcrit

    group_n = 0
    cluster_n = 0
    group_mean = nan_dp()
    rho = nan_dp()
    correction = nan_dp()
    variance_effect = nan_dp()
    variance_unadjusted = nan_dp()
    design_effect = nan_dp()
    effect = nan_dp()
    se = nan_dp()
    conf_limits = nan_dp()
    z = nan_dp()
    p_value = nan_dp()
    status = 0
    if (size(cluster) /= size(y) .or. size(group) /= size(y)) then
      status = 1
      return
    end if
    if (conf_level <= 0.0_dp .or. conf_level >= 1.0_dp .or. ieee_is_nan(conf_level)) then
      status = 2
      return
    end if
    nc = count(.not. ieee_is_nan(y))
    if (nc < 2) then
      status = 3
      return
    end if
    allocate(yy(nc), cc(nc), gg(nc))
    j = 0
    do i = 1, size(y)
      if (ieee_is_nan(y(i))) cycle
      j = j + 1
      yy(j) = y(i)
      cc(j) = cluster(i)
      gg(j) = group(i)
    end do
    allocate(labels(nc))
    labels = 0
    ng = 0
    do i = 1, nc
      if (ng == 0 .or. .not. any(labels(1:ng) == gg(i))) then
        ng = ng + 1
        labels(ng) = gg(i)
      end if
    end do
    if (ng /= 2) then
      status = 4
      return
    end if
    if (labels(2) < labels(1)) then
      k = labels(1)
      labels(1) = labels(2)
      labels(2) = k
    end if
    ssc = 0.0_dp
    ssw = 0.0_dp
    do gidx = 1, 2
      group_n(gidx) = count(gg == labels(gidx))
      if (group_n(gidx) < 2) then
        status = 5
        return
      end if
      group_mean(gidx) = sum(yy, mask=gg == labels(gidx)) / real(group_n(gidx), dp)
      allocate(clabels(group_n(gidx)), counts(group_n(gidx)), cmeans(group_n(gidx)), cvars(group_n(gidx)))
      clabels = 0
      counts = 0
      cmeans = 0.0_dp
      cvars = 0.0_dp
      ccount = 0
      do i = 1, nc
        if (gg(i) /= labels(gidx)) cycle
        k = 0
        do j = 1, ccount
          if (clabels(j) == cc(i)) then
            k = j
            exit
          end if
        end do
        if (k == 0) then
          ccount = ccount + 1
          k = ccount
          clabels(k) = cc(i)
        end if
        counts(k) = counts(k) + 1
        cmeans(k) = cmeans(k) + yy(i)
      end do
      cluster_n(gidx) = ccount
      do j = 1, ccount
        if (counts(j) < 2) then
          status = 6
          return
        end if
        cmeans(j) = cmeans(j) / real(counts(j), dp)
        ss = 0.0_dp
        do i = 1, nc
          if (gg(i) == labels(gidx) .and. cc(i) == clabels(j)) ss = ss + (yy(i) - cmeans(j))**2
        end do
        cvars(j) = ss / real(counts(j) - 1, dp)
        ssc(gidx) = ssc(gidx) + real(counts(j), dp) * (cmeans(j) - group_mean(gidx))**2
        ssw(gidx) = ssw(gidx) + real(counts(j) - 1, dp) * cvars(j)
      end do
      deallocate(clabels, counts, cmeans, cvars)
    end do
    df_msc = sum(cluster_n) - 2
    df_mse = sum(group_n) - sum(cluster_n)
    if (df_msc <= 0 .or. df_mse <= 0) then
      status = 7
      return
    end if
    msc = sum(ssc) / real(df_msc, dp)
    mse = sum(ssw) / real(df_mse, dp)
    csum = 0.0_dp
    do gidx = 1, 2
      allocate(clabels(group_n(gidx)), counts(group_n(gidx)))
      clabels = 0
      counts = 0
      ccount = 0
      do i = 1, nc
        if (gg(i) /= labels(gidx)) cycle
        k = 0
        do j = 1, ccount
          if (clabels(j) == cc(i)) then
            k = j
            exit
          end if
        end do
        if (k == 0) then
          ccount = ccount + 1
          k = ccount
          clabels(k) = cc(i)
        end if
        counts(k) = counts(k) + 1
      end do
      csum = csum + sum(real(counts(1:ccount), dp)**2) / real(group_n(gidx), dp)
      deallocate(clabels, counts)
    end do
    na = (real(sum(group_n), dp) - csum) / real(sum(cluster_n) - 1, dp)
    if (msc + (na - 1.0_dp) * mse == 0.0_dp) then
      status = 8
      return
    end if
    rho = (msc - mse) / (msc + (na - 1.0_dp) * mse)
    ruse = max(rho, 0.0_dp)
    do gidx = 1, 2
      allocate(clabels(group_n(gidx)), counts(group_n(gidx)))
      clabels = 0
      counts = 0
      ccount = 0
      do i = 1, nc
        if (gg(i) /= labels(gidx)) cycle
        k = 0
        do j = 1, ccount
          if (clabels(j) == cc(i)) then
            k = j
            exit
          end if
        end do
        if (k == 0) then
          ccount = ccount + 1
          k = ccount
          clabels(k) = cc(i)
        end if
        counts(k) = counts(k) + 1
      end do
      correction(gidx) = sum(real(counts(1:ccount), dp) * &
                              (1.0_dp + real(counts(1:ccount) - 1, dp) * ruse)) / real(group_n(gidx), dp)
      deallocate(clabels, counts)
    end do
    variance_effect = mse * (correction(1) / real(group_n(1), dp) + correction(2) / real(group_n(2), dp))
    variance_unadjusted = mse * (1.0_dp / real(group_n(1), dp) + 1.0_dp / real(group_n(2), dp))
    if (variance_effect <= 0.0_dp .or. variance_unadjusted <= 0.0_dp) then
      status = 9
      return
    end if
    design_effect = variance_effect / variance_unadjusted
    effect = group_mean(2) - group_mean(1)
    se = sqrt(variance_effect)
    zcrit = inv_norm_cdf(0.5_dp * (1.0_dp + conf_level))
    conf_limits = [effect - zcrit * se, effect + zcrit * se]
    z = effect / se
    p_value = 2.0_dp * norm_cdf(-abs(z))
  end subroutine t_test_cluster



  pure elemental logical function all_digits(string) result(ok)
    character(*), intent(in) :: string !! Character value to test; every character must be an ASCII decimal digit.
    integer :: i

    ok = len(string) > 0
    do i = 1, len(string)
      if (string(i:i) < '0' .or. string(i:i) > '9') then
        ok = .false.
        return
      end if
    end do
  end function all_digits

  logical function all_is_numeric(x, extras) result(ok)
    character(*), intent(in) :: x(:) !! Character values to test after trimming leading and trailing spaces.
    character(*), intent(in), optional :: extras(:) !! Values ignored during testing; defaults to '.' and 'NA'.
    character(len=:), allocatable :: s
    real(dp) :: value
    integer :: i, ios, ncandidate
    logical :: skip

    ok = .false.
    ncandidate = 0
    do i = 1, size(x)
      s = trim(adjustl(x(i)))
      skip = len(s) == 0
      if (present(extras)) then
        if (.not. skip) skip = any(s == extras)
      else
        if (.not. skip) skip = s == '.' .or. s == 'NA'
      end if
      if (skip) cycle
      ncandidate = ncandidate + 1
      read(s, *, iostat=ios) value
      if (ios /= 0) return
    end do
    ok = ncandidate > 0
  end function all_is_numeric

  pure subroutine km_quick(time, event, curve_time, curve_surv, status, query_times, query_surv, interval_ge)
    real(dp), intent(in) :: time(:) !! Right-censored follow-up times; NaN times are omitted with their event indicators.
    logical, intent(in) :: event(:) !! Event indicators paired with time; true denotes an observed failure.
    real(dp), allocatable, intent(out) :: curve_time(:) !! Distinct observed event times in ascending order.
    real(dp), allocatable, intent(out) :: curve_surv(:) !! Kaplan-Meier survival immediately after each event time.
    integer, intent(out) :: status !! Zero on success; nonzero for shape errors or an empty retained sample.
    real(dp), intent(in), optional :: query_times(:) !! Times at which to evaluate the Kaplan-Meier step function.
    real(dp), allocatable, intent(out), optional :: query_surv(:) !! Survival values corresponding to query_times.
    logical, intent(in), optional :: interval_ge !! If true, evaluate survival just before events at equal query times.
    real(dp), allocatable :: t(:), u(:), st(:)
    logical, allocatable :: e(:)
    logical :: before_equal
    integer :: i, j, k, n, nevent, nrisk, d
    real(dp) :: surv

    status = 0
    allocate(curve_time(0), curve_surv(0))
    if (present(query_surv)) allocate(query_surv(0))
    if (size(event) /= size(time)) then
      status = 1
      return
    end if
    n = count(.not. ieee_is_nan(time))
    if (n == 0) then
      status = 2
      return
    end if
    allocate(t(n), e(n))
    j = 0
    do i = 1, size(time)
      if (ieee_is_nan(time(i))) cycle
      j = j + 1
      t(j) = time(i)
      e(j) = event(i)
    end do
    call sort_time_event(t, e)
    nevent = 0
    do i = 1, n
      if (.not. e(i)) cycle
      if (i == 1) then
        nevent = nevent + 1
      else if (t(i) /= t(i - 1)) then
        nevent = nevent + 1
      end if
    end do
    allocate(u(nevent), st(nevent))
    surv = 1.0_dp
    k = 0
    i = 1
    do while (i <= n)
      j = i
      do while (j < n)
        if (t(j + 1) /= t(i)) exit
        j = j + 1
      end do
      d = count(e(i:j))
      nrisk = n - i + 1
      if (d > 0) then
        surv = surv * (1.0_dp - real(d, dp) / real(nrisk, dp))
        k = k + 1
        u(k) = t(i)
        st(k) = surv
      end if
      i = j + 1
    end do
    call move_alloc(u, curve_time)
    call move_alloc(st, curve_surv)

    if (present(query_times) .neqv. present(query_surv)) then
      status = 3
      return
    end if
    if (present(query_times)) then
      before_equal = .false.
      if (present(interval_ge)) before_equal = interval_ge
      deallocate(query_surv)
      allocate(query_surv(size(query_times)))
      do i = 1, size(query_times)
        query_surv(i) = 1.0_dp
        do j = 1, size(curve_time)
          if (before_equal) then
            if (curve_time(j) >= query_times(i)) exit
          else
            if (curve_time(j) > query_times(i)) exit
          end if
          query_surv(i) = curve_surv(j)
        end do
      end do
    end if
  end subroutine km_quick

  pure subroutine sort_time_event(time, event)
    real(dp), intent(inout) :: time(:) !! Times to sort in ascending order.
    logical, intent(inout) :: event(:) !! Event flags permuted with time; events precede censoring within exact ties.
    real(dp) :: tv
    logical :: ev
    integer :: i, j

    do i = 2, size(time)
      tv = time(i)
      ev = event(i)
      j = i - 1
      do while (j >= 1)
        if (time(j) < tv) exit
        if (time(j) == tv .and. (event(j) .or. .not. ev)) exit
        time(j + 1) = time(j)
        event(j + 1) = event(j)
        j = j - 1
      end do
      time(j + 1) = tv
      event(j + 1) = ev
    end do
  end subroutine sort_time_event

  pure subroutine bystats_mean(y, group, labels, n_complete, n_missing, means, status)
    real(dp), intent(in) :: y(:,:) !! Numeric response matrix; a row is complete only when all columns are non-NaN.
    integer, intent(in) :: group(:) !! Integer grouping labels paired with rows of y.
    integer, allocatable, intent(out) :: labels(:) !! Sorted distinct group labels followed by no synthetic ALL label.
    integer, allocatable, intent(out) :: n_complete(:) !! Complete-row counts by group, with ALL as the final element.
    integer, allocatable, intent(out) :: n_missing(:) !! Incomplete-row counts by group, with ALL as the final element.
    real(dp), allocatable, intent(out) :: means(:,:) !! Column means by group; final row contains overall means.
    integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or no response columns.
    integer, allocatable :: tmp(:)
    logical, allocatable :: complete(:)
    integer :: i, j, k, ng, nr, nc, countg

    status = 0
    nr = size(y, 1)
    nc = size(y, 2)
    if (size(group) /= nr .or. nc < 1) then
      allocate(labels(0), n_complete(0), n_missing(0), means(0,0))
      status = 1
      return
    end if
    allocate(tmp(max(1, nr)))
    ng = 0
    do i = 1, nr
      if (ng == 0 .or. .not. any(tmp(1:ng) == group(i))) then
        ng = ng + 1
        tmp(ng) = group(i)
      end if
    end do
    if (ng > 1) call sort_int(tmp(1:ng))
    allocate(labels(ng), n_complete(ng + 1), n_missing(ng + 1), means(ng + 1, nc), complete(nr))
    labels = tmp(1:ng)
    complete = .true.
    do i = 1, nr
      complete(i) = .not. any(ieee_is_nan(y(i,:)))
    end do
    means = nan_dp()
    do k = 1, ng
      n_complete(k) = count(complete .and. group == labels(k))
      countg = count(group == labels(k))
      n_missing(k) = countg - n_complete(k)
      if (n_complete(k) > 0) then
        do j = 1, nc
          means(k,j) = sum(y(:,j), mask=complete .and. group == labels(k)) / real(n_complete(k), dp)
        end do
      end if
    end do
    n_complete(ng + 1) = count(complete)
    n_missing(ng + 1) = nr - n_complete(ng + 1)
    if (n_complete(ng + 1) > 0) then
      do j = 1, nc
        means(ng + 1,j) = sum(y(:,j), mask=complete) / real(n_complete(ng + 1), dp)
      end do
    end if
  end subroutine bystats_mean

  pure subroutine bystats2_mean(y, v, h, vlabels, hlabels, n_complete, n_missing, means, status)
    real(dp), intent(in) :: y(:,:) !! Numeric response matrix; rows with any NaN response are incomplete.
    integer, intent(in) :: v(:) !! First integer grouping variable paired with rows of y.
    integer, intent(in) :: h(:) !! Second integer grouping variable paired with rows of y.
    integer, allocatable, intent(out) :: vlabels(:) !! Sorted distinct labels for the first grouping variable.
    integer, allocatable, intent(out) :: hlabels(:) !! Sorted distinct labels for the second grouping variable.
    integer, allocatable, intent(out) :: n_complete(:,:) !! Complete counts; last row/column are ALL margins.
    integer, allocatable, intent(out) :: n_missing(:,:) !! Missing counts; last row/column are ALL margins.
    real(dp), allocatable, intent(out) :: means(:,:,:) !! Means by v, h, response column; last indices are ALL margins.
    integer, intent(out) :: status !! Zero on success; nonzero for incompatible shapes or no response columns.
    integer, allocatable :: tv(:), th(:)
    logical, allocatable :: complete(:), select_row(:)
    integer :: i, j, a, b, nv, nh, nr, nc, ntotal

    status = 0
    nr = size(y, 1)
    nc = size(y, 2)
    if (size(v) /= nr .or. size(h) /= nr .or. nc < 1) then
      allocate(vlabels(0), hlabels(0), n_complete(0,0), n_missing(0,0), means(0,0,0))
      status = 1
      return
    end if
    allocate(tv(max(1,nr)), th(max(1,nr)))
    nv = 0
    nh = 0
    do i = 1, nr
      if (nv == 0 .or. .not. any(tv(1:nv) == v(i))) then
        nv = nv + 1
        tv(nv) = v(i)
      end if
      if (nh == 0 .or. .not. any(th(1:nh) == h(i))) then
        nh = nh + 1
        th(nh) = h(i)
      end if
    end do
    if (nv > 1) call sort_int(tv(1:nv))
    if (nh > 1) call sort_int(th(1:nh))
    allocate(vlabels(nv), hlabels(nh), n_complete(nv+1,nh+1), n_missing(nv+1,nh+1))
    allocate(means(nv+1,nh+1,nc), complete(nr), select_row(nr))
    vlabels = tv(1:nv)
    hlabels = th(1:nh)
    do i = 1, nr
      complete(i) = .not. any(ieee_is_nan(y(i,:)))
    end do
    means = nan_dp()
    do a = 1, nv + 1
      do b = 1, nh + 1
        select_row = .true.
        if (a <= nv) select_row = select_row .and. v == vlabels(a)
        if (b <= nh) select_row = select_row .and. h == hlabels(b)
        ntotal = count(select_row)
        n_complete(a,b) = count(select_row .and. complete)
        n_missing(a,b) = ntotal - n_complete(a,b)
        if (n_complete(a,b) > 0) then
          do j = 1, nc
            means(a,b,j) = sum(y(:,j), mask=select_row .and. complete) / real(n_complete(a,b), dp)
          end do
        end if
      end do
    end do
  end subroutine bystats2_mean


  pure subroutine find_matches(x, y, tol, scale, maxmatch, matches, distance, nmatch, status)
    real(dp), intent(in) :: x(:,:) !! Query observations, one observation per row and matching variable per column.
    real(dp), intent(in) :: y(:,:) !! Candidate observations; columns must match x and rows identify candidate records.
    real(dp), intent(in) :: tol(:) !! Per-variable absolute matching tolerances; nonnegative and length equal to column count.
    real(dp), intent(in) :: scale(:) !! Upstream scaling selector; zero uses unit scale, nonzero uses tol.
    integer, intent(in) :: maxmatch !! Maximum number of candidate matches retained for each query row; must be positive.
    integer, allocatable, intent(out) :: matches(:,:) !! Candidate row indices ordered by increasing scaled squared distance.
    real(dp), allocatable, intent(out) :: distance(:,:) !! Scaled squared distances corresponding to matches; NaN when absent.
    integer, allocatable, intent(out) :: nmatch(:) !! Number of retained matches for each query row, bounded by maxmatch.
    integer, intent(out) :: status !! Zero on success; nonzero for incompatible dimensions or invalid tolerances/maxmatch.
    real(dp), allocatable :: distall(:), scale_eff(:)
    integer, allocatable :: cand(:)
    integer :: i, j, k, n, ny, p, nc, keep, best, tmpi
    real(dp) :: dif, d, bestd, tmpd
    status = 0
    n = size(x,1)
    p = size(x,2)
    ny = size(y,1)
    if (size(y,2) /= p .or. size(tol) /= p .or. size(scale) /= p .or. &
        maxmatch <= 0 .or. any(tol < 0.0_dp)) then
      status = 1
      allocate(matches(0,0), distance(0,0), nmatch(0))
      return
    end if
    allocate(matches(n,maxmatch), distance(n,maxmatch), nmatch(n))
    allocate(distall(ny), cand(ny), scale_eff(p))
    matches = 0
    distance = nan_dp()
    nmatch = 0
    do j = 1, p
      if (scale(j) == 0.0_dp) then
        scale_eff(j) = 1.0_dp
      else
        scale_eff(j) = tol(j)
      end if
    end do
    do i = 1, n
      nc = 0
      do k = 1, ny
        d = 0.0_dp
        do j = 1, p
          dif = abs(y(k,j) - x(i,j))
          if (ieee_is_nan(dif) .or. dif > tol(j)) exit
          if (scale_eff(j) == 0.0_dp) then
            if (dif /= 0.0_dp) exit
          else
            d = d + (dif / scale_eff(j))**2
          end if
        end do
        if (j > p) then
          nc = nc + 1
          cand(nc) = k
          distall(nc) = d
        end if
      end do
      keep = min(nc, maxmatch)
      nmatch(i) = keep
      do j = 1, keep
        best = j
        bestd = distall(j)
        do k = j + 1, nc
          if (distall(k) < bestd) then
            best = k
            bestd = distall(k)
          end if
        end do
        if (best /= j) then
          tmpd = distall(j)
          distall(j) = distall(best)
          distall(best) = tmpd
          tmpi = cand(j)
          cand(j) = cand(best)
          cand(best) = tmpi
        end if
        matches(i,j) = cand(j)
        distance(i,j) = distall(j)
      end do
    end do
  end subroutine find_matches

  pure subroutine seq_freq(positive, assignment, order, obs_per_numcond, status)
    integer, intent(in) :: positive(:,:) !! Observation-by-condition indicators; values must be zero or one.
    integer, allocatable, intent(out) :: assignment(:) !! Selected condition index for each observation, or zero if none.
    integer, allocatable, intent(out) :: order(:) !! Conditions in hierarchical selection order, omitting unused conditions.
    integer, allocatable, intent(out) :: obs_per_numcond(:) !! Counts with 1st element for zero positives, then 1..K positives.
    integer, intent(out) :: status !! Zero on success; nonzero if an indicator is outside zero/one.
    integer, allocatable :: z(:,:), tmporder(:), freq(:), po(:)
    integer :: n, k, i, j, imax, nsel, maxf
    status = 0
    n = size(positive,1)
    k = size(positive,2)
    if (any(positive < 0) .or. any(positive > 1)) then
      status = 1
      allocate(assignment(0), order(0), obs_per_numcond(0))
      return
    end if
    allocate(assignment(n), tmporder(k), obs_per_numcond(k+1), z(n,k), freq(k), po(n))
    assignment = 0
    tmporder = 0
    obs_per_numcond = 0
    z = positive
    do i = 1, n
      po(i) = sum(z(i,:))
      obs_per_numcond(po(i)+1) = obs_per_numcond(po(i)+1) + 1
    end do
    nsel = 0
    do i = 1, k
      do j = 1, k
        freq(j) = sum(z(:,j))
      end do
      maxf = maxval(freq)
      if (maxf == 0) exit
      imax = 1
      do j = 2, k
        if (freq(j) > freq(imax)) imax = j
      end do
      nsel = nsel + 1
      tmporder(nsel) = imax
      do j = 1, n
        if (assignment(j) == 0 .and. z(j,imax) == 1) assignment(j) = imax
        if (z(j,imax) == 1) z(j,:) = 0
      end do
    end do
    allocate(order(nsel))
    if (nsel > 0) order = tmporder(1:nsel)
  end subroutine seq_freq

  pure subroutine pair_up_diff(x, major, minor, group, refgroup, pair_major, pair_minor, &
                               difference, midpoint, sd_difference, lower_difference, &
                               upper_difference, lower_midpoint, upper_midpoint, status, &
                               lower, upper, minkeep, conf_int)
    real(dp), intent(in) :: x(:) !! Numeric estimates to pair across the two group values.
    integer, intent(in) :: major(:) !! Major-category identifiers corresponding one-to-one with x.
    integer, intent(in) :: minor(:) !! Minor-category identifiers corresponding one-to-one with x.
    integer, intent(in) :: group(:) !! Two-valued group identifiers corresponding one-to-one with x.
    integer, intent(in) :: refgroup !! Group identifier subtracted from the alternate group.
    integer, allocatable, intent(out) :: pair_major(:) !! Major identifier for each retained unique major/minor pair.
    integer, allocatable, intent(out) :: pair_minor(:) !! Minor identifier for each retained unique major/minor pair.
    real(dp), allocatable, intent(out) :: difference(:) !! Alternate-minus-reference estimate for each retained pair.
    real(dp), allocatable, intent(out) :: midpoint(:) !! Midpoint of alternate and reference estimates for each pair.
    real(dp), allocatable, intent(out) :: sd_difference(:) !! Approximate standard deviation of each difference, or NaN.
    real(dp), allocatable, intent(out) :: lower_difference(:) !! Lower confidence bound for each difference, or NaN.
    real(dp), allocatable, intent(out) :: upper_difference(:) !! Upper confidence bound for each difference, or NaN.
    real(dp), allocatable, intent(out) :: lower_midpoint(:) !! Midpoint minus half CI width, or NaN when limits absent.
    real(dp), allocatable, intent(out) :: upper_midpoint(:) !! Midpoint plus half CI width, or NaN when limits absent.
    integer, intent(out) :: status !! Zero on success; nonzero for shape, group, or confidence-level errors.
    real(dp), intent(in), optional :: lower(:) !! Lower confidence limits for x; must be supplied together with upper.
    real(dp), intent(in), optional :: upper(:) !! Upper confidence limits for x; must be supplied together with lower.
    real(dp), intent(in), optional :: minkeep !! Drop a complete pair when both estimates are below this threshold.
    real(dp), intent(in), optional :: conf_int !! Confidence level used for supplied limits; defaults to 0.95.
    integer, allocatable :: umajor(:), uminor(:), keep(:), glev(:), tmj(:), tmn(:)
    real(dp), allocatable :: xa(:), xb(:), sda(:), sdb(:), td(:), tm(:), tsd(:), tl(:), tu(:), tlm(:), tum(:)
    integer :: n, i, j, nu, ng, altgroup, nkeep, a, b, best
    real(dp) :: conf, zcrit, tmp
    logical :: lowup
    status = 0
    n = size(x)
    if (size(major) /= n .or. size(minor) /= n .or. size(group) /= n) then
      status = 1
      goto 900
    end if
    lowup = present(lower) .and. present(upper)
    if (present(lower) .neqv. present(upper)) then
      status = 2
      goto 900
    end if
    if (lowup) then
      if (size(lower) /= n .or. size(upper) /= n) then
        status = 1
        goto 900
      end if
    end if
    conf = 0.95_dp
    if (present(conf_int)) conf = conf_int
    if (conf <= 0.0_dp .or. conf >= 1.0_dp) then
      status = 3
      goto 900
    end if
    allocate(glev(2))
    ng = 0
    do i = 1, n
      if (ng == 0 .or. .not. any(glev(1:ng) == group(i))) then
        ng = ng + 1
        if (ng > 2) then
          status = 4
          goto 900
        end if
        glev(ng) = group(i)
      end if
    end do
    if (ng /= 2 .or. .not. any(glev == refgroup)) then
      status = 4
      goto 900
    end if
    if (glev(1) == refgroup) then
      altgroup = glev(2)
    else
      altgroup = glev(1)
    end if
    allocate(umajor(n), uminor(n), xa(n), xb(n), sda(n), sdb(n), keep(n))
    umajor = 0
    uminor = 0
    xa = nan_dp()
    xb = nan_dp()
    sda = nan_dp()
    sdb = nan_dp()
    nu = 0
    zcrit = inv_norm_cdf((1.0_dp + conf) / 2.0_dp)
    do i = 1, n
      j = 0
      do a = 1, nu
        if (umajor(a) == major(i) .and. uminor(a) == minor(i)) then
          j = a
          exit
        end if
      end do
      if (j == 0) then
        nu = nu + 1
        j = nu
        umajor(j) = major(i)
        uminor(j) = minor(i)
      end if
      if (group(i) == refgroup) then
        xa(j) = x(i)
        if (lowup) sda(j) = 0.5_dp * (upper(i) - lower(i)) / zcrit
      else if (group(i) == altgroup) then
        xb(j) = x(i)
        if (lowup) sdb(j) = 0.5_dp * (upper(i) - lower(i)) / zcrit
      end if
    end do
    nkeep = 0
    do j = 1, nu
      if (present(minkeep)) then
        if (.not. ieee_is_nan(xa(j)) .and. .not. ieee_is_nan(xb(j))) then
          if (xa(j) < minkeep .and. xb(j) < minkeep) cycle
        end if
      end if
      nkeep = nkeep + 1
      keep(nkeep) = j
    end do
    allocate(tmj(nkeep), tmn(nkeep), td(nkeep), tm(nkeep), tsd(nkeep), tl(nkeep), tu(nkeep), tlm(nkeep), tum(nkeep))
    do i = 1, nkeep
      j = keep(i)
      tmj(i) = umajor(j)
      tmn(i) = uminor(j)
      td(i) = xb(j) - xa(j)
      tm(i) = 0.5_dp * (xa(j) + xb(j))
      if (lowup) then
        tsd(i) = sqrt(sda(j)**2 + sdb(j)**2)
        tl(i) = td(i) - zcrit * tsd(i)
        tu(i) = td(i) + zcrit * tsd(i)
        tlm(i) = tm(i) - 0.5_dp * zcrit * tsd(i)
        tum(i) = tm(i) + 0.5_dp * zcrit * tsd(i)
      else
        tsd(i) = nan_dp()
        tl(i) = nan_dp()
        tu(i) = nan_dp()
        tlm(i) = nan_dp()
        tum(i) = nan_dp()
      end if
    end do
    do i = 1, nkeep - 1
      best = i
      do j = i + 1, nkeep
        if (tmj(j) < tmj(best) .or. (tmj(j) == tmj(best) .and. td(j) > td(best))) best = j
      end do
      if (best /= i) then
        a = tmj(i)
        tmj(i) = tmj(best)
        tmj(best) = a
        b = tmn(i)
        tmn(i) = tmn(best)
        tmn(best) = b
        tmp = td(i)
        td(i) = td(best)
        td(best) = tmp
        tmp = tm(i)
        tm(i) = tm(best)
        tm(best) = tmp
        tmp = tsd(i)
        tsd(i) = tsd(best)
        tsd(best) = tmp
        tmp = tl(i)
        tl(i) = tl(best)
        tl(best) = tmp
        tmp = tu(i)
        tu(i) = tu(best)
        tu(best) = tmp
        tmp = tlm(i)
        tlm(i) = tlm(best)
        tlm(best) = tmp
        tmp = tum(i)
        tum(i) = tum(best)
        tum(best) = tmp
      end if
    end do
    call move_alloc(tmj, pair_major)
    call move_alloc(tmn, pair_minor)
    call move_alloc(td, difference)
    call move_alloc(tm, midpoint)
    call move_alloc(tsd, sd_difference)
    call move_alloc(tl, lower_difference)
    call move_alloc(tu, upper_difference)
    call move_alloc(tlm, lower_midpoint)
    call move_alloc(tum, upper_midpoint)
    return
900 continue
    allocate(pair_major(0), pair_minor(0), difference(0), midpoint(0), sd_difference(0))
    allocate(lower_difference(0), upper_difference(0), lower_midpoint(0), upper_midpoint(0))
  end subroutine pair_up_diff


  pure subroutine smean_cl_boot(x, bootstrap_index, conf_int, mean_value, lower, upper, reps, status)
    real(dp), intent(in) :: x(:) !! Numeric observations; NaNs are omitted before applying bootstrap indices.
    integer, intent(in) :: bootstrap_index(:,:) !! Resample indices with shape (n_complete,B), each in 1..n_complete.
    real(dp), intent(in) :: conf_int !! Central confidence level strictly between zero and one.
    real(dp), intent(out) :: mean_value !! Mean of the complete observations.
    real(dp), intent(out) :: lower !! Lower percentile bootstrap confidence limit.
    real(dp), intent(out) :: upper !! Upper percentile bootstrap confidence limit.
    real(dp), allocatable, intent(out) :: reps(:) !! Bootstrap replicate means, one per column of bootstrap_index.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, confidence level, or indices.
    real(dp), allocatable :: xx(:), work(:)
    integer :: i, j, n, b
    status = 0
    mean_value = nan_dp()
    lower = nan_dp()
    upper = nan_dp()
    n = count(.not. ieee_is_nan(x))
    b = size(bootstrap_index, 2)
    allocate(reps(b))
    reps = nan_dp()
    if (conf_int <= 0.0_dp .or. conf_int >= 1.0_dp) then
      status = 1
      return
    end if
    if (size(bootstrap_index, 1) /= n .or. b < 1) then
      status = 2
      return
    end if
    allocate(xx(n))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i))) then
        j = j + 1
        xx(j) = x(i)
      end if
    end do
    if (n == 0) return
    mean_value = sum(xx) / real(n, dp)
    if (n < 2) return
    do j = 1, b
      if (any(bootstrap_index(:,j) < 1) .or. any(bootstrap_index(:,j) > n)) then
        status = 3
        return
      end if
      reps(j) = sum(xx(bootstrap_index(:,j))) / real(n, dp)
    end do
    allocate(work(b))
    work = reps
    call sort_real(work)
    lower = quantile_type7_sorted(work, (1.0_dp-conf_int)/2.0_dp)
    upper = quantile_type7_sorted(work, (1.0_dp+conf_int)/2.0_dp)
  end subroutine smean_cl_boot

  pure subroutine smearing_est_tabulated(trans_est, inverse_x, inverse_y, residuals, mode, q, estimate, status)
    real(dp), intent(in) :: trans_est(:) !! Estimated values on the transformed scale.
    real(dp), intent(in) :: inverse_x(:) !! Transformed-scale abscissae defining the tabulated inverse transformation.
    real(dp), intent(in) :: inverse_y(:) !! Original-scale ordinates paired with inverse_x.
    real(dp), intent(in) :: residuals(:) !! Residuals on the transformed scale; NaNs are omitted for mean/quantile modes.
    integer, intent(in) :: mode !! 1=mean smearing, 2=quantile smearing, 3=fitted inverse, 4=unchanged linear predictor.
    real(dp), intent(in) :: q !! Residual quantile probability for mode 2; ignored for other modes.
    real(dp), intent(out) :: estimate(:) !! Original-scale or transformed-scale estimates corresponding to trans_est.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, mode, q, or insufficient data.
    real(dp), allocatable :: rr(:), query(:), vals(:)
    real(dp) :: rq
    integer :: i, j, nr
    status = 0
    estimate = nan_dp()
    if (size(estimate) /= size(trans_est) .or. size(inverse_x) /= size(inverse_y)) then
      status = 1
      return
    end if
    if (mode == 4) then
      estimate = trans_est
      return
    end if
    if (size(inverse_x) < 2) then
      status = 2
      return
    end if
    if (mode == 3) then
      call approx_extrap(inverse_x, inverse_y, trans_est, estimate)
      return
    end if
    nr = count(.not. ieee_is_nan(residuals))
    if (nr < 1) then
      status = 3
      return
    end if
    allocate(rr(nr))
    j = 0
    do i = 1, size(residuals)
      if (.not. ieee_is_nan(residuals(i))) then
        j = j + 1
        rr(j) = residuals(i)
      end if
    end do
    if (mode == 2) then
      if (q < 0.0_dp .or. q > 1.0_dp) then
        status = 4
        return
      end if
      call sort_real(rr)
      rq = quantile_type7_sorted(rr, q)
      allocate(query(size(trans_est)))
      query = trans_est + rq
      call approx_extrap(inverse_x, inverse_y, query, estimate)
    else if (mode == 1) then
      allocate(query(nr), vals(nr))
      do i = 1, size(trans_est)
        query = trans_est(i) + rr
        call approx_extrap(inverse_x, inverse_y, query, vals)
        estimate(i) = sum(vals) / real(nr, dp)
      end do
    else
      status = 5
    end if
  end subroutine smearing_est_tabulated


  pure subroutine bpower_sim_counts(d1, d2, n1, n2, alpha, power, lower, upper, status)
    integer, intent(in) :: d1(:) !! Simulated event counts for group 1, one value per replicate in the range 0..n1.
    integer, intent(in) :: d2(:) !! Simulated event counts for group 2, one value per replicate in the range 0..n2.
    integer, intent(in) :: n1 !! Number of Bernoulli observations in group 1 for every replicate; must be positive.
    integer, intent(in) :: n2 !! Number of Bernoulli observations in group 2 for every replicate; must be positive.
    real(dp), intent(in) :: alpha !! Two-sided type-I error probability strictly between zero and one.
    real(dp), intent(out) :: power !! Fraction of supplied replicates whose Pearson 2x2 chi-square exceeds the critical value.
    real(dp), intent(out) :: lower !! Normal-approximation lower confidence bound power - 1.96*SE, matching upstream.
    real(dp), intent(out) :: upper !! Normal-approximation upper confidence bound power + 1.96*SE, matching upstream.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, counts, sample sizes, or alpha.
    real(dp) :: chisq, crit, den, num, se
    integer :: i, nsim, ntot, nexceed
    status = 0
    power = nan_dp()
    lower = nan_dp()
    upper = nan_dp()
    nsim = size(d1)
    if (size(d2) /= nsim .or. nsim < 1) then
      status = 1
      return
    end if
    if (n1 < 1 .or. n2 < 1) then
      status = 2
      return
    end if
    if (any(d1 < 0) .or. any(d1 > n1) .or. any(d2 < 0) .or. any(d2 > n2)) then
      status = 3
      return
    end if
    if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      status = 4
      return
    end if
    crit = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)**2
    ntot = n1 + n2
    nexceed = 0
    do i = 1, nsim
      den = real(d1(i) + d2(i), dp) * real(ntot - d1(i) - d2(i), dp) * &
            real(n1, dp) * real(n2, dp)
      if (den > 0.0_dp) then
        num = real(ntot, dp) * (real(d1(i), dp) * real(n2 - d2(i), dp) - &
              real(n1 - d1(i), dp) * real(d2(i), dp))**2
        chisq = num / den
        if (chisq > crit) nexceed = nexceed + 1
      end if
    end do
    power = real(nexceed, dp) / real(nsim, dp)
    se = sqrt(power * (1.0_dp - power) / real(nsim, dp))
    lower = power - 1.96_dp * se
    upper = power + 1.96_dp * se
  end subroutine bpower_sim_counts

  pure subroutine spower_simulated(control_failure, intervention_failure, censor_time, alpha, power, &
                                   max_failure, max_censor, status)
    real(dp), intent(in) :: control_failure(:,:) !! Control failure times, shape (nc,nsim), already simulated by caller.
    real(dp), intent(in) :: intervention_failure(:,:) !! Intervention failure times, shape (ni,nsim), caller-supplied.
    real(dp), intent(in) :: censor_time(:,:) !! Censoring times, shape (nc+ni,nsim), paired by subject and replicate.
    real(dp), intent(in) :: alpha !! Type-I error probability for the one-df log-rank chi-square test.
    real(dp), intent(out) :: power !! Fraction of supplied simulations exceeding the chi-square critical value.
    real(dp), intent(out) :: max_failure !! Maximum supplied failure time across both treatment groups and simulations.
    real(dp), intent(out) :: max_censor !! Maximum supplied censoring time across all subjects and simulations.
    integer, intent(out) :: status !! Zero on success; 1 invalid shapes/alpha, 2 log-rank failure, 3 censor support warning.
    real(dp), allocatable :: time(:)
    logical, allocatable :: event(:)
    integer, allocatable :: group(:)
    real(dp) :: chisq, hr, crit, y
    integer :: nc, ni, nsim, i, j, lr_status, nexceed
    status = 0
    power = nan_dp()
    max_failure = nan_dp()
    max_censor = nan_dp()
    nc = size(control_failure, 1)
    ni = size(intervention_failure, 1)
    nsim = size(control_failure, 2)
    if (nc < 1 .or. ni < 1 .or. nsim < 1 .or. size(intervention_failure, 2) /= nsim .or. &
        size(censor_time, 1) /= nc + ni .or. size(censor_time, 2) /= nsim .or. &
        alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      status = 1
      return
    end if
    if (any(ieee_is_nan(control_failure)) .or. any(ieee_is_nan(intervention_failure)) .or. &
        any(ieee_is_nan(censor_time))) then
      status = 1
      return
    end if
    allocate(time(nc + ni), event(nc + ni), group(nc + ni))
    group(:nc) = 1
    group(nc + 1:) = 2
    crit = inv_norm_cdf(1.0_dp - alpha / 2.0_dp)**2
    max_failure = max(maxval(control_failure), maxval(intervention_failure))
    max_censor = maxval(censor_time)
    nexceed = 0
    do j = 1, nsim
      do i = 1, nc
        y = control_failure(i,j)
        time(i) = min(y, censor_time(i,j))
        event(i) = y <= censor_time(i,j)
      end do
      do i = 1, ni
        y = intervention_failure(i,j)
        time(nc+i) = min(y, censor_time(nc+i,j))
        event(nc+i) = y <= censor_time(nc+i,j)
      end do
      call logrank(time, event, group, chisq, hr, lr_status)
      if (lr_status /= 0) then
        status = 2
        return
      end if
      if (chisq > crit) nexceed = nexceed + 1
    end do
    power = real(nexceed, dp) / real(nsim, dp)
    if (max_failure < 0.99_dp * max_censor) status = 3
  end subroutine spower_simulated

  pure subroutine inverse_function_all(x, y, query, roots, nroots, status)
    real(dp), intent(in) :: x(:) !! Abscissae defining the sampled function; must be strictly increasing and finite.
    real(dp), intent(in) :: y(:) !! Ordinates paired with x; adjacent equal ordinates are permitted as flat segments.
    real(dp), intent(in) :: query(:) !! Function values to invert; one output row is produced for each query value.
    real(dp), allocatable, intent(out) :: roots(:,:) !! Roots by query and monotone branch; unused entries are NaN.
    integer, allocatable, intent(out) :: nroots(:) !! Number of roots found for each query value.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid input dimensions, order, or NaNs.
    integer, allocatable :: first(:), last(:)
    integer :: i, j, k, m, n, nb, sgn, prev_sgn, start_idx
    real(dp) :: q, y1, y2, root
    logical :: hit
    status = 0
    n = size(x)
    m = size(query)
    if (size(y) /= n .or. n < 2) then
      status = 1
      allocate(roots(m,1), nroots(m))
      roots = nan_dp()
      nroots = 0
      return
    end if
    if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y)) .or. any(ieee_is_nan(query))) then
      status = 2
      allocate(roots(m,1), nroots(m))
      roots = nan_dp()
      nroots = 0
      return
    end if
    do i = 2, n
      if (x(i) <= x(i-1)) then
        status = 3
        allocate(roots(m,1), nroots(m))
        roots = nan_dp()
        nroots = 0
        return
      end if
    end do
    allocate(first(n), last(n))
    nb = 0
    start_idx = 1
    prev_sgn = 0
    do i = 1, n - 1
      if (y(i+1) > y(i)) then
        sgn = 1
      else if (y(i+1) < y(i)) then
        sgn = -1
      else
        sgn = 0
      end if
      if (sgn /= 0) then
        if (prev_sgn /= 0 .and. sgn /= prev_sgn) then
          nb = nb + 1
          first(nb) = start_idx
          last(nb) = i
          start_idx = i
        end if
        prev_sgn = sgn
      end if
    end do
    nb = nb + 1
    first(nb) = start_idx
    last(nb) = n
    allocate(roots(m,nb), nroots(m))
    roots = nan_dp()
    nroots = 0
    do j = 1, m
      q = query(j)
      do k = 1, nb
        hit = .false.
        do i = first(k), last(k) - 1
          y1 = y(i)
          y2 = y(i+1)
          if (q == y1) then
            root = x(i)
            hit = .true.
            exit
          end if
          if ((q > min(y1,y2) .and. q < max(y1,y2)) .or. q == y2) then
            if (y2 == y1) then
              root = x(i)
            else
              root = x(i) + (q-y1) * (x(i+1)-x(i)) / (y2-y1)
            end if
            hit = .true.
            exit
          end if
        end do
        if (hit) then
          nroots(j) = nroots(j) + 1
          roots(j,nroots(j)) = root
        end if
      end do
      if (nroots(j) == 0) then
        if (abs(q-y(1)) <= abs(q-y(n))) then
          roots(j,1) = x(1)
        else
          roots(j,1) = x(n)
        end if
        nroots(j) = 1
      end if
    end do
  end subroutine inverse_function_all


  pure subroutine bezier_curve(x, y, evaluation, xout, yout, status)
    real(dp), intent(in) :: x(:) !! Control-point x coordinates; must have the same positive length as y.
    real(dp), intent(in) :: y(:) !! Control-point y coordinates corresponding one-to-one with x.
    integer, intent(in) :: evaluation !! Number of equally spaced parameter values on [0,1]; must be at least two.
    real(dp), allocatable, intent(out) :: xout(:) !! Evaluated Bezier-curve x coordinates from first through last control point.
    real(dp), allocatable, intent(out) :: yout(:) !! Evaluated Bezier-curve y coordinates corresponding to xout.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or evaluation count.
    integer :: i, j, n
    real(dp) :: z, bx, by, basis, one_minus

    status = 0
    n = size(x)
    if (n < 1 .or. size(y) /= n .or. evaluation < 2) then
      status = 1
      allocate(xout(0), yout(0))
      return
    end if
    allocate(xout(evaluation), yout(evaluation))
    do i = 1, evaluation
      if (evaluation == 1) then
        z = 0.0_dp
      else
        z = real(i - 1, dp) / real(evaluation - 1, dp)
      end if
      if (i == 1) then
        xout(i) = x(1)
        yout(i) = y(1)
      else if (i == evaluation) then
        xout(i) = x(n)
        yout(i) = y(n)
      else
        one_minus = 1.0_dp - z
        basis = one_minus ** (n - 1)
        bx = 0.0_dp
        by = 0.0_dp
        do j = 0, n - 1
          bx = bx + basis * x(j + 1)
          by = by + basis * y(j + 1)
          if (j < n - 1) then
            basis = basis * real(n - 1 - j, dp) / real(j + 1, dp) * z / one_minus
          end if
        end do
        xout(i) = bx
        yout(i) = by
      end if
    end do
  end subroutine bezier_curve


  pure subroutine cut2_explicit(x, cuts, minmax, groups, effective_cuts, status)
    real(dp), intent(in) :: x(:) !! Numeric observations to categorize; NaNs receive group zero.
    real(dp), intent(in) :: cuts(:) !! Strictly increasing requested cut points defining left-inclusive intervals.
    logical, intent(in) :: minmax !! If true, extend cuts to the observed finite minimum and maximum when needed.
    integer, allocatable, intent(out) :: groups(:) !! Group indices 1..K-1; zero denotes missing or outside supplied cuts.
    real(dp), allocatable, intent(out) :: effective_cuts(:) !! Cut points after optional observed-range extension.
    integer, intent(out) :: status !! Zero on success; nonzero for insufficient/unsorted cuts or no finite observations.
    real(dp), allocatable :: tmp(:)
    real(dp) :: xmin, xmax
    integer :: i, j, ncut, finite_n

    status = 0
    finite_n = count(.not. ieee_is_nan(x))
    if (size(cuts) < 2) then
      status = 1
      allocate(groups(size(x)), effective_cuts(0))
      groups = 0
      return
    end if
    do i = 2, size(cuts)
      if (cuts(i) <= cuts(i - 1) .or. ieee_is_nan(cuts(i))) then
        status = 2
        allocate(groups(size(x)), effective_cuts(0))
        groups = 0
        return
      end if
    end do
    if (ieee_is_nan(cuts(1))) then
      status = 2
      allocate(groups(size(x)), effective_cuts(0))
      groups = 0
      return
    end if
    if (finite_n == 0) then
      status = 3
      allocate(groups(size(x)), effective_cuts(0))
      groups = 0
      return
    end if

    xmin = minval(x, mask=.not. ieee_is_nan(x))
    xmax = maxval(x, mask=.not. ieee_is_nan(x))
    ncut = size(cuts)
    if (minmax) then
      allocate(tmp(ncut + merge(1, 0, xmin < cuts(1)) + merge(1, 0, xmax > cuts(ncut))))
      j = 0
      if (xmin < cuts(1)) then
        j = j + 1
        tmp(j) = xmin
      end if
      tmp(j + 1:j + ncut) = cuts
      j = j + ncut
      if (xmax > cuts(ncut)) then
        j = j + 1
        tmp(j) = xmax
      end if
      effective_cuts = tmp(:j)
    else
      effective_cuts = cuts
    end if

    allocate(groups(size(x)))
    groups = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i))) cycle
      do j = 1, size(effective_cuts) - 1
        if (j < size(effective_cuts) - 1) then
          if (x(i) >= effective_cuts(j) .and. x(i) < effective_cuts(j + 1)) then
            groups(i) = j
            exit
          end if
        else
          if (x(i) >= effective_cuts(j) .and. x(i) <= effective_cuts(j + 1)) then
            groups(i) = j
            exit
          end if
        end if
      end do
    end do
  end subroutine cut2_explicit


  pure subroutine sim_po_cuts(n, p, odds_ratio, uniforms0, uniforms1, odds_ratios, status)
    integer, intent(in) :: n !! Total simulated sample size; must be positive and even, split equally between groups.
    real(dp), intent(in) :: p(:) !! Baseline category probabilities; nonnegative and summing to one.
    real(dp), intent(in) :: odds_ratio !! Proportional-odds ratio used to transform baseline probabilities.
    real(dp), intent(in) :: uniforms0(:,:) !! Baseline-group uniforms in [0,1], shape (n/2, nsim).
    real(dp), intent(in) :: uniforms1(:,:) !! Treatment-group uniforms in [0,1], same shape as uniforms0.
    real(dp), allocatable, intent(out) :: odds_ratios(:,:) !! Simulated cumulative odds ratios, shape (nsim, K-1).
    integer, intent(out) :: status !! Zero on success; nonzero for invalid sizes, probabilities, or uniforms.
    real(dp), allocatable :: p1(:), cum0(:), cum1(:)
    integer, allocatable :: y0(:), y1(:)
    integer :: half, i, j, k, nsim, cat
    real(dp) :: prop0, prop1, odds0, odds1

    status = 0
    half = n / 2
    nsim = size(uniforms0, 2)
    if (n <= 0 .or. mod(n, 2) /= 0 .or. size(p) < 2) then
      status = 1
      allocate(odds_ratios(0,0))
      return
    end if
    if (size(uniforms0, 1) /= half .or. size(uniforms1, 1) /= half .or. &
        size(uniforms1, 2) /= nsim) then
      status = 2
      allocate(odds_ratios(0,0))
      return
    end if
    if (any(p < 0.0_dp) .or. abs(sum(p) - 1.0_dp) > 1.0e-5_dp) then
      status = 3
      allocate(odds_ratios(0,0))
      return
    end if
    if (any(uniforms0 < 0.0_dp) .or. any(uniforms0 > 1.0_dp) .or. &
        any(uniforms1 < 0.0_dp) .or. any(uniforms1 > 1.0_dp)) then
      status = 4
      allocate(odds_ratios(0,0))
      return
    end if

    allocate(p1(size(p)), cum0(size(p)), cum1(size(p)))
    call pomodm(p, odds_ratio, p1)
    cum0(1) = p(1)
    cum1(1) = p1(1)
    do k = 2, size(p)
      cum0(k) = cum0(k - 1) + p(k)
      cum1(k) = cum1(k - 1) + p1(k)
    end do
    allocate(y0(half), y1(half), odds_ratios(nsim, size(p) - 1))
    odds_ratios = nan_dp()
    do i = 1, nsim
      do j = 1, half
        y0(j) = size(p)
        y1(j) = size(p)
        do cat = 1, size(p)
          if (uniforms0(j,i) <= cum0(cat)) then
            y0(j) = cat
            exit
          end if
        end do
        do cat = 1, size(p)
          if (uniforms1(j,i) <= cum1(cat)) then
            y1(j) = cat
            exit
          end if
        end do
      end do
      do k = 2, size(p)
        prop0 = real(count(y0 >= k), dp) / real(half, dp)
        prop1 = real(count(y1 >= k), dp) / real(half, dp)
        if (prop0 <= 0.0_dp .or. prop0 >= 1.0_dp .or. &
            prop1 <= 0.0_dp .or. prop1 >= 1.0_dp) cycle
        odds0 = prop0 / (1.0_dp - prop0)
        odds1 = prop1 / (1.0_dp - prop1)
        odds_ratios(i,k - 1) = odds1 / odds0
      end do
    end do
  end subroutine sim_po_cuts


  pure subroutine rcspline_function(x, knots, coef, norm, integral, y, status)
    real(dp), intent(in) :: x(:) !! Evaluation points for the restricted cubic spline function.
    real(dp), intent(in) :: knots(:) !! Strictly increasing spline knots; at least three are required.
    real(dp), intent(in) :: coef(:) !! Coefficients: K values include intercept, or K-1 values omit intercept.
    integer, intent(in) :: norm !! Normalization mode: 0 none, 1 last-knot spacing, otherwise overall-range two-thirds power.
    logical, intent(in) :: integral !! If true, evaluate the integrated restricted cubic spline form.
    real(dp), allocatable, intent(out) :: y(:) !! Spline values corresponding to x.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid knots, coefficient length, or normalization.
    real(dp), allocatable :: c(:)
    real(dp) :: kd, term, t1, t2, t3
    integer :: i, j, k

    status = 0
    k = size(knots)
    if (k < 3) then
      status = 1
      allocate(y(0))
      return
    end if
    do i = 2, k
      if (knots(i) <= knots(i - 1)) then
        status = 2
        allocate(y(0))
        return
      end if
    end do
    if (size(coef) == k) then
      c = coef
    else if (size(coef) == k - 1) then
      allocate(c(k))
      c(1) = 0.0_dp
      c(2:) = coef
    else
      status = 3
      allocate(y(0))
      return
    end if
    select case (norm)
    case (0)
      kd = 1.0_dp
    case (1)
      kd = knots(k) - knots(k - 1)
    case default
      kd = (knots(k) - knots(1)) ** (2.0_dp / 3.0_dp)
    end select
    if (kd <= 0.0_dp) then
      status = 4
      allocate(y(0))
      return
    end if

    allocate(y(size(x)))
    if (.not. integral) then
      y = c(1) + c(2) * x
      do j = 1, k - 2
        do i = 1, size(x)
          t1 = max((x(i) - knots(j)) / kd, 0.0_dp) ** 3
          t2 = max((x(i) - knots(k)) / kd, 0.0_dp) ** 3
          t3 = max((x(i) - knots(k - 1)) / kd, 0.0_dp) ** 3
          term = t1 + ((knots(k - 1) - knots(j)) * t2 - &
                 (knots(k) - knots(j)) * t3) / (knots(k) - knots(k - 1))
          y(i) = y(i) + c(j + 2) * term
        end do
      end do
    else
      y = c(1) * x + 0.5_dp * c(2) * x * x
      do j = 1, k - 2
        do i = 1, size(x)
          t1 = max((x(i) - knots(j)) / kd, 0.0_dp) ** 4
          t2 = max((x(i) - knots(k)) / kd, 0.0_dp) ** 4
          t3 = max((x(i) - knots(k - 1)) / kd, 0.0_dp) ** 4
          term = t1 + ((knots(k - 1) - knots(j)) * t2 - &
                 (knots(k) - knots(j)) * t3) / (knots(k) - knots(k - 1))
          y(i) = y(i) + 0.25_dp * c(j + 2) * kd * term
        end do
      end do
    end if
  end subroutine rcspline_function


  pure subroutine rcspline_restate_coefficients(knots, coef, norm, integral, expanded, status)
    real(dp), intent(in) :: knots(:) !! Strictly increasing restricted-cubic-spline knots; at least three required.
    real(dp), intent(in) :: coef(:) !! K coefficients include intercept, or K-1 coefficients omit it.
    integer, intent(in) :: norm !! Normalization mode matching rcspline.restate: 0 none, 1 last interval cubed, 2 range squared.
    logical, intent(in) :: integral !! If true, transform coefficients for the integrated truncated-power representation.
    real(dp), allocatable, intent(out) :: expanded(:) !! Numeric truncated-power coefficients, with intercept when supplied.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid knot or coefficient dimensions.
    real(dp), allocatable :: work(:), core(:), dvec(:)
    real(dp) :: intercept, kd, coefk, coefk1
    integer :: i, k, p, offset
    logical :: has_intercept

    status = 0
    k = size(knots)
    if (k < 3) then
      status = 1
      allocate(expanded(0))
      return
    end if
    do i = 2, k
      if (knots(i) <= knots(i - 1)) then
        status = 2
        allocate(expanded(0))
        return
      end if
    end do
    has_intercept = size(coef) == k
    if (.not. has_intercept .and. size(coef) /= k - 1) then
      status = 3
      allocate(expanded(0))
      return
    end if
    if (has_intercept) then
      intercept = coef(1)
      allocate(work(k - 1))
      work = coef(2:)
    else
      intercept = 0.0_dp
      allocate(work(k - 1))
      work = coef
    end if
    p = size(work)
    select case (norm)
    case (0)
      kd = 1.0_dp
    case (1)
      kd = (knots(k) - knots(k - 1)) ** 3
    case default
      kd = (knots(k) - knots(1)) ** 2
    end select
    if (p > 1) work(2:) = work(2:) / kd
    allocate(dvec(p))
    dvec = 0.0_dp
    do i = 2, p
      dvec(i) = knots(i - 1) - knots(k)
    end do
    coefk = sum(work * dvec) / (knots(k) - knots(k - 1))
    dvec = 0.0_dp
    do i = 2, p
      dvec(i) = knots(i - 1) - knots(k - 1)
    end do
    coefk1 = sum(work * dvec) / (knots(k - 1) - knots(k))
    allocate(core(p + 2))
    core(1:p) = work
    core(p + 1) = coefk
    core(p + 2) = coefk1
    if (integral) then
      core(1) = 0.5_dp * core(1)
      if (size(core) > 1) core(2:) = 0.25_dp * core(2:)
    end if
    offset = merge(1, 0, has_intercept)
    allocate(expanded(size(core) + offset))
    if (has_intercept) then
      expanded(1) = intercept
      expanded(2:) = core
    else
      expanded = core
    end if
  end subroutine rcspline_restate_coefficients


  pure subroutine spearman_test(x, y, power, rsquare, fstat, df1, df2, p_value, n_complete, status)
    real(dp), intent(in) :: x(:) !! Numeric predictor values; NaN pairs are omitted before ranking.
    real(dp), intent(in) :: y(:) !! Numeric response values paired one-to-one with x; NaN pairs are omitted.
    integer, intent(in) :: power !! Rank-polynomial degree, either 1 or 2 as in Hmisc spearman.test.
    real(dp), intent(out) :: rsquare !! Rank-regression coefficient of determination.
    real(dp), intent(out) :: fstat !! F statistic for testing the rank association.
    real(dp), intent(out) :: df1 !! Numerator degrees of freedom, equal to power.
    real(dp), intent(out) :: df2 !! Denominator degrees of freedom, n - power - 1.
    real(dp), intent(out) :: p_value !! Upper-tail F-test probability.
    integer, intent(out) :: n_complete !! Number of complete x/y pairs used.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shapes, power, or insufficient data.
    real(dp) :: adjusted

    call spearman2_numeric(x, y, power, rsquare, fstat, df1, df2, p_value, &
                           adjusted, n_complete, status)
  end subroutine spearman_test


  pure subroutine chi_square_codes(x, y, statistic, df, excess, p_value, status)
    integer, intent(in) :: x(:) !! Positive integer category codes for the predictor; nonpositive codes are treated as missing.
    integer, intent(in) :: y(:) !! Positive integer category codes for the outcome; nonpositive codes are treated as missing.
    real(dp), intent(out) :: statistic !! Pearson chi-square statistic without continuity correction.
    real(dp), intent(out) :: df !! Pearson chi-square degrees of freedom after empty categories are discarded.
    real(dp), intent(out) :: excess !! Hmisc chiSquare summary statistic chi-square minus degrees of freedom.
    real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
    integer, intent(out) :: status !! Zero on success; nonzero for shape or degenerate-table errors.
    real(dp), allocatable :: tab(:,:)
    integer :: i, nx, ny

    statistic = nan_dp()
    df = nan_dp()
    excess = nan_dp()
    p_value = nan_dp()
    status = 0
    if (size(x) /= size(y) .or. size(x) < 1) then
      status = 1
      return
    end if
    nx = max(0, maxval(x))
    ny = max(0, maxval(y))
    if (nx < 2 .or. ny < 2) then
      status = 2
      return
    end if
    allocate(tab(nx,ny))
    tab = 0.0_dp
    do i = 1, size(x)
      if (x(i) > 0 .and. y(i) > 0) tab(x(i),y(i)) = tab(x(i),y(i)) + 1.0_dp
    end do
    call cat_test_chisq(tab, statistic, df, p_value, status)
    if (status == 0) excess = statistic - df
  end subroutine chi_square_codes


  pure subroutine combine_levels_codes(x, min_count, ordered, combined, map, status)
    integer, intent(in) :: x(:) !! Integer category codes 1..K; nonpositive values are preserved as missing code zero.
    integer, intent(in) :: min_count !! Minimum required frequency for a retained category, corresponding to Hmisc m.
    logical, intent(in) :: ordered !! If true, only adjacent categories are combined; otherwise rare categories are pooled.
    integer, allocatable, intent(out) :: combined(:) !! Re-coded observations with consecutive positive category codes.
    integer, allocatable, intent(out) :: map(:) !! Mapping from each original category 1..K to its new category code.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid minimum count or impossible ordered pooling.
    integer, allocatable :: counts(:), workmap(:), freq(:), uniq(:)
    integer :: i, j, ncat, nactive, low, neighbor, target, tmp, rare, nr

    status = 0
    if (min_count < 1) then
      status = 1
      allocate(combined(0), map(0))
      return
    end if
    ncat = max(0, maxval(x))
    allocate(combined(size(x)), map(ncat))
    combined = 0
    if (ncat == 0) then
      map = 0
      return
    end if
    allocate(counts(ncat), workmap(ncat))
    counts = 0
    do i = 1, size(x)
      if (x(i) > 0 .and. x(i) <= ncat) counts(x(i)) = counts(x(i)) + 1
    end do
    workmap = [(i, i=1,ncat)]
    if (count(counts > 0) < 2) then
      map = workmap
      do i = 1, size(x)
        if (x(i) > 0 .and. x(i) <= ncat) combined(i) = map(x(i))
      end do
      return
    end if

    if (ordered) then
      if (sum(counts) < 2 * min_count) then
        status = 2
        map = workmap
        return
      end if
      do
        allocate(uniq(ncat), freq(ncat))
        uniq = 0
        freq = 0
        nactive = 0
        do i = 1, ncat
          if (counts(i) <= 0) cycle
          if (nactive == 0 .or. all(uniq(1:nactive) /= workmap(i))) then
            nactive = nactive + 1
            uniq(nactive) = workmap(i)
          end if
          do j = 1, nactive
            if (uniq(j) == workmap(i)) then
              freq(j) = freq(j) + counts(i)
              exit
            end if
          end do
        end do
        low = 1
        do j = 2, nactive
          if (freq(j) < freq(low)) low = j
        end do
        if (freq(low) >= min_count .or. nactive < 2) then
          deallocate(uniq, freq)
          exit
        end if
        if (low == 1) then
          neighbor = 2
        else if (low == nactive) then
          neighbor = nactive - 1
        else if (freq(low-1) < freq(low+1)) then
          neighbor = low - 1
        else
          neighbor = low + 1
        end if
        target = min(uniq(low), uniq(neighbor))
        do i = 1, ncat
          if (workmap(i) == uniq(low) .or. workmap(i) == uniq(neighbor)) workmap(i) = target
        end do
        deallocate(uniq, freq)
      end do
    else
      rare = 0
      do i = 1, ncat
        if (counts(i) > 0 .and. counts(i) < min_count) rare = rare + 1
      end do
      if (rare == 1) then
        low = 0
        neighbor = 0
        do i = 1, ncat
          if (counts(i) <= 0) cycle
          if (low == 0 .or. counts(i) < counts(low)) then
            neighbor = low
            low = i
          else if (neighbor == 0 .or. counts(i) < counts(neighbor)) then
            neighbor = i
          end if
        end do
        if (neighbor > 0) workmap(low) = workmap(neighbor)
      else if (rare > 1) then
        target = 0
        do i = 1, ncat
          if (counts(i) > 0 .and. counts(i) < min_count) then
            if (target == 0) target = i
            workmap(i) = target
          end if
        end do
      end if
    end if

    allocate(uniq(ncat))
    uniq = 0
    nr = 0
    do i = 1, ncat
      if (counts(i) <= 0) cycle
      if (nr == 0 .or. all(uniq(1:nr) /= workmap(i))) then
        nr = nr + 1
        uniq(nr) = workmap(i)
      end if
    end do
    do i = 1, nr - 1
      do j = i + 1, nr
        if (uniq(j) < uniq(i)) then
          tmp = uniq(i)
          uniq(i) = uniq(j)
          uniq(j) = tmp
        end if
      end do
    end do
    map = 0
    do i = 1, ncat
      if (counts(i) <= 0) cycle
      do j = 1, nr
        if (workmap(i) == uniq(j)) then
          map(i) = j
          exit
        end if
      end do
    end do
    do i = 1, size(x)
      if (x(i) > 0 .and. x(i) <= ncat) combined(i) = map(x(i))
    end do
  end subroutine combine_levels_codes


  pure subroutine jacobi_eigensystem(a, values, vectors, status)
    real(dp), intent(in) :: a(:,:) !! Real symmetric matrix whose complete eigensystem is required.
    real(dp), allocatable, intent(out) :: values(:) !! Eigenvalues sorted in descending order.
    real(dp), allocatable, intent(out) :: vectors(:,:) !! Corresponding orthonormal eigenvectors in columns.
    integer, intent(out) :: status !! Zero on convergence; nonzero for non-square input or iteration failure.
    real(dp), allocatable :: d(:,:), v(:,:)
    real(dp) :: app, aqq, apq, phi, c, sn, dip, diq, vip, viq, off, tmp
    integer :: n, i, j, p, q, sweep, maxi

    n = size(a,1)
    status = 0
    if (size(a,2) /= n .or. n < 1) then
      status = 1
      allocate(values(0), vectors(0,0))
      return
    end if
    allocate(d(n,n), v(n,n), values(n), vectors(n,n))
    d = a
    v = 0.0_dp
    do i = 1, n
      v(i,i) = 1.0_dp
    end do
    do sweep = 1, max(50, 30*n*n)
      off = 0.0_dp
      p = 1
      q = min(2,n)
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(d(i,j)) > off) then
            off = abs(d(i,j))
            p = i
            q = j
          end if
        end do
      end do
      if (off <= 1.0e-13_dp * max(1.0_dp, maxval(abs(d)))) exit
      app = d(p,p)
      aqq = d(q,q)
      apq = d(p,q)
      phi = 0.5_dp * atan2(2.0_dp*apq, aqq-app)
      c = cos(phi)
      sn = sin(phi)
      do i = 1, n
        if (i == p .or. i == q) cycle
        dip = d(i,p)
        diq = d(i,q)
        d(i,p) = c*dip - sn*diq
        d(p,i) = d(i,p)
        d(i,q) = sn*dip + c*diq
        d(q,i) = d(i,q)
      end do
      d(p,p) = c*c*app - 2.0_dp*c*sn*apq + sn*sn*aqq
      d(q,q) = sn*sn*app + 2.0_dp*c*sn*apq + c*c*aqq
      d(p,q) = 0.0_dp
      d(q,p) = 0.0_dp
      do i = 1, n
        vip = v(i,p)
        viq = v(i,q)
        v(i,p) = c*vip - sn*viq
        v(i,q) = sn*vip + c*viq
      end do
    end do
    if (sweep > max(50, 30*n*n)) then
      status = 2
      values = nan_dp()
      vectors = nan_dp()
      return
    end if
    do i = 1, n
      values(i) = d(i,i)
      vectors(:,i) = v(:,i)
    end do
    do i = 1, n - 1
      maxi = i
      do j = i + 1, n
        if (values(j) > values(maxi)) maxi = j
      end do
      if (maxi /= i) then
        tmp = values(i)
        values(i) = values(maxi)
        values(maxi) = tmp
        do j = 1, n
          tmp = vectors(j,i)
          vectors(j,i) = vectors(j,maxi)
          vectors(j,maxi) = tmp
        end do
      end if
    end do
  end subroutine jacobi_eigensystem


  pure subroutine princmp_regular(x, cor, k, scores, standardized_loadings, &
                                  original_loadings, variances, scale, status)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable numeric matrix; rows containing NaN are omitted for fitting.
    logical, intent(in) :: cor !! If true, standardize variables before PCA as in princmp(method='regular', cor=TRUE).
    integer, intent(in) :: k !! Number of leading components requested, from one through the number of variables.
    real(dp), allocatable, intent(out) :: scores(:,:) !! N-by-k component scores, with NaNs on incomplete input rows.
    real(dp), allocatable, intent(out) :: standardized_loadings(:,:) !! Variable-by-k loadings on the PCA analysis scale.
    real(dp), allocatable, intent(out) :: original_loadings(:,:) !! Variable-by-k loadings on the original variable scale.
    real(dp), allocatable, intent(out) :: variances(:) !! Variances of all principal components in descending order.
    real(dp), allocatable, intent(out) :: scale(:) !! Per-variable scale factors used for standardization, or ones if cor is false.
    integer, intent(out) :: status !! Zero on success; nonzero for shape, constant-variable, or eigensolver errors.
    real(dp), allocatable :: complete(:,:), z(:,:), means(:), cov(:,:), eigvec(:,:), eigval(:)
    integer :: n, p, nc, i, j, row, st

    n = size(x,1)
    p = size(x,2)
    status = 0
    if (p < 1 .or. k < 1 .or. k > p) then
      status = 1
      allocate(scores(0,0), standardized_loadings(0,0), original_loadings(0,0), variances(0), scale(0))
      return
    end if
    nc = 0
    do i = 1, n
      if (.not. any(ieee_is_nan(x(i,:)))) nc = nc + 1
    end do
    if (nc < 2) then
      status = 2
      allocate(scores(0,0), standardized_loadings(0,0), original_loadings(0,0), variances(0), scale(0))
      return
    end if
    allocate(complete(nc,p), z(nc,p), means(p), scale(p))
    row = 0
    do i = 1, n
      if (.not. any(ieee_is_nan(x(i,:)))) then
        row = row + 1
        complete(row,:) = x(i,:)
      end if
    end do
    do j = 1, p
      means(j) = sum(complete(:,j)) / real(nc,dp)
      if (cor) then
        scale(j) = sqrt(sum((complete(:,j)-means(j))**2) / real(nc,dp))
        if (scale(j) <= tiny(1.0_dp)) then
          status = 3
          allocate(scores(0,0), standardized_loadings(0,0), original_loadings(0,0), variances(0))
          return
        end if
      else
        scale(j) = 1.0_dp
      end if
      z(:,j) = (complete(:,j) - means(j)) / scale(j)
    end do
    allocate(cov(p,p))
    cov = matmul(transpose(z), z) / real(nc,dp)
    call jacobi_eigensystem(cov, eigval, eigvec, st)
    if (st /= 0) then
      status = 4
      allocate(scores(0,0), standardized_loadings(0,0), original_loadings(0,0), variances(0))
      return
    end if
    allocate(variances(p), standardized_loadings(p,k), original_loadings(p,k), scores(n,k))
    variances = max(eigval, 0.0_dp)
    standardized_loadings = eigvec(:,1:k)
    do j = 1, k
      original_loadings(:,j) = standardized_loadings(:,j) / scale
    end do
    scores = nan_dp()
    do i = 1, n
      if (.not. any(ieee_is_nan(x(i,:)))) then
        do j = 1, k
          scores(i,j) = sum(((x(i,:)-means)/scale) * standardized_loadings(:,j))
        end do
      end if
    end do
  end subroutine princmp_regular


  pure subroutine bootkm_resampled(time, event, bootstrap_index, use_time, q, eval_time, estimates, status)
    real(dp), intent(in) :: time(:) !! Right-censored follow-up times for the original sample.
    logical, intent(in) :: event(:) !! Event indicators paired with time; true denotes an observed failure.
    integer, intent(in) :: bootstrap_index(:,:) !! N-by-B one-based resampling indices, one bootstrap sample per column.
    logical, intent(in) :: use_time !! If true estimate survival at eval_time; otherwise estimate the q survival quantile.
    real(dp), intent(in) :: q !! Survival probability defining the requested quantile when use_time is false.
    real(dp), intent(in) :: eval_time !! Follow-up time at which survival is evaluated when use_time is true.
    real(dp), allocatable, intent(out) :: estimates(:) !! One bootstrap survival or quantile estimate per resample.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shapes, indices, or Kaplan-Meier failures.
    real(dp), allocatable :: bt(:), ct(:), cs(:), qs(:)
    logical, allocatable :: be(:)
    integer :: b, i, j, n, st

    status = 0
    allocate(estimates(0))
    n = size(time)
    if (size(event) /= n .or. size(bootstrap_index,1) /= n .or. n < 1) then
      status = 1
      return
    end if
    if (.not. use_time .and. (q < 0.0_dp .or. q > 1.0_dp)) then
      status = 2
      return
    end if
    if (any(bootstrap_index < 1) .or. any(bootstrap_index > n)) then
      status = 3
      return
    end if
    deallocate(estimates)
    allocate(estimates(size(bootstrap_index,2)), bt(n), be(n))
    do b = 1, size(bootstrap_index,2)
      do i = 1, n
        j = bootstrap_index(i,b)
        bt(i) = time(j)
        be(i) = event(j)
      end do
      if (use_time) then
        call km_quick(bt, be, ct, cs, st, [eval_time], qs, .false.)
        if (st /= 0) then
          status = 10 + st
          return
        end if
        estimates(b) = qs(1)
      else
        call km_quick(bt, be, ct, cs, st)
        if (st /= 0) then
          status = 10 + st
          return
        end if
        estimates(b) = nan_dp()
        do i = 1, size(cs)
          if (cs(i) <= q) then
            estimates(b) = ct(i)
            exit
          end if
        end do
      end if
    end do
  end subroutine bootkm_resampled

  pure subroutine match_cases_core(xcase, xcontrol, tol, maxmatch, use_closest, random_score, &
                                   matches, nmatch, status)
    real(dp), intent(in) :: xcase(:) !! Scalar matching variable for each case; NaN cases receive no matches.
    real(dp), intent(in) :: xcontrol(:) !! Scalar matching variable for controls; NaN controls are never selected.
    real(dp), intent(in) :: tol !! Maximum absolute case-control difference allowed; must be nonnegative.
    integer, intent(in) :: maxmatch !! Maximum number of controls retained for each case; must be positive.
    logical, intent(in) :: use_closest !! If true retain smallest absolute differences; otherwise rank by random_score.
    real(dp), intent(in) :: random_score(:,:) !! Case-by-control selection scores used when use_closest is false.
    integer, allocatable, intent(out) :: matches(:,:) !! Case-by-maxmatch one-based control indices; zeros are unused slots.
    integer, allocatable, intent(out) :: nmatch(:) !! Number of retained controls for each case.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid limits or random-score dimensions.
    integer, allocatable :: candidate(:)
    real(dp), allocatable :: key(:)
    integer :: i, j, k, m, nc, tmpi
    real(dp) :: tmpk

    status = 0
    allocate(matches(0,0), nmatch(0))
    if (tol < 0.0_dp .or. maxmatch < 1) then
      status = 1
      return
    end if
    if (.not. use_closest) then
      if (size(random_score,1) /= size(xcase) .or. size(random_score,2) /= size(xcontrol)) then
        status = 2
        return
      end if
    end if
    deallocate(matches, nmatch)
    allocate(matches(size(xcase),maxmatch), nmatch(size(xcase)))
    matches = 0
    nmatch = 0
    allocate(candidate(size(xcontrol)), key(size(xcontrol)))
    do i = 1, size(xcase)
      if (ieee_is_nan(xcase(i))) cycle
      nc = 0
      do j = 1, size(xcontrol)
        if (ieee_is_nan(xcontrol(j))) cycle
        if (abs(xcontrol(j) - xcase(i)) <= tol) then
          nc = nc + 1
          candidate(nc) = j
          if (use_closest) then
            key(nc) = abs(xcontrol(j) - xcase(i))
          else
            key(nc) = random_score(i,j)
          end if
        end if
      end do
      do j = 2, nc
        tmpi = candidate(j)
        tmpk = key(j)
        k = j - 1
        do while (k >= 1)
          if (key(k) <= tmpk) exit
          candidate(k+1) = candidate(k)
          key(k+1) = key(k)
          k = k - 1
        end do
        candidate(k+1) = tmpi
        key(k+1) = tmpk
      end do
      m = min(nc, maxmatch)
      nmatch(i) = m
      if (m > 0) matches(i,1:m) = candidate(1:m)
    end do
  end subroutine match_cases_core

  pure subroutine largest_empty_rexhaustive(x, y, xlim, ylim, width, height, center, rect, area, status)
    real(dp), intent(in) :: x(:) !! Point x coordinates; pairs containing NaN are omitted.
    real(dp), intent(in) :: y(:) !! Point y coordinates; pairs containing NaN are omitted.
    real(dp), intent(in) :: xlim(2) !! Lower and upper horizontal bounds of the search rectangle.
    real(dp), intent(in) :: ylim(2) !! Lower and upper vertical bounds of the search rectangle.
    real(dp), intent(in) :: width !! Strict minimum accepted rectangle width, matching largest.empty rexhaustive.
    real(dp), intent(in) :: height !! Strict minimum accepted rectangle height, matching largest.empty rexhaustive.
    real(dp), intent(out) :: center(2) !! Center coordinates of the selected empty rectangle, or NaNs on failure.
    real(dp), intent(out) :: rect(4) !! Rectangle bounds as xmin, ymin, xmax, ymax, or NaNs on failure.
    real(dp), intent(out) :: area !! Area of the selected rectangle, or NaN when no acceptable rectangle exists.
    integer, intent(out) :: status !! Zero on success; nonzero for shapes/bounds or no rectangle satisfying limits.
    real(dp), allocatable :: xx(:), yy(:), d(:)
    real(dp) :: maxr, trial, tl, tr, ri, li, tx, ty
    integer :: i, j, k, n, m, nkeep

    status = 0
    center = nan_dp()
    rect = nan_dp()
    area = nan_dp()
    if (size(x) /= size(y) .or. xlim(2) <= xlim(1) .or. ylim(2) <= ylim(1) .or. &
        width < 0.0_dp .or. height < 0.0_dp) then
      status = 1
      return
    end if
    nkeep = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    allocate(xx(nkeep), yy(nkeep))
    k = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
      k = k + 1
      xx(k) = x(i)
      yy(k) = y(i)
    end do
    n = nkeep
    do i = 2, n
      tx = xx(i)
      ty = yy(i)
      j = i - 1
      do while (j >= 1)
        if (yy(j) <= ty) exit
        xx(j+1) = xx(j)
        yy(j+1) = yy(j)
        j = j - 1
      end do
      xx(j+1) = tx
      yy(j+1) = ty
    end do
    allocate(d(n+2))
    d(1) = xlim(1)
    if (n > 0) d(2:n+1) = xx
    d(n+2) = xlim(2)
    call sort_real(d)
    m = 1
    maxr = (d(2)-d(1)) * (ylim(2)-ylim(1))
    do i = 2, n+1
      if (d(i+1)-d(i) > d(m+1)-d(m)) m = i
    end do
    rect = [d(m), ylim(1), d(m+1), ylim(2)]
    maxr = (rect(3)-rect(1)) * (rect(4)-rect(2))
    do i = 1, n
      tl = xlim(1)
      tr = xlim(2)
      if (i < n) then
        do j = i+1, n
          if (xx(j) > tl .and. xx(j) < tr) then
            trial = (tr-tl) * (yy(j)-yy(i))
            if (trial > maxr .and. tr-tl > width .and. yy(j)-yy(i) > height) then
              maxr = trial
              rect = [tl, yy(i), tr, yy(j)]
            end if
            if (xx(j) > xx(i)) then
              tr = xx(j)
            else
              tl = xx(j)
            end if
          end if
        end do
      end if
      trial = (tr-tl) * (ylim(2)-yy(i))
      if (trial > maxr .and. tr-tl > width .and. ylim(2)-yy(i) > height) then
        maxr = trial
        rect = [tl, yy(i), tr, ylim(2)]
      end if
      ri = xlim(2)
      li = xlim(1)
      do j = 1, n
        if (yy(j) < yy(i) .and. xx(j) > xx(i)) ri = min(ri, xx(j))
        if (yy(j) < yy(i) .and. xx(j) < xx(i)) li = max(li, xx(j))
      end do
      trial = (ri-li) * (yy(i)-ylim(1))
      if (trial > maxr .and. ri-li > width .and. yy(i)-ylim(1) > height) then
        maxr = trial
        rect = [li, ylim(1), ri, yy(i)]
      end if
    end do
    if (rect(3)-rect(1) < width .or. rect(4)-rect(2) < height) then
      status = 2
      center = nan_dp()
      rect = nan_dp()
      area = nan_dp()
      return
    end if
    center = [0.5_dp*(rect(1)+rect(3)), 0.5_dp*(rect(2)+rect(4))]
    area = (rect(3)-rect(1)) * (rect(4)-rect(2))
  end subroutine largest_empty_rexhaustive


  pure real(dp) function norm_pdf(x, mean, sd) result(density)
    real(dp), intent(in) :: x !! Evaluation point for the normal density.
    real(dp), intent(in) :: mean !! Mean of the normal distribution.
    real(dp), intent(in) :: sd !! Positive standard deviation of the normal distribution.
    real(dp), parameter :: sqrt_two_pi = 2.5066282746310005024_dp
    if (sd <= 0.0_dp .or. ieee_is_nan(x) .or. ieee_is_nan(mean) .or. ieee_is_nan(sd)) then
      density = nan_dp()
    else
      density = exp(-0.5_dp*((x-mean)/sd)**2) / (sqrt_two_pi*sd)
    end if
  end function norm_pdf


  pure subroutine gbayes_update(mean_prior, var_prior, stat, var_stat, mean_post, var_post, status, &
                                var_future, mean_pred, var_pred)
    real(dp), intent(in) :: mean_prior !! Prior mean for the scalar Gaussian parameter.
    real(dp), intent(in) :: var_prior !! Positive prior variance for the scalar Gaussian parameter.
    real(dp), intent(in) :: stat !! Observed approximately Gaussian statistic.
    real(dp), intent(in) :: var_stat !! Positive sampling variance of the observed statistic.
    real(dp), intent(out) :: mean_post !! Posterior mean after the Gaussian conjugate update.
    real(dp), intent(out) :: var_post !! Posterior variance after the Gaussian conjugate update.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid variances.
    real(dp), optional, intent(in) :: var_future !! Positive future sampling variance for predictive output.
    real(dp), optional, intent(out) :: mean_pred !! Predictive mean; equals posterior mean when requested.
    real(dp), optional, intent(out) :: var_pred !! Predictive variance including future sampling variance.
    status = 0
    mean_post = nan_dp()
    var_post = nan_dp()
    if (present(mean_pred)) mean_pred = nan_dp()
    if (present(var_pred)) var_pred = nan_dp()
    if (var_prior <= 0.0_dp .or. var_stat <= 0.0_dp) then
      status = 1
      return
    end if
    if (ieee_is_nan(mean_prior) .or. ieee_is_nan(var_prior) .or. ieee_is_nan(stat) .or. &
        ieee_is_nan(var_stat)) then
      status = 2
      return
    end if
    var_post = 1.0_dp / (1.0_dp/var_prior + 1.0_dp/var_stat)
    mean_post = (mean_prior/var_prior + stat/var_stat) * var_post
    if (present(var_future)) then
      if (var_future <= 0.0_dp .or. ieee_is_nan(var_future)) then
        status = 3
        return
      end if
      if (present(mean_pred)) mean_pred = mean_post
      if (present(var_pred)) var_pred = var_post + var_future
    end if
  end subroutine gbayes_update


  pure real(dp) function gbayes_mix_pred_density(delta, v, mix, d0, v0, d1, v1) result(ans)
    real(dp), intent(in) :: delta !! Parameter value at which to evaluate the future posterior density.
    real(dp), intent(in) :: v !! Positive variance of the future statistic.
    real(dp), intent(in) :: mix !! Prior mixture weight for component zero, in the closed interval [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp) :: pv0, pv1
    if (v <= 0.0_dp .or. v0 <= 0.0_dp .or. v1 <= 0.0_dp .or. mix < 0.0_dp .or. mix > 1.0_dp) then
      ans = nan_dp()
      return
    end if
    pv0 = 1.0_dp / (1.0_dp/v0 + 1.0_dp/v)
    pv1 = 1.0_dp / (1.0_dp/v1 + 1.0_dp/v)
    ans = mix*norm_pdf(delta, d0, sqrt(pv0)) + (1.0_dp-mix)*norm_pdf(delta, d1, sqrt(pv1))
  end function gbayes_mix_pred_density


  pure real(dp) function gbayes_mix_pred_cdf(delta, v, mix, d0, v0, d1, v1) result(ans)
    real(dp), intent(in) :: delta !! Parameter value at which to evaluate the future posterior CDF.
    real(dp), intent(in) :: v !! Positive variance of the future statistic.
    real(dp), intent(in) :: mix !! Prior mixture weight for component zero, in the closed interval [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp) :: pv0, pv1
    if (v <= 0.0_dp .or. v0 <= 0.0_dp .or. v1 <= 0.0_dp .or. mix < 0.0_dp .or. mix > 1.0_dp) then
      ans = nan_dp()
      return
    end if
    pv0 = 1.0_dp / (1.0_dp/v0 + 1.0_dp/v)
    pv1 = 1.0_dp / (1.0_dp/v1 + 1.0_dp/v)
    ans = mix*norm_cdf((delta-d0)/sqrt(pv0)) + (1.0_dp-mix)*norm_cdf((delta-d1)/sqrt(pv1))
  end function gbayes_mix_pred_cdf


  pure subroutine gbayes_mix_post_components(x, v, mix, d0, v0, d1, v1, mixp, pm0, pv0, pm1, pv1, status)
    real(dp), intent(in) :: x !! Observed approximately Gaussian statistic.
    real(dp), intent(in) :: v !! Positive sampling variance of the observed statistic.
    real(dp), intent(in) :: mix !! Prior weight for component zero, in the closed interval [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp), intent(out) :: mixp !! Posterior mixture weight for component zero.
    real(dp), intent(out) :: pm0 !! Posterior mean for component zero.
    real(dp), intent(out) :: pv0 !! Posterior variance for component zero.
    real(dp), intent(out) :: pm1 !! Posterior mean for component one.
    real(dp), intent(out) :: pv1 !! Posterior variance for component one.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid variance or mixture inputs.
    real(dp) :: prior_odds, likelihood_ratio, post_odds
    status = 0
    if (v <= 0.0_dp .or. v0 <= 0.0_dp .or. v1 <= 0.0_dp .or. mix < 0.0_dp .or. mix > 1.0_dp) then
      status = 1
      mixp = nan_dp()
      pm0 = nan_dp()
      pm1 = nan_dp()
      pv0 = nan_dp()
      pv1 = nan_dp()
      return
    end if
    pv0 = 1.0_dp / (1.0_dp/v0 + 1.0_dp/v)
    pv1 = 1.0_dp / (1.0_dp/v1 + 1.0_dp/v)
    pm0 = (d0/v0 + x/v) * pv0
    pm1 = (d1/v1 + x/v) * pv1
    if (mix <= 0.0_dp) then
      mixp = 0.0_dp
    else if (mix >= 1.0_dp) then
      mixp = 1.0_dp
    else
      prior_odds = mix / (1.0_dp-mix)
      likelihood_ratio = norm_pdf(x, d0, sqrt(v+v0)) / norm_pdf(x, d1, sqrt(v+v1))
      post_odds = prior_odds * likelihood_ratio
      mixp = post_odds / (1.0_dp+post_odds)
    end if
  end subroutine gbayes_mix_post_components


  pure real(dp) function gbayes_mix_post_density(delta, x, v, mix, d0, v0, d1, v1) result(ans)
    real(dp), intent(in) :: delta !! Parameter value at which to evaluate the posterior density.
    real(dp), intent(in) :: x !! Observed approximately Gaussian statistic.
    real(dp), intent(in) :: v !! Positive sampling variance of the observed statistic.
    real(dp), intent(in) :: mix !! Prior mixture weight for component zero, in [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp) :: mixp, pm0, pv0, pm1, pv1
    integer :: status
    call gbayes_mix_post_components(x, v, mix, d0, v0, d1, v1, mixp, pm0, pv0, pm1, pv1, status)
    if (status /= 0) then
      ans = nan_dp()
    else
      ans = mixp*norm_pdf(delta, pm0, sqrt(pv0)) + (1.0_dp-mixp)*norm_pdf(delta, pm1, sqrt(pv1))
    end if
  end function gbayes_mix_post_density


  pure real(dp) function gbayes_mix_post_cdf(delta, x, v, mix, d0, v0, d1, v1) result(ans)
    real(dp), intent(in) :: delta !! Parameter value at which to evaluate the posterior CDF.
    real(dp), intent(in) :: x !! Observed approximately Gaussian statistic.
    real(dp), intent(in) :: v !! Positive sampling variance of the observed statistic.
    real(dp), intent(in) :: mix !! Prior mixture weight for component zero, in [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp) :: mixp, pm0, pv0, pm1, pv1
    integer :: status
    call gbayes_mix_post_components(x, v, mix, d0, v0, d1, v1, mixp, pm0, pv0, pm1, pv1, status)
    if (status /= 0) then
      ans = nan_dp()
    else
      ans = mixp*norm_cdf((delta-pm0)/sqrt(pv0)) + (1.0_dp-mixp)*norm_cdf((delta-pm1)/sqrt(pv1))
    end if
  end function gbayes_mix_post_cdf


  pure real(dp) function gbayes_mix_post_mean(x, v, mix, d0, v0, d1, v1) result(ans)
    real(dp), intent(in) :: x !! Observed approximately Gaussian statistic.
    real(dp), intent(in) :: v !! Positive sampling variance of the observed statistic.
    real(dp), intent(in) :: mix !! Prior mixture weight for component zero, in [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp) :: mixp, pm0, pv0, pm1, pv1
    integer :: status
    call gbayes_mix_post_components(x, v, mix, d0, v0, d1, v1, mixp, pm0, pv0, pm1, pv1, status)
    if (status /= 0) then
      ans = nan_dp()
    else
      ans = mixp*pm0 + (1.0_dp-mixp)*pm1
    end if
  end function gbayes_mix_post_mean


  pure real(dp) function gbayes1_power_np(d0, v0, delta, v, delta_w, alpha) result(power)
    real(dp), intent(in) :: d0 !! Mean of the single Gaussian prior.
    real(dp), intent(in) :: v0 !! Positive variance of the single Gaussian prior.
    real(dp), intent(in) :: delta !! True effect used to evaluate frequentist power.
    real(dp), intent(in) :: v !! Positive sampling variance of the future statistic.
    real(dp), intent(in) :: delta_w !! Effect threshold defining the posterior tail event.
    real(dp), intent(in) :: alpha !! Two-sided posterior probability threshold, strictly between zero and one.
    real(dp) :: pv, z, critical
    if (v0 <= 0.0_dp .or. v <= 0.0_dp .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      power = nan_dp()
      return
    end if
    pv = 1.0_dp / (1.0_dp/v0 + 1.0_dp/v)
    z = inv_norm_cdf(alpha/2.0_dp)
    critical = v * ((delta_w-sqrt(pv)*z)/pv - d0/v0)
    power = 1.0_dp - norm_cdf((critical-delta)/sqrt(v))
  end function gbayes1_power_np


  pure subroutine gbayes_mix_power_np(delta, v, delta_w, mix, d0, v0, d1, v1, interval, alpha, &
                                      critical, power, status)
    real(dp), intent(in) :: delta !! True effect used to evaluate frequentist power.
    real(dp), intent(in) :: v !! Positive sampling variance of the future statistic.
    real(dp), intent(in) :: delta_w !! Effect threshold defining the posterior lower-tail event.
    real(dp), intent(in) :: mix !! Prior mixture weight for component zero, in [0,1].
    real(dp), intent(in) :: d0 !! Mean of Gaussian prior component zero.
    real(dp), intent(in) :: v0 !! Positive variance of Gaussian prior component zero.
    real(dp), intent(in) :: d1 !! Mean of Gaussian prior component one.
    real(dp), intent(in) :: v1 !! Positive variance of Gaussian prior component one.
    real(dp), intent(in) :: interval(2) !! Bracketing interval for the critical observed statistic.
    real(dp), intent(in) :: alpha !! Two-sided posterior probability threshold, strictly between zero and one.
    real(dp), intent(out) :: critical !! Statistic value where posterior CDF at delta_w equals alpha/2.
    real(dp), intent(out) :: power !! Probability that a N(delta,v) statistic exceeds the critical value.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or a missing root bracket.
    real(dp) :: lo, hi, mid, flo, fhi, fmid, target
    integer :: iter
    critical = nan_dp()
    power = nan_dp()
    status = 0
    if (v <= 0.0_dp .or. v0 <= 0.0_dp .or. v1 <= 0.0_dp .or. mix < 0.0_dp .or. mix > 1.0_dp .or. &
        alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. interval(1) >= interval(2)) then
      status = 1
      return
    end if
    lo = interval(1)
    hi = interval(2)
    target = alpha / 2.0_dp
    flo = gbayes_mix_post_cdf(delta_w, lo, v, mix, d0, v0, d1, v1) - target
    fhi = gbayes_mix_post_cdf(delta_w, hi, v, mix, d0, v0, d1, v1) - target
    if (ieee_is_nan(flo) .or. ieee_is_nan(fhi) .or. flo*fhi > 0.0_dp) then
      status = 2
      return
    end if
    do iter = 1, 200
      mid = 0.5_dp * (lo + hi)
      fmid = gbayes_mix_post_cdf(delta_w, mid, v, mix, d0, v0, d1, v1) - target
      if (abs(fmid) <= 1.0e-13_dp .or. abs(hi-lo) <= 1.0e-12_dp*(1.0_dp+abs(mid))) exit
      if (flo*fmid <= 0.0_dp) then
        hi = mid
        fhi = fmid
      else
        lo = mid
        flo = fmid
      end if
    end do
    critical = 0.5_dp * (lo + hi)
    power = 1.0_dp - norm_cdf((critical-delta)/sqrt(v))
  end subroutine gbayes_mix_power_np




  pure subroutine gbayes2_tabulated(delta, prior, sd, delta_w, alpha, probability, status)
    real(dp), intent(in) :: delta(:) !! Strictly increasing grid for the prior density; its endpoints define integration limits.
    real(dp), intent(in) :: prior(:) !! Nonnegative prior-density values corresponding one-to-one with delta.
    real(dp), intent(in) :: sd !! Positive standard deviation of the future approximately Gaussian estimator.
    real(dp), intent(in) :: delta_w !! Effect threshold; the numerator integrates posterior success probability above this value.
    real(dp), intent(in) :: alpha !! Two-sided significance level strictly between zero and one.
    real(dp), intent(out) :: probability !! Integrated Bayesian probability, normalized by prior mass represented on the grid.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid grids, probabilities, or represented prior mass.
    real(dp) :: denom, numer, z, xa, xb, pa, pb, fa, fb, x0, p0, f0
    integer :: i, n

    status = 0
    probability = nan_dp()
    n = size(delta)
    if (n < 2 .or. size(prior) /= n .or. sd <= 0.0_dp) then
      status = 1
      return
    end if
    if (alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      status = 2
      return
    end if
    if (any(ieee_is_nan(delta)) .or. any(ieee_is_nan(prior)) .or. any(prior < 0.0_dp)) then
      status = 3
      return
    end if
    if (any(delta(2:n) <= delta(1:n-1))) then
      status = 4
      return
    end if

    denom = 0.0_dp
    do i = 1, n - 1
      denom = denom + 0.5_dp * (delta(i+1)-delta(i)) * (prior(i)+prior(i+1))
    end do
    if (denom <= 0.0_dp) then
      status = 5
      return
    end if

    z = inv_norm_cdf(1.0_dp-alpha/2.0_dp)
    numer = 0.0_dp
    do i = 1, n - 1
      xa = delta(i)
      xb = delta(i+1)
      if (xb <= delta_w) cycle
      pa = prior(i)
      pb = prior(i+1)
      fa = pa * (1.0_dp-norm_cdf((delta_w-xa)/sd+z))
      fb = pb * (1.0_dp-norm_cdf((delta_w-xb)/sd+z))
      if (xa >= delta_w) then
        numer = numer + 0.5_dp * (xb-xa) * (fa+fb)
      else
        x0 = delta_w
        p0 = pa + (pb-pa) * (x0-xa) / (xb-xa)
        f0 = p0 * (1.0_dp-norm_cdf(z))
        numer = numer + 0.5_dp * (xb-x0) * (f0+fb)
      end if
    end do
    probability = numer / denom
  end subroutine gbayes2_tabulated


  pure subroutine props_po_counts(counts, odds_ratio, ref_group, observed, predicted, status)
    integer, intent(in) :: counts(:,:) !! Nonnegative category counts with rows as groups and columns as ordered outcome categories.
    real(dp), intent(in) :: odds_ratio(:) !! Positive group-wise proportional-odds ratios versus the reference profile.
    integer, intent(in) :: ref_group !! One-based row of counts supplying the reference outcome distribution.
    real(dp), allocatable, intent(out) :: observed(:,:) !! Row-wise observed category proportions; NaN for empty groups.
    real(dp), allocatable, intent(out) :: predicted(:,:) !! Category proportions implied by the reference row and odds ratios.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, counts, reference row, or odds ratios.
    real(dp), allocatable :: pref(:), ptmp(:)
    integer :: i, ng, nc, total

    status = 0
    ng = size(counts,1)
    nc = size(counts,2)
    allocate(observed(ng,nc), predicted(ng,nc))
    observed = nan_dp()
    predicted = nan_dp()
    if (ng < 1 .or. nc < 2 .or. size(odds_ratio) /= ng) then
      status = 1
      return
    end if
    if (ref_group < 1 .or. ref_group > ng .or. any(counts < 0)) then
      status = 2
      return
    end if
    if (any(odds_ratio <= 0.0_dp) .or. any(ieee_is_nan(odds_ratio))) then
      status = 3
      return
    end if
    do i = 1, ng
      total = sum(counts(i,:))
      if (total > 0) observed(i,:) = real(counts(i,:),dp) / real(total,dp)
    end do
    total = sum(counts(ref_group,:))
    if (total <= 0) then
      status = 4
      return
    end if
    allocate(pref(nc), ptmp(nc))
    pref = real(counts(ref_group,:),dp) / real(total,dp)
    do i = 1, ng
      call pomodm(pref, odds_ratio(i), ptmp)
      if (any(ieee_is_nan(ptmp))) then
        status = 5
        return
      end if
      predicted(i,:) = ptmp
    end do
  end subroutine props_po_counts


  pure subroutine props_trans_counts(states, nstates, counts, proportions, status)
    integer, intent(in) :: states(:,:) !! Subject-by-time state codes 1:nstates; nonpositive entries are missing.
    integer, intent(in) :: nstates !! Number of ordered or unordered states represented by positive state codes.
    integer, allocatable, intent(out) :: counts(:,:,:) !! Counts indexed by previous state, current state, and time transition.
    real(dp), allocatable, intent(out) :: proportions(:,:,:) !! Row-conditional transition proportions; NaN for empty rows.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid state count or positive state code above nstates.
    integer :: i, t, a, b, den, nt

    status = 0
    nt = size(states,2)
    if (nstates < 1 .or. nt < 2) then
      allocate(counts(0,0,0), proportions(0,0,0))
      status = 1
      return
    end if
    if (any(states > nstates)) then
      allocate(counts(0,0,0), proportions(0,0,0))
      status = 2
      return
    end if
    allocate(counts(nstates,nstates,nt-1), proportions(nstates,nstates,nt-1))
    counts = 0
    proportions = nan_dp()
    do t = 1, nt - 1
      do i = 1, size(states,1)
        a = states(i,t)
        b = states(i,t+1)
        if (a <= 0 .or. b <= 0) cycle
        counts(a,b,t) = counts(a,b,t) + 1
      end do
      do a = 1, nstates
        den = sum(counts(a,:,t))
        if (den > 0) proportions(a,:,t) = real(counts(a,:,t),dp) / real(den,dp)
      end do
    end do
  end subroutine props_trans_counts


  pure subroutine ynbind_numeric(x, sort_columns, y, column_order, proportions, status)
    real(dp), intent(in) :: x(:,:) !! Numeric binary variables coded 0 or 1; NaNs represent missing observations.
    logical, intent(in) :: sort_columns !! Sort output columns in ascending order of nonmissing one proportions when true.
    real(dp), allocatable, intent(out) :: y(:,:) !! Bound binary matrix with NaNs preserved and optional column reordering.
    integer, allocatable, intent(out) :: column_order(:) !! One-based source-column index for each output column.
    real(dp), allocatable, intent(out) :: proportions(:) !! Nonmissing proportion of ones for each output column, in output order.
    integer, intent(out) :: status !! Zero on success; nonzero when a nonmissing numeric value is not exactly zero or one.
    real(dp), allocatable :: p(:)
    integer :: i, j, k, nvalid, tmpi
    real(dp) :: tmpr

    status = 0
    allocate(y(size(x,1),size(x,2)), column_order(size(x,2)), proportions(size(x,2)), p(size(x,2)))
    y = x
    do j = 1, size(x,2)
      column_order(j) = j
      nvalid = 0
      p(j) = nan_dp()
      do i = 1, size(x,1)
        if (ieee_is_nan(x(i,j))) cycle
        if (x(i,j) /= 0.0_dp .and. x(i,j) /= 1.0_dp) then
          status = 1
          return
        end if
        nvalid = nvalid + 1
      end do
      if (nvalid > 0) p(j) = sum(merge(x(:,j),0.0_dp,.not.ieee_is_nan(x(:,j)))) / real(nvalid,dp)
    end do
    if (sort_columns) then
      do i = 2, size(x,2)
        j = i
        do while (j > 1)
          if (.not. ieee_is_nan(p(j-1)) .and. ieee_is_nan(p(j))) exit
          if (.not. ieee_is_nan(p(j-1)) .and. .not. ieee_is_nan(p(j))) then
            if (p(j-1) <= p(j)) exit
          end if
          tmpr = p(j-1)
          p(j-1) = p(j)
          p(j) = tmpr
          tmpi = column_order(j-1)
          column_order(j-1) = column_order(j)
          column_order(j) = tmpi
          j = j - 1
        end do
      end do
      y = x(:,column_order)
    end if
    do k = 1, size(x,2)
      proportions(k) = p(k)
    end do
  end subroutine ynbind_numeric


  pure subroutine reformm_missing_order(missing, order, missing_per_variable, missing_per_observation, &
                                        recommended_imputations, status)
    logical, intent(in) :: missing(:,:) !! Observation-by-variable missingness mask; true marks a missing value.
    integer, allocatable, intent(out) :: order(:) !! One-based variable indices sorted by decreasing missing count.
    integer, allocatable, intent(out) :: missing_per_variable(:) !! Number of missing observations for each variable.
    integer, allocatable, intent(out) :: missing_per_observation(:) !! Number of missing variables for each observation.
    integer, intent(out) :: recommended_imputations !! Upstream heuristic max(5, ceiling(percent of rows with any missingness)).
    integer, intent(out) :: status !! Zero on success; nonzero when there are no observations or variables.
    integer :: i, j, key, n_any, tmp

    status = 0
    if (size(missing,1) < 1 .or. size(missing,2) < 1) then
      allocate(order(0), missing_per_variable(0), missing_per_observation(0))
      recommended_imputations = 0
      status = 1
      return
    end if
    allocate(order(size(missing,2)), missing_per_variable(size(missing,2)))
    allocate(missing_per_observation(size(missing,1)))
    do j = 1, size(missing,2)
      missing_per_variable(j) = count(missing(:,j))
      order(j) = j
    end do
    do i = 1, size(missing,1)
      missing_per_observation(i) = count(missing(i,:))
    end do
    n_any = count(missing_per_observation > 0)
    recommended_imputations = max(5, ceiling(100.0_dp * real(n_any,dp) / real(size(missing,1),dp)))

    do i = 2, size(order)
      key = order(i)
      j = i - 1
      do while (j >= 1)
        if (missing_per_variable(order(j)) >= missing_per_variable(key)) exit
        order(j+1) = order(j)
        j = j - 1
      end do
      order(j+1) = key
    end do

    ! Preserve source order among variables having identical missing counts.
    do i = 2, size(order)
      j = i
      do while (j > 1)
        if (missing_per_variable(order(j-1)) /= missing_per_variable(order(j))) exit
        if (order(j-1) < order(j)) exit
        tmp = order(j-1)
        order(j-1) = order(j)
        order(j) = tmp
        j = j - 1
      end do
    end do
  end subroutine reformm_missing_order


  pure subroutine mapply_group_mean(x, group, means, group_levels, counts, status)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable numeric matrix to summarize within integer-coded groups.
    integer, intent(in) :: group(:) !! Group code per observation; positive codes are included and others are omitted.
    real(dp), allocatable, intent(out) :: means(:,:) !! Group-by-variable arithmetic means; NaNs are omitted per variable.
    integer, allocatable, intent(out) :: group_levels(:) !! Sorted positive group codes represented in the output rows.
    integer, allocatable, intent(out) :: counts(:,:) !! Nonmissing observation counts by group and variable.
    integer, intent(out) :: status !! Zero on success; nonzero when dimensions disagree or no positive groups exist.
    integer :: i, j, g, ng, row
    integer, allocatable :: tmp_levels(:)
    real(dp) :: s

    status = 0
    if (size(x,1) /= size(group)) then
      allocate(means(0,0), group_levels(0), counts(0,0))
      status = 1
      return
    end if
    allocate(tmp_levels(size(group)))
    ng = 0
    do i = 1, size(group)
      if (group(i) <= 0) cycle
      if (ng == 0 .or. .not. any(tmp_levels(1:ng) == group(i))) then
        ng = ng + 1
        tmp_levels(ng) = group(i)
      end if
    end do
    if (ng == 0) then
      allocate(means(0,size(x,2)), group_levels(0), counts(0,size(x,2)))
      status = 2
      return
    end if
    call sort_integer_prefix(tmp_levels, ng)
    allocate(group_levels(ng), means(ng,size(x,2)), counts(ng,size(x,2)))
    group_levels = tmp_levels(1:ng)
    means = nan_dp()
    counts = 0
    do row = 1, ng
      g = group_levels(row)
      do j = 1, size(x,2)
        s = 0.0_dp
        do i = 1, size(x,1)
          if (group(i) /= g) cycle
          if (ieee_is_nan(x(i,j))) cycle
          s = s + x(i,j)
          counts(row,j) = counts(row,j) + 1
        end do
        if (counts(row,j) > 0) means(row,j) = s / real(counts(row,j),dp)
      end do
    end do
  end subroutine mapply_group_mean


  pure subroutine soprob_markov_ord(intercepts, initial_offset, transition_offset, initial_state, absorbing, &
                                    probabilities, status)
    real(dp), intent(in) :: intercepts(:) !! Ordered-logit intercepts; length is one fewer than the number of states.
    real(dp), intent(in) :: initial_offset(:) !! Non-intercept linear-predictor offsets for the first-time distribution by cutpoint.
    real(dp), intent(in) :: transition_offset(:,:,:) !! Offsets by later time, previous state, and cutpoint; shape (nt-1,k,k-1).
    integer, intent(in) :: initial_state !! One-based initial conditioning state; used only to reject an absorbing baseline state.
    logical, intent(in) :: absorbing(:) !! True for absorbing states; length equals the number of states.
    real(dp), allocatable, intent(out) :: probabilities(:,:) !! Time-by-state exact occupancy probabilities.
    integer, intent(out) :: status !! Zero on success; nonzero for inconsistent dimensions or invalid probabilities/state.
    integer :: k, nt, i, j, t
    real(dp), allocatable :: cp(:,:), cell(:), cdf(:)
    real(dp) :: z

    status = 0
    k = size(intercepts) + 1
    if (k < 2 .or. size(initial_offset) /= k-1 .or. size(absorbing) /= k) then
      allocate(probabilities(0,0))
      status = 1
      return
    end if
    if (initial_state < 1 .or. initial_state > k) then
      allocate(probabilities(0,0))
      status = 2
      return
    end if
    if (absorbing(initial_state)) then
      allocate(probabilities(0,0))
      status = 3
      return
    end if
    if (size(transition_offset,2) /= k .or. size(transition_offset,3) /= k-1) then
      allocate(probabilities(0,0))
      status = 4
      return
    end if
    nt = size(transition_offset,1) + 1
    allocate(probabilities(nt,k), cp(k,k), cell(k), cdf(k-1))

    do j = 1, k-1
      z = intercepts(j) + initial_offset(j)
      if (z >= 0.0_dp) then
        cdf(j) = 1.0_dp / (1.0_dp + exp(-z))
      else
        cdf(j) = exp(z) / (1.0_dp + exp(z))
      end if
    end do
    cell(1) = 1.0_dp - cdf(1)
    do j = 2, k-1
      cell(j) = cdf(j-1) - cdf(j)
    end do
    cell(k) = cdf(k-1)
    if (any(cell < -100.0_dp * epsilon(1.0_dp))) then
      probabilities = nan_dp()
      status = 5
      return
    end if
    cell = max(cell, 0.0_dp)
    if (sum(cell) <= 0.0_dp) then
      probabilities = nan_dp()
      status = 5
      return
    end if
    probabilities(1,:) = cell / sum(cell)

    do t = 2, nt
      cp = 0.0_dp
      do i = 1, k
        if (absorbing(i)) then
          cp(i,i) = 1.0_dp
          cycle
        end if
        do j = 1, k-1
          z = intercepts(j) + transition_offset(t-1,i,j)
          if (z >= 0.0_dp) then
            cdf(j) = 1.0_dp / (1.0_dp + exp(-z))
          else
            cdf(j) = exp(z) / (1.0_dp + exp(z))
          end if
        end do
        cell(1) = 1.0_dp - cdf(1)
        do j = 2, k-1
          cell(j) = cdf(j-1) - cdf(j)
        end do
        cell(k) = cdf(k-1)
        if (any(cell < -100.0_dp * epsilon(1.0_dp))) then
          probabilities = nan_dp()
          status = 6
          return
        end if
        cell = max(cell,0.0_dp)
        if (sum(cell) <= 0.0_dp) then
          probabilities = nan_dp()
          status = 6
          return
        end if
        cp(i,:) = cell / sum(cell)
      end do
      probabilities(t,:) = matmul(transpose(cp), probabilities(t-1,:))
    end do
  end subroutine soprob_markov_ord


  pure subroutine ord_group_boot_mean(y, m_candidates, aprob, grouped, m_used, status, bootstrap_indices)
    real(dp), intent(in) :: y(:) !! Ordinal numeric response; NaNs are ignored when forming groups and preserved in output.
    integer, intent(in) :: m_candidates(:) !! Candidate minimum group sizes, tested in the supplied order; all must be positive.
    real(dp), intent(in) :: aprob !! Minimum approximate complete-bootstrap-coverage probability when bootstrap indices are absent.
    real(dp), allocatable, intent(out) :: grouped(:) !! Group-mean representation aligned with y; NaNs are preserved.
    integer, intent(out) :: m_used !! First candidate minimum group size satisfying the requested coverage criterion.
    integer, intent(out) :: status !! Zero on success; 1 invalid inputs, 2 invalid bootstrap indices, 3 no candidate succeeds.
    integer, intent(in), optional :: &
      bootstrap_indices(:,:) !! One-based original-row bootstrap indices; each column is one resample.
    integer, allocatable :: ord(:), grp_sorted(:), grp_orig(:), counts(:)
    real(dp), allocatable :: vals(:), means(:)
    integer :: i, j, k, n, nn, ng, mm, b, idx, keyi
    real(dp) :: keyv, prob
    logical :: ok

    n = size(y)
    allocate(grouped(n))
    grouped = nan_dp()
    m_used = 0
    status = 0
    if (n == 0 .or. size(m_candidates) == 0 .or. any(m_candidates <= 0) .or. &
        aprob < 0.0_dp .or. aprob > 1.0_dp) then
      status = 1
      return
    end if
    nn = count(.not. ieee_is_nan(y))
    if (nn == 0) then
      status = 1
      return
    end if
    allocate(ord(nn), vals(nn), grp_sorted(nn), grp_orig(n))
    grp_orig = 0
    k = 0
    do i = 1, n
      if (.not. ieee_is_nan(y(i))) then
        k = k + 1
        ord(k) = i
        vals(k) = y(i)
      end if
    end do
    do i = 2, nn
      keyv = vals(i)
      keyi = ord(i)
      j = i - 1
      do while (j >= 1)
        if (vals(j) <= keyv) exit
        vals(j+1) = vals(j)
        ord(j+1) = ord(j)
        j = j - 1
      end do
      vals(j+1) = keyv
      ord(j+1) = keyi
    end do

    do k = 1, size(m_candidates)
      mm = m_candidates(k)
      call cut_gn(vals, mm, grp_sorted)
      ng = maxval(grp_sorted)
      if (ng <= 0) cycle
      allocate(counts(ng), means(ng))
      counts = 0
      means = 0.0_dp
      grp_orig = 0
      do i = 1, nn
        j = grp_sorted(i)
        counts(j) = counts(j) + 1
        means(j) = means(j) + vals(i)
        grp_orig(ord(i)) = j
      end do
      do j = 1, ng
        means(j) = means(j) / real(counts(j), dp)
      end do
      do i = 1, n
        if (grp_orig(i) > 0) grouped(i) = means(grp_orig(i))
      end do

      ok = .true.
      if (present(bootstrap_indices)) then
        do b = 1, size(bootstrap_indices, 2)
          do j = 1, ng
            ok = .false.
            do i = 1, size(bootstrap_indices, 1)
              idx = bootstrap_indices(i,b)
              if (idx < 1 .or. idx > n) then
                status = 2
                deallocate(counts, means)
                return
              end if
              if (grp_orig(idx) == j) then
                ok = .true.
                exit
              end if
            end do
            if (.not. ok) exit
          end do
          if (.not. ok) exit
        end do
      else
        prob = 1.0_dp - sum(exp(-real(counts, dp)))
        ok = prob > aprob
      end if
      if (ok) then
        m_used = mm
        deallocate(counts, means)
        return
      end if
      deallocate(counts, means)
    end do
    grouped = nan_dp()
    status = 3
  end subroutine ord_group_boot_mean


  pure subroutine clowess_smooth(x, y, fraction, iterations, x_smooth, y_smooth, status)
    real(dp), intent(in) :: x(:) !! Predictor coordinates; rows containing NaN x or y are omitted.
    real(dp), intent(in) :: y(:) !! Response values corresponding one-to-one with x.
    real(dp), intent(in) :: fraction !! Fraction of complete observations used in each local linear fit, in (0,1].
    integer, intent(in) :: iterations !! Number of robust bisquare reweighting iterations; must be nonnegative.
    real(dp), allocatable, intent(out) :: x_smooth(:) !! Complete predictor values sorted in ascending order.
    real(dp), allocatable, intent(out) :: y_smooth(:) !! LOWESS fitted values corresponding to x_smooth.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, parameters, or no complete rows.
    real(dp), allocatable :: xs(:), ys(:), robust(:), fit(:), dist(:), absres(:)
    real(dp) :: keyx, keyy, h, dx, wt, s0, s1, s2, t0, t1, den, med, u
    real(dp) :: ymin, ymax
    integer :: i, j, nc, k, r, it

    status = 0
    if (size(x) /= size(y) .or. fraction <= 0.0_dp .or. fraction > 1.0_dp .or. iterations < 0) then
      allocate(x_smooth(0), y_smooth(0))
      status = 1
      return
    end if
    nc = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (nc == 0) then
      allocate(x_smooth(0), y_smooth(0))
      status = 2
      return
    end if
    allocate(xs(nc), ys(nc), robust(nc), fit(nc), dist(nc), absres(nc))
    k = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        k = k + 1
        xs(k) = x(i)
        ys(k) = y(i)
      end if
    end do
    do i = 2, nc
      keyx = xs(i)
      keyy = ys(i)
      j = i - 1
      do while (j >= 1)
        if (xs(j) <= keyx) exit
        xs(j+1) = xs(j)
        ys(j+1) = ys(j)
        j = j - 1
      end do
      xs(j+1) = keyx
      ys(j+1) = keyy
    end do
    r = max(2, min(nc, ceiling(fraction * real(nc, dp))))
    robust = 1.0_dp
    do it = 0, iterations
      do i = 1, nc
        do j = 1, nc
          dist(j) = abs(xs(j) - xs(i))
        end do
        call sort_real(dist)
        h = dist(r)
        s0 = 0.0_dp
        s1 = 0.0_dp
        s2 = 0.0_dp
        t0 = 0.0_dp
        t1 = 0.0_dp
        do j = 1, nc
          dx = xs(j) - xs(i)
          if (h == 0.0_dp) then
            wt = merge(robust(j), 0.0_dp, dx == 0.0_dp)
          else if (abs(dx) < h) then
            u = abs(dx) / h
            wt = robust(j) * (1.0_dp - u**3)**3
          else if (abs(dx) == h) then
            wt = 0.0_dp
          else
            wt = 0.0_dp
          end if
          s0 = s0 + wt
          s1 = s1 + wt * dx
          s2 = s2 + wt * dx * dx
          t0 = t0 + wt * ys(j)
          t1 = t1 + wt * ys(j) * dx
        end do
        if (s0 <= 0.0_dp) then
          fit(i) = ys(i)
        else
          den = s0*s2 - s1*s1
          if (abs(den) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, s0*s2)) then
            fit(i) = t0 / s0
          else
            fit(i) = (t0*s2 - t1*s1) / den
          end if
        end if
      end do
      if (it == iterations) exit
      absres = abs(ys - fit)
      med = median_sorted_copy(absres)
      if (med <= 10.0_dp*epsilon(1.0_dp)) exit
      do j = 1, nc
        u = absres(j) / (6.0_dp*med)
        if (u < 1.0_dp) then
          robust(j) = (1.0_dp - u*u)**2
        else
          robust(j) = 0.0_dp
        end if
      end do
    end do

    if (iterations > 0) then
      ymin = minval(ys)
      ymax = maxval(ys)
      if (any(fit < ymin) .or. any(fit > ymax)) then
        robust = 1.0_dp
        do i = 1, nc
          do j = 1, nc
            dist(j) = abs(xs(j) - xs(i))
          end do
          call sort_real(dist)
          h = dist(r)
          s0 = 0.0_dp
          s1 = 0.0_dp
          s2 = 0.0_dp
          t0 = 0.0_dp
          t1 = 0.0_dp
          do j = 1, nc
            dx = xs(j) - xs(i)
            if (h == 0.0_dp) then
              wt = merge(1.0_dp, 0.0_dp, dx == 0.0_dp)
            else if (abs(dx) < h) then
              u = abs(dx) / h
              wt = (1.0_dp - u**3)**3
            else
              wt = 0.0_dp
            end if
            s0 = s0 + wt
            s1 = s1 + wt*dx
            s2 = s2 + wt*dx*dx
            t0 = t0 + wt*ys(j)
            t1 = t1 + wt*ys(j)*dx
          end do
          den = s0*s2 - s1*s1
          if (s0 <= 0.0_dp) then
            fit(i) = ys(i)
          else if (abs(den) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, s0*s2)) then
            fit(i) = t0 / s0
          else
            fit(i) = (t0*s2 - t1*s1) / den
          end if
        end do
      end if
    end if
    allocate(x_smooth(nc), y_smooth(nc))
    x_smooth = xs
    y_smooth = fit
  end subroutine clowess_smooth


  pure subroutine curve_smooth_observed(x, y, id, fraction, iterations, x_out, y_out, id_out, status)
    real(dp), intent(in) :: x(:) !! Curve predictor coordinates; incomplete x/y rows are omitted from outputs.
    real(dp), intent(in) :: y(:) !! Curve responses corresponding one-to-one with x.
    integer, intent(in) :: id(:) !! Integer curve identifiers corresponding one-to-one with x and y.
    real(dp), intent(in) :: fraction !! LOWESS neighborhood fraction used for curves having at least three unique x values.
    integer, intent(in) :: iterations !! Number of robust LOWESS iterations.
    real(dp), allocatable, intent(out) :: x_out(:) !! Complete predictor coordinates in original row order.
    real(dp), allocatable, intent(out) :: y_out(:) !! Smoothed responses in original row order for complete rows.
    integer, allocatable, intent(out) :: id_out(:) !! Curve identifiers aligned with x_out and y_out.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or an internal LOWESS failure.
    integer, allocatable :: subpos(:)
    real(dp), allocatable :: sx(:), sy(:), xs(:), ys(:)
    integer :: n, nc, i, j, k, m, uid, st
    logical :: seen

    status = 0
    n = size(x)
    if (size(y) /= n .or. size(id) /= n .or. fraction <= 0.0_dp .or. fraction > 1.0_dp .or. iterations < 0) then
      allocate(x_out(0), y_out(0), id_out(0))
      status = 1
      return
    end if
    nc = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    allocate(x_out(nc), y_out(nc), id_out(nc))
    k = 0
    do i = 1, n
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        k = k + 1
        x_out(k) = x(i)
        y_out(k) = y(i)
        id_out(k) = id(i)
      end if
    end do
    do i = 1, nc
      uid = id_out(i)
      seen = .false.
      do j = 1, i-1
        if (id_out(j) == uid) then
          seen = .true.
          exit
        end if
      end do
      if (seen) cycle
      m = count(id_out == uid)
      allocate(subpos(m), sx(m), sy(m))
      k = 0
      do j = 1, nc
        if (id_out(j) == uid) then
          k = k + 1
          subpos(k) = j
          sx(k) = x_out(j)
          sy(k) = y_out(j)
        end if
      end do
      if (count_unique_real(sx) >= 3) then
        call clowess_smooth(sx, sy, fraction, iterations, xs, ys, st)
        if (st /= 0) then
          status = 2
          deallocate(subpos, sx, sy)
          return
        end if
        do j = 1, m
          y_out(subpos(j)) = linear_interp_sorted(xs, ys, sx(j))
        end do
        deallocate(xs, ys)
      end if
      deallocate(subpos, sx, sy)
    end do
  end subroutine curve_smooth_observed


  pure integer function count_unique_real(x) result(nunique)
    real(dp), intent(in) :: x(:) !! Complete real values whose exact distinct count is requested.
    real(dp), allocatable :: z(:)
    integer :: i
    if (size(x) == 0) then
      nunique = 0
      return
    end if
    z = x
    call sort_real(z)
    nunique = 1
    do i = 2, size(z)
      if (z(i) /= z(i-1)) nunique = nunique + 1
    end do
  end function count_unique_real


  pure real(dp) function linear_interp_sorted(x, y, xout) result(value)
    real(dp), intent(in) :: x(:) !! Ascending interpolation abscissae.
    real(dp), intent(in) :: y(:) !! Interpolation ordinates aligned with x.
    real(dp), intent(in) :: xout !! Point within or at the bounds of x at which to interpolate.
    integer :: i, n
    n = min(size(x), size(y))
    if (n == 0) then
      value = nan_dp()
      return
    end if
    if (xout <= x(1)) then
      value = y(1)
      return
    end if
    if (xout >= x(n)) then
      value = y(n)
      return
    end if
    do i = 1, n-1
      if (xout <= x(i+1)) then
        if (x(i+1) == x(i)) then
          value = y(i+1)
        else
          value = y(i) + (y(i+1)-y(i)) * (xout-x(i)) / (x(i+1)-x(i))
        end if
        return
      end if
    end do
    value = y(n)
  end function linear_interp_sorted


  pure subroutine sort_integer_prefix(x, n)
    integer, intent(inout) :: x(:) !! Integer work array whose first n entries are sorted in ascending order.
    integer, intent(in) :: n !! Number of leading entries to sort; must be between zero and size(x).
    integer :: i, j, key

    do i = 2, n
      key = x(i)
      j = i - 1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j - 1
      end do
      x(j+1) = key
    end do
  end subroutine sort_integer_prefix



  pure subroutine wtd_loess_noiter(x, y, weights, span, degree, fitted, status)
    real(dp), intent(in) :: x(:) !! Predictor values; NaN rows are omitted from each local fit and return NaN fitted values.
    real(dp), intent(in) :: y(:) !! Response values aligned with x; NaN rows are omitted from each local fit.
    real(dp), intent(in) :: weights(:) !! Nonnegative observation weights aligned with x and y; zero weights do not contribute.
    real(dp), intent(in) :: span !! Fraction of complete observations used in each local neighborhood; valid range (0,1].
    integer, intent(in) :: degree !! Local polynomial degree; supported values are 0, 1, or 2.
    real(dp), allocatable, intent(out) :: fitted(:) !! Fitted values aligned with the original observations; invalid rows are NaN.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, span/degree, or singular local fits.
    integer :: n, nc, i, j, q, kk, a, b, st
    integer, allocatable :: good(:), order(:)
    real(dp), allocatable :: dist(:), z(:,:), rhs(:), beta(:), ww(:)
    real(dp) :: h, u, tw, xi

    status = 0
    n = size(x)
    allocate(fitted(n))
    fitted = nan_dp()
    if (size(y) /= n .or. size(weights) /= n .or. span <= 0.0_dp .or. span > 1.0_dp .or. &
        degree < 0 .or. degree > 2 .or. any(weights < 0.0_dp)) then
      status = 1
      return
    end if
    nc = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y) .and. .not. ieee_is_nan(weights) .and. weights > 0.0_dp)
    if (nc < degree + 1) then
      status = 2
      return
    end if
    allocate(good(nc))
    j = 0
    do i = 1, n
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i)) .and. .not. ieee_is_nan(weights(i)) .and. &
          weights(i) > 0.0_dp) then
        j = j + 1
        good(j) = i
      end if
    end do
    q = max(degree + 1, ceiling(span * real(nc, dp)))
    q = min(q, nc)
    allocate(dist(nc), order(nc), ww(nc), z(nc,degree+1), rhs(degree+1), beta(degree+1))

    do i = 1, n
      if (ieee_is_nan(x(i))) cycle
      xi = x(i)
      do j = 1, nc
        dist(j) = abs(x(good(j)) - xi)
        order(j) = j
      end do
      do a = 1, nc-1
        kk = a
        do b = a+1, nc
          if (dist(order(b)) < dist(order(kk))) kk = b
        end do
        if (kk /= a) then
          j = order(a)
          order(a) = order(kk)
          order(kk) = j
        end if
      end do
      h = dist(order(q))
      if (h <= 0.0_dp) then
        h = maxval(dist)
        if (h <= 0.0_dp) then
          tw = sum(weights(good))
          if (tw > 0.0_dp) fitted(i) = sum(weights(good) * y(good)) / tw
          cycle
        end if
      end if
      ww = 0.0_dp
      z = 0.0_dp
      do j = 1, nc
        u = abs(x(good(j)) - xi) / h
        if (u <= 1.0_dp) ww(j) = weights(good(j)) * (1.0_dp - u**3)**3
        z(j,1) = 1.0_dp
        if (degree >= 1) z(j,2) = x(good(j)) - xi
        if (degree >= 2) z(j,3) = (x(good(j)) - xi)**2
      end do
      rhs = matmul(transpose(z), ww * y(good))
      call solvet_vector(matmul(transpose(z), spread(ww,2,degree+1) * z), rhs, beta, status=st)
      if (st == 0) then
        fitted(i) = beta(1)
      else
        status = max(status, 3)
      end if
    end do
  end subroutine wtd_loess_noiter


  pure subroutine sim_markov_ord(intercepts, transition_offset, initial_state, absorbing, uniforms, carry, states, status)
    real(dp), intent(in) :: intercepts(:) !! Ordered-logit intercepts; length is one fewer than the number of ordinal states.
    real(dp), intent(in) :: transition_offset(:,:,:,:) !! Offsets by subject, time, previous state, and cutpoint.
    integer, intent(in) :: initial_state(:) !! One-based baseline state for each subject; length n and not absorbing.
    logical, intent(in) :: absorbing(:) !! True for absorbing states; length k.
    real(dp), intent(in) :: uniforms(:,:) !! Caller-supplied U(0,1) draws by subject and time; shape (n,nt).
    logical, intent(in) :: carry !! If true, carry an absorbing state forward; otherwise later states are zero after absorption.
    integer, allocatable, intent(out) :: states(:,:) !! Simulated states; zero marks post-absorption truncation.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, states, uniforms, or probabilities.
    integer :: n, nt, k, id, t, j, prev, chosen
    real(dp), allocatable :: exceed(:), cell(:)
    real(dp) :: z, total, cumu, u

    status = 0
    n = size(initial_state)
    k = size(intercepts) + 1
    if (k < 2 .or. size(absorbing) /= k .or. size(transition_offset,1) /= n .or. &
        size(transition_offset,3) /= k .or. size(transition_offset,4) /= k-1 .or. &
        size(uniforms,1) /= n .or. size(uniforms,2) /= size(transition_offset,2)) then
      allocate(states(0,0))
      status = 1
      return
    end if
    nt = size(uniforms,2)
    if (any(initial_state < 1) .or. any(initial_state > k)) then
      allocate(states(0,0))
      status = 2
      return
    end if
    do id = 1, n
      if (absorbing(initial_state(id))) then
        allocate(states(0,0))
        status = 3
        return
      end if
    end do
    if (any(uniforms < 0.0_dp) .or. any(uniforms >= 1.0_dp) .or. any(ieee_is_nan(uniforms))) then
      allocate(states(0,0))
      status = 4
      return
    end if
    allocate(states(n,nt), exceed(k-1), cell(k))
    states = 0
    do id = 1, n
      prev = initial_state(id)
      do t = 1, nt
        if (t > 1 .and. absorbing(prev)) then
          if (carry) then
            states(id,t) = prev
            cycle
          else
            exit
          end if
        end if
        do j = 1, k-1
          z = intercepts(j) + transition_offset(id,t,prev,j)
          if (z >= 0.0_dp) then
            exceed(j) = 1.0_dp / (1.0_dp + exp(-z))
          else
            exceed(j) = exp(z) / (1.0_dp + exp(z))
          end if
        end do
        cell(1) = 1.0_dp - exceed(1)
        do j = 2, k-1
          cell(j) = exceed(j-1) - exceed(j)
        end do
        cell(k) = exceed(k-1)
        if (any(cell < -100.0_dp * epsilon(1.0_dp))) then
          states = 0
          status = 5
          return
        end if
        cell = max(cell, 0.0_dp)
        total = sum(cell)
        if (total <= 0.0_dp) then
          states = 0
          status = 5
          return
        end if
        cell = cell / total
        u = uniforms(id,t)
        cumu = 0.0_dp
        chosen = k
        do j = 1, k
          cumu = cumu + cell(j)
          if (u < cumu) then
            chosen = j
            exit
          end if
        end do
        states(id,t) = chosen
        prev = chosen
      end do
    end do
  end subroutine sim_markov_ord


  pure subroutine impute_median(x, y, imputed, status)
    real(dp), intent(in) :: x(:) !! Numeric vector; NaNs are treated as missing values to be imputed.
    real(dp), allocatable, intent(out) :: y(:) !! Copy of x with each NaN replaced by the median of observed values.
    integer, allocatable, intent(out) :: imputed(:) !! One-based positions that were missing and therefore imputed.
    integer, intent(out) :: status !! Zero on success; 1 when all values are missing and no median can be computed.
    real(dp), allocatable :: observed(:)
    real(dp) :: fill
    integer :: i, nobs, nmiss, jo, jm

    nobs = count(.not. ieee_is_nan(x))
    nmiss = size(x) - nobs
    allocate(y(size(x)), imputed(nmiss))
    y = x
    status = 0
    if (nmiss == 0) return
    if (nobs == 0) then
      status = 1
      return
    end if
    allocate(observed(nobs))
    jo = 0
    jm = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i))) then
        jm = jm + 1
        imputed(jm) = i
      else
        jo = jo + 1
        observed(jo) = x(i)
      end if
    end do
    fill = median_sorted_copy(observed)
    do i = 1, nmiss
      y(imputed(i)) = fill
    end do
  end subroutine impute_median


  pure subroutine impute_constant(x, fill, y, imputed, status)
    real(dp), intent(in) :: x(:) !! Numeric vector; NaNs are treated as missing values to be imputed.
    real(dp), intent(in) :: fill !! Constant replacement value used for every missing element; may itself be NaN.
    real(dp), allocatable, intent(out) :: y(:) !! Copy of x with each NaN replaced by fill.
    integer, allocatable, intent(out) :: imputed(:) !! One-based positions that were missing and therefore replaced.
    integer, intent(out) :: status !! Zero on success; retained for API consistency and future validation.
    integer :: i, jm, nmiss

    nmiss = count(ieee_is_nan(x))
    allocate(y(size(x)), imputed(nmiss))
    y = x
    jm = 0
    do i = 1, size(x)
      if (ieee_is_nan(x(i))) then
        jm = jm + 1
        imputed(jm) = i
        y(i) = fill
      end if
    end do
    status = 0
  end subroutine impute_constant


  pure subroutine mov_stats_n_raw(x, y, eps, nignore, xinc, x_center, mean_y, q1, median_y, q3, n_window, status)
    real(dp), intent(in) :: x(:) !! Continuous ordering variable; rows with NaN x or y are omitted before windowing.
    real(dp), intent(in) :: y(:) !! Numeric response whose moving mean and quartiles are computed.
    integer, intent(in) :: eps !! Half-window size in observations; each target uses indices target-eps through target+eps.
    integer, intent(in) :: nignore !! One-based first target index and symmetric tail exclusion, as in space='n'.
    integer, intent(in) :: xinc !! Positive increment in sorted-observation indices between successive target windows.
    real(dp), allocatable, intent(out) :: x_center(:) !! Mean x value within each moving window, matching upstream .xmean.
    real(dp), allocatable, intent(out) :: mean_y(:) !! Arithmetic mean response within each moving window.
    real(dp), allocatable, intent(out) :: q1(:) !! Type-7 first quartile of response within each moving window.
    real(dp), allocatable, intent(out) :: median_y(:) !! Type-7 median of response within each moving window.
    real(dp), allocatable, intent(out) :: q3(:) !! Type-7 third quartile of response within each moving window.
    integer, allocatable, intent(out) :: n_window(:) !! Number of complete observations contributing to each window.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions, controls, or too few complete rows.
    real(dp), allocatable :: xx(:), yy(:), ys(:)
    integer, allocatable :: ord(:)
    integer :: i, j, n, ntarget, target, lo, hi, m, nc

    status = 0
    if (size(x) /= size(y) .or. eps < 0 .or. nignore < 1 .or. xinc < 1) then
      allocate(x_center(0), mean_y(0), q1(0), median_y(0), q3(0), n_window(0))
      status = 1
      return
    end if
    nc = count(.not. ieee_is_nan(x) .and. .not. ieee_is_nan(y))
    if (nc < 1 .or. nc - nignore + 1 < nignore) then
      allocate(x_center(0), mean_y(0), q1(0), median_y(0), q3(0), n_window(0))
      status = 2
      return
    end if
    allocate(xx(nc), yy(nc), ord(nc))
    j = 0
    do i = 1, size(x)
      if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) then
        j = j + 1
        xx(j) = x(i)
        yy(j) = y(i)
      end if
    end do
    call sort_index_real(xx, ord)
    xx = xx(ord)
    yy = yy(ord)
    n = nc
    ntarget = (n - nignore + 1 - nignore) / xinc + 1
    allocate(x_center(ntarget), mean_y(ntarget), q1(ntarget), median_y(ntarget), q3(ntarget), n_window(ntarget))
    do i = 1, ntarget
      target = nignore + (i - 1) * xinc
      lo = max(1, target - eps)
      hi = min(n, target + eps)
      m = hi - lo + 1
      n_window(i) = m
      x_center(i) = sum(xx(lo:hi)) / real(m, dp)
      mean_y(i) = sum(yy(lo:hi)) / real(m, dp)
      allocate(ys(m))
      ys = yy(lo:hi)
      call sort_real(ys)
      q1(i) = quantile_type7_sorted(ys, 0.25_dp)
      median_y(i) = quantile_type7_sorted(ys, 0.5_dp)
      q3(i) = quantile_type7_sorted(ys, 0.75_dp)
      deallocate(ys)
    end do
  end subroutine mov_stats_n_raw


  pure subroutine summarize_mean_codes(x, by, levels, means, counts, status)
    real(dp), intent(in) :: x(:,:) !! Numeric observations by variables; NaNs are omitted separately for each variable.
    integer, intent(in) :: by(:,:) !! Integer grouping codes by observation and grouping variable; shape (n, nby).
    integer, allocatable, intent(out) :: levels(:,:) !! Unique grouping-code combinations in lexicographic ascending order.
    real(dp), allocatable, intent(out) :: means(:,:) !! Group-by-variable arithmetic means; NaN when a cell has no data.
    integer, allocatable, intent(out) :: counts(:,:) !! Nonmissing counts corresponding to each grouped mean.
    integer, intent(out) :: status !! Zero on success; nonzero for row mismatch or an empty grouping specification.
    integer, allocatable :: order(:), raw_levels(:,:)
    logical :: is_new
    integer :: n, p, nb, ng, i, j, g, k, pos, tmpi
    real(dp) :: tmpx

    n = size(x,1)
    p = size(x,2)
    nb = size(by,2)
    if (size(by,1) /= n .or. nb < 1) then
      allocate(levels(0,0), means(0,0), counts(0,0))
      status = 1
      return
    end if
    allocate(raw_levels(max(1,n),nb))
    ng = 0
    do i = 1, n
      is_new = .true.
      do g = 1, ng
        if (all(raw_levels(g,:) == by(i,:))) then
          is_new = .false.
          exit
        end if
      end do
      if (is_new) then
        ng = ng + 1
        raw_levels(ng,:) = by(i,:)
      end if
    end do
    allocate(order(ng))
    order = [(i, i=1,ng)]
    do i = 2, ng
      tmpi = order(i)
      pos = i - 1
      do while (pos >= 1)
        is_new = .false.
        do k = 1, nb
          if (raw_levels(tmpi,k) < raw_levels(order(pos),k)) then
            is_new = .true.
            exit
          else if (raw_levels(tmpi,k) > raw_levels(order(pos),k)) then
            exit
          end if
        end do
        if (.not. is_new) exit
        order(pos+1) = order(pos)
        pos = pos - 1
      end do
      order(pos+1) = tmpi
    end do
    allocate(levels(ng,nb), means(ng,p), counts(ng,p))
    means = nan_dp()
    counts = 0
    do g = 1, ng
      levels(g,:) = raw_levels(order(g),:)
      do j = 1, p
        tmpx = 0.0_dp
        do i = 1, n
          if (all(by(i,:) == levels(g,:)) .and. .not. ieee_is_nan(x(i,j))) then
            tmpx = tmpx + x(i,j)
            counts(g,j) = counts(g,j) + 1
          end if
        end do
        if (counts(g,j) > 0) means(g,j) = tmpx / real(counts(g,j), dp)
      end do
    end do
    status = 0
  end subroutine summarize_mean_codes



  pure subroutine areg_tran_numeric(z, type_code, knots, nlevels, transformed, status)
    real(dp), intent(in) :: z(:) !! Numeric predictor values; categorical values use integer codes starting at one.
    character(len=1), intent(in) :: type_code !! Transformation code: l linear, s restricted cubic spline, or c categorical.
    real(dp), intent(in) :: knots(:) !! Explicit restricted-cubic-spline knots for type s; ignored otherwise.
    integer, intent(in) :: nlevels !! Number of categorical levels for type c; ignored otherwise.
    real(dp), allocatable, intent(out) :: transformed(:,:) !! Transformed design matrix matching Hmisc aregTran numerical columns.
    integer, intent(out) :: status !! Zero on success; nonzero for unsupported type, invalid knots/codes, or dimensions.
    real(dp), allocatable :: b(:,:)
    integer :: i, code

    status = 0
    select case (type_code)
    case ('l', 'L')
      allocate(transformed(size(z),1))
      transformed(:,1) = z
    case ('s', 'S')
      call rcspline_eval(z, knots, .true., .false., 2, b)
      if (size(b,2) == 0) then
        allocate(transformed(size(z),0))
        status = 1
      else
        call move_alloc(b, transformed)
      end if
    case ('c', 'C')
      if (nlevels < 2) then
        allocate(transformed(size(z),0))
        status = 2
        return
      end if
      allocate(transformed(size(z),nlevels-1))
      transformed = 0.0_dp
      do i = 1, size(z)
        if (ieee_is_nan(z(i))) then
          transformed(i,:) = nan_dp()
          cycle
        end if
        code = nint(z(i))
        if (abs(z(i)-real(code,dp)) > 1.0e-10_dp .or. code < 1 .or. code > nlevels) then
          transformed(i,:) = nan_dp()
          status = 3
        else if (code > 1) then
          transformed(i,code-1) = 1.0_dp
        end if
      end do
    case default
      allocate(transformed(size(z),0))
      status = 4
    end select
  end subroutine areg_tran_numeric


  pure subroutine areg_linear_fit(x, y, coefficients, fitted_values, residuals, &
                                  rsquared, mean_abs_error, median_abs_error, status)
    real(dp), intent(in) :: x(:,:) !! Numeric predictor matrix for the all-linear xtype path of Hmisc areg.
    real(dp), intent(in) :: y(:) !! Numeric response for the linear ytype path; incomplete rows are omitted for fitting.
    real(dp), allocatable, intent(out) :: coefficients(:) !! Intercept followed by linear predictor coefficients.
    real(dp), allocatable, intent(out) :: fitted_values(:) !! Fitted values aligned to original rows; incomplete rows are NaN.
    real(dp), allocatable, intent(out) :: residuals(:) !! Residuals aligned to original rows; incomplete rows are NaN.
    real(dp), intent(out) :: rsquared !! Ordinary coefficient of determination on complete rows.
    real(dp), intent(out) :: mean_abs_error !! Mean absolute residual on complete rows.
    real(dp), intent(out) :: median_abs_error !! Median absolute residual on complete rows.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid shape, too few complete rows, or singular fit.
    real(dp), allocatable :: xc(:,:), yc(:), coef(:), res(:), fit(:), ae(:)
    integer, allocatable :: row(:)
    integer :: i, n, p, nc, st

    n = size(x,1)
    p = size(x,2)
    rsquared = nan_dp()
    mean_abs_error = nan_dp()
    median_abs_error = nan_dp()
    if (size(y) /= n) then
      status = 1
      allocate(coefficients(0), fitted_values(0), residuals(0))
      return
    end if
    allocate(row(n))
    nc = 0
    do i = 1, n
      if (.not. ieee_is_nan(y(i)) .and. all(.not. ieee_is_nan(x(i,:)))) then
        nc = nc + 1
        row(nc) = i
      end if
    end do
    allocate(fitted_values(n), residuals(n))
    fitted_values = nan_dp()
    residuals = nan_dp()
    if (nc <= p) then
      status = 2
      allocate(coefficients(0))
      return
    end if
    allocate(xc(nc,p), yc(nc))
    do i = 1, nc
      xc(i,:) = x(row(i),:)
      yc(i) = y(row(i))
    end do
    call lm_fit_qr_bare(xc, yc, .true., coef, res, fit, rsquared, st)
    if (st /= 0) then
      status = 3
      allocate(coefficients(0))
      return
    end if
    coefficients = coef
    do i = 1, nc
      fitted_values(row(i)) = fit(i)
      residuals(row(i)) = res(i)
    end do
    allocate(ae(nc))
    ae = abs(res)
    mean_abs_error = sum(ae) / real(nc,dp)
    median_abs_error = median_sorted_copy(ae)
    status = 0
  end subroutine areg_linear_fit


  pure subroutine dataframe_reduce_numeric(data, fracmiss, minprev, reduced, keep_indices, removed_reason, status)
    real(dp), intent(in) :: data(:,:) !! Numeric observation-by-variable matrix; NaNs denote missing values.
    real(dp), intent(in) :: fracmiss !! Maximum allowed missing fraction before a variable is removed, in [0,1].
    real(dp), intent(in) :: minprev !! Minimum prevalence for either value of a binary numeric variable, in [0,0.5].
    real(dp), allocatable, intent(out) :: reduced(:,:) !! Matrix containing retained variables in original order.
    integer, allocatable, intent(out) :: keep_indices(:) !! One-based original indices of retained variables.
    integer, allocatable, intent(out) :: removed_reason(:) !! Code per variable: 0 retained, 1 missing, 2 prevalence.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid thresholds.
    integer, allocatable :: keep_tmp(:)
    real(dp), allocatable :: vals(:)
    integer :: i, j, n, p, nk, nv, nnon, count1, count2
    real(dp) :: missfrac

    n = size(data,1)
    p = size(data,2)
    if (fracmiss < 0.0_dp .or. fracmiss > 1.0_dp .or. minprev < 0.0_dp .or. minprev > 0.5_dp) then
      status = 1
      allocate(reduced(n,0), keep_indices(0), removed_reason(p))
      removed_reason = 0
      return
    end if
    allocate(removed_reason(p), keep_tmp(p), vals(max(1,n)))
    removed_reason = 0
    nk = 0
    do j = 1, p
      nnon = count(.not. ieee_is_nan(data(:,j)))
      missfrac = merge(real(n-nnon,dp)/real(n,dp), 0.0_dp, n > 0)
      if (missfrac > fracmiss) then
        removed_reason(j) = 1
        cycle
      end if
      nv = 0
      do i = 1, n
        if (ieee_is_nan(data(i,j))) cycle
        if (nv == 0) then
          nv = 1
          vals(1) = data(i,j)
        else if (all(abs(vals(1:nv)-data(i,j)) > 1.0e-12_dp)) then
          nv = nv + 1
          if (nv <= size(vals)) vals(nv) = data(i,j)
          if (nv > 2) exit
        end if
      end do
      if (nv == 2 .and. minprev > 0.0_dp .and. nnon > 0) then
        count1 = count((.not. ieee_is_nan(data(:,j))) .and. abs(data(:,j)-vals(1)) <= 1.0e-12_dp)
        count2 = nnon - count1
        if (real(min(count1,count2),dp)/real(nnon,dp) < minprev) then
          removed_reason(j) = 2
          cycle
        end if
      end if
      nk = nk + 1
      keep_tmp(nk) = j
    end do
    allocate(keep_indices(nk), reduced(n,nk))
    if (nk > 0) then
      keep_indices = keep_tmp(:nk)
      do j = 1, nk
        reduced(:,j) = data(:,keep_indices(j))
      end do
    end if
    status = 0
  end subroutine dataframe_reduce_numeric


  pure subroutine data_rep_numeric(data, unique_rows, counts, original_rows, status)
    real(dp), intent(in) :: data(:,:) !! Numeric descriptor matrix; rows containing NaNs are omitted.
    real(dp), allocatable, intent(out) :: unique_rows(:,:) !! Distinct complete descriptor combinations in first-occurrence order.
    integer, allocatable, intent(out) :: counts(:) !! Frequencies of each distinct descriptor combination.
    integer, intent(out) :: original_rows !! Number of complete rows represented in the output.
    integer, intent(out) :: status !! Zero on success; nonzero only for an empty variable dimension.
    real(dp), allocatable :: u(:,:)
    integer, allocatable :: c(:)
    integer :: i, j, n, p, nu
    logical :: found

    n = size(data,1)
    p = size(data,2)
    original_rows = 0
    if (p < 1) then
      status = 1
      allocate(unique_rows(0,0), counts(0))
      return
    end if
    allocate(u(max(1,n),p), c(max(1,n)))
    c = 0
    nu = 0
    do i = 1, n
      if (any(ieee_is_nan(data(i,:)))) cycle
      original_rows = original_rows + 1
      found = .false.
      do j = 1, nu
        if (all(abs(u(j,:)-data(i,:)) <= 1.0e-12_dp)) then
          c(j) = c(j) + 1
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        nu = nu + 1
        u(nu,:) = data(i,:)
        c(nu) = 1
      end if
    end do
    allocate(unique_rows(nu,p), counts(nu))
    if (nu > 0) then
      unique_rows = u(:nu,:)
      counts = c(:nu)
    end if
    status = 0
  end subroutine data_rep_numeric


  pure subroutine varclus_numeric(data, similarity_code, transform_code, method_code, similarity, npair, &
                                  merge_left, merge_right, height, status)
    real(dp), intent(in) :: data(:,:) !! Numeric data matrix; NaNs are handled pairwise in similarity calculations.
    integer, intent(in) :: similarity_code !! Similarity: 1 Pearson, 2 Spearman, or 3 thirty times Hoeffding D.
    integer, intent(in) :: transform_code !! Correlation transform: 0 none, 1 square, or 2 absolute value; ignored for Hoeffding.
    integer, intent(in) :: method_code !! Agglomeration linkage: 1 complete, 2 single, or 3 average.
    real(dp), allocatable, intent(out) :: similarity(:,:) !! Variable similarity matrix used by clustering.
    integer, allocatable, intent(out) :: npair(:,:) !! Pairwise complete observation counts.
    integer, allocatable, intent(out) :: merge_left(:) !! Left cluster id for each merge; initial leaves are 1..p.
    integer, allocatable, intent(out) :: merge_right(:) !! Right cluster id for each merge; new clusters are p+step.
    real(dp), allocatable, intent(out) :: height(:) !! Dissimilarity 1-similarity at each agglomeration.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid options or undefined similarities.
    real(dp), allocatable :: aad(:,:), mad(:,:), dist(:,:)
    logical, allocatable :: members(:,:), active(:)
    real(dp) :: best, dlink, total
    integer :: p, a, b, i, u, v, step, ca, cb, idnew, nlink

    p = size(data,2)
    if (p < 1 .or. similarity_code < 1 .or. similarity_code > 3 .or. &
        transform_code < 0 .or. transform_code > 2 .or. method_code < 1 .or. method_code > 3) then
      status = 1
      allocate(similarity(0,0), npair(0,0), merge_left(0), merge_right(0), height(0))
      return
    end if
    allocate(similarity(p,p), npair(p,p))
    if (similarity_code == 3) then
      allocate(aad(p,p), mad(p,p))
      call hoeffd(data, similarity, aad, mad, npair)
      similarity = 30.0_dp * similarity
    else
      call rcorr(data, similarity_code == 2, similarity, npair)
      if (transform_code == 1) similarity = similarity**2
      if (transform_code == 2) similarity = abs(similarity)
    end if
    if (any(ieee_is_nan(similarity))) then
      status = 2
      allocate(merge_left(0), merge_right(0), height(0))
      return
    end if
    allocate(merge_left(max(0,p-1)), merge_right(max(0,p-1)), height(max(0,p-1)))
    if (p == 1) then
      status = 0
      return
    end if
    allocate(dist(p,p), members(p,2*p-1), active(2*p-1))
    dist = 1.0_dp - similarity
    members = .false.
    active = .false.
    do i = 1, p
      members(i,i) = .true.
      active(i) = .true.
    end do
    do step = 1, p-1
      best = huge(1.0_dp)
      ca = 0
      cb = 0
      do a = 1, p+step-1
        if (.not. active(a)) cycle
        do b = a+1, p+step-1
          if (.not. active(b)) cycle
          total = 0.0_dp
          nlink = 0
          if (method_code == 1) dlink = -huge(1.0_dp)
          if (method_code == 2) dlink = huge(1.0_dp)
          do u = 1, p
            if (.not. members(u,a)) cycle
            do v = 1, p
              if (.not. members(v,b)) cycle
              select case (method_code)
              case (1)
                dlink = max(dlink, dist(u,v))
              case (2)
                dlink = min(dlink, dist(u,v))
              case default
                total = total + dist(u,v)
                nlink = nlink + 1
              end select
            end do
          end do
          if (method_code == 3) dlink = total / real(nlink,dp)
          if (dlink < best) then
            best = dlink
            ca = a
            cb = b
          end if
        end do
      end do
      merge_left(step) = ca
      merge_right(step) = cb
      height(step) = best
      idnew = p + step
      members(:,idnew) = members(:,ca) .or. members(:,cb)
      active(ca) = .false.
      active(cb) = .false.
      active(idnew) = .true.
    end do
    status = 0
  end subroutine varclus_numeric


  pure subroutine redun_numeric(data, r2_threshold, adjusted, kept, removed, removed_r2, status)
    real(dp), intent(in) :: data(:,:) !! Complete numeric variable matrix for ordinary linear redundancy analysis.
    real(dp), intent(in) :: r2_threshold !! R-squared threshold at or above which the most redundant variable is removed.
    logical, intent(in) :: adjusted !! Apply the adjusted R-squared correction used by Hmisc redun type='adjusted'.
    integer, allocatable, intent(out) :: kept(:) !! One-based indices of variables remaining after iterative removal.
    integer, allocatable, intent(out) :: removed(:) !! One-based indices of removed variables in removal order.
    real(dp), allocatable, intent(out) :: removed_r2(:) !! R-squared value that triggered each removal.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid data or threshold.
    logical, allocatable :: active(:)
    integer, allocatable :: remtmp(:), keeptmp(:)
    real(dp), allocatable :: rtmp(:), rhs(:,:), coef(:), res(:), fit(:)
    real(dp) :: r2, ar2, best
    integer :: n, p, m, i, j, k, nr, nk, st, bestj

    n = size(data,1)
    p = size(data,2)
    if (n < 2 .or. p < 1 .or. any(ieee_is_nan(data)) .or. r2_threshold < 0.0_dp .or. r2_threshold > 1.0_dp) then
      status = 1
      allocate(kept(0), removed(0), removed_r2(0))
      return
    end if
    allocate(active(p), remtmp(p), rtmp(p), keeptmp(p))
    active = .true.
    nr = 0
    do
      m = count(active)
      if (m < 2) exit
      best = -1.0_dp
      bestj = 0
      do j = 1, p
        if (.not. active(j)) cycle
        allocate(rhs(n,m-1))
        k = 0
        do i = 1, p
          if (active(i) .and. i /= j) then
            k = k + 1
            rhs(:,k) = data(:,i)
          end if
        end do
        call lm_fit_qr_bare(rhs, data(:,j), .true., coef, res, fit, r2, st)
        deallocate(rhs)
        if (st == 0 .and. .not. ieee_is_nan(r2)) then
          ar2 = r2
          if (adjusted .and. n > m) ar2 = max(0.0_dp, 1.0_dp-(1.0_dp-r2)*real(n-1,dp)/real(n-m,dp))
          if (ar2 > best) then
            best = ar2
            bestj = j
          end if
        end if
        if (allocated(coef)) deallocate(coef)
        if (allocated(res)) deallocate(res)
        if (allocated(fit)) deallocate(fit)
      end do
      if (bestj == 0 .or. best < r2_threshold) exit
      nr = nr + 1
      remtmp(nr) = bestj
      rtmp(nr) = best
      active(bestj) = .false.
    end do
    nk = 0
    do j = 1, p
      if (active(j)) then
        nk = nk + 1
        keeptmp(nk) = j
      end if
    end do
    allocate(kept(nk), removed(nr), removed_r2(nr))
    if (nk > 0) kept = keeptmp(:nk)
    if (nr > 0) then
      removed = remtmp(:nr)
      removed_r2 = rtmp(:nr)
    end if
    status = 0
  end subroutine redun_numeric


  pure subroutine gbayes_seq(est, vest, direction, cutoff_low, cutoff_high, prior_mu, prior_sigma, &
                             probability, post_mean, post_sd, status)
    real(dp), intent(in) :: est(:) !! Sequential parameter estimates, one per look/simulation row.
    real(dp), intent(in) :: vest(:) !! Sampling variances corresponding one-to-one with est and strictly positive.
    integer, intent(in) :: direction(:) !! Assertion codes: -1 parameter below cutoff, 1 above cutoff, 0 inside interval.
    real(dp), intent(in) :: cutoff_low(:) !! Scalar cutoff for -1 or 1; lower interval endpoint for direction zero.
    real(dp), intent(in) :: cutoff_high(:) !! Upper interval endpoint for direction zero; ignored otherwise.
    real(dp), intent(in) :: prior_mu(:) !! Gaussian prior mean for each assertion.
    real(dp), intent(in) :: prior_sigma(:) !! Positive Gaussian prior standard deviation for each assertion.
    real(dp), allocatable, intent(out) :: probability(:,:) !! Posterior assertion probabilities, rows by assertions.
    real(dp), allocatable, intent(out) :: post_mean(:,:) !! Posterior Gaussian means, rows by assertions.
    real(dp), allocatable, intent(out) :: post_sd(:,:) !! Posterior Gaussian standard deviations, rows by assertions.
    integer, intent(out) :: status !! Zero on success; nonzero for inconsistent dimensions or invalid variances/assertions.
    real(dp) :: pv, vv, m, sd
    integer :: i, j, n, na

    n = size(est)
    na = size(direction)
    if (size(vest) /= n .or. size(cutoff_low) /= na .or. size(cutoff_high) /= na .or. &
        size(prior_mu) /= na .or. size(prior_sigma) /= na .or. any(vest <= 0.0_dp) .or. &
        any(prior_sigma <= 0.0_dp) .or. any(abs(direction) > 1)) then
      status = 1
      allocate(probability(0,0), post_mean(0,0), post_sd(0,0))
      return
    end if
    allocate(probability(n,na), post_mean(n,na), post_sd(n,na))
    do j = 1, na
      pv = prior_sigma(j)**2
      do i = 1, n
        if (ieee_is_nan(est(i)) .or. ieee_is_nan(vest(i))) then
          probability(i,j) = nan_dp()
          post_mean(i,j) = nan_dp()
          post_sd(i,j) = nan_dp()
          cycle
        end if
        vv = 1.0_dp / (1.0_dp/pv + 1.0_dp/vest(i))
        m = (prior_mu(j)/pv + est(i)/vest(i)) * vv
        sd = sqrt(vv)
        post_mean(i,j) = m
        post_sd(i,j) = sd
        select case (direction(j))
        case (-1)
          probability(i,j) = norm_cdf((cutoff_low(j)-m)/sd)
        case (1)
          probability(i,j) = 1.0_dp - norm_cdf((cutoff_low(j)-m)/sd)
        case default
          if (cutoff_high(j) < cutoff_low(j)) then
            probability(i,j) = nan_dp()
          else
            probability(i,j) = norm_cdf((cutoff_high(j)-m)/sd) - norm_cdf((cutoff_low(j)-m)/sd)
          end if
        end select
      end do
    end do
    status = 0
  end subroutine gbayes_seq


  pure subroutine soprob_markov_ordm_core(initial, transition, occupancy, status)
    real(dp), intent(in) :: initial(:) !! State probabilities at the first measurement time, summing to one.
    real(dp), intent(in) :: transition(:,:,:) !! Conditional transition matrices: previous state by current state by interval.
    real(dp), allocatable, intent(out) :: occupancy(:,:) !! Unconditional occupancy probabilities by time and state.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid dimensions or probability matrices.
    integer :: k, nt, t, i

    k = size(initial)
    nt = size(transition,3) + 1
    if (k < 1 .or. size(transition,1) /= k .or. size(transition,2) /= k .or. &
        any(initial < 0.0_dp) .or. abs(sum(initial)-1.0_dp) > 1.0e-8_dp .or. &
        any(transition < 0.0_dp)) then
      status = 1
      allocate(occupancy(0,0))
      return
    end if
    do t = 1, nt-1
      do i = 1, k
        if (abs(sum(transition(i,:,t))-1.0_dp) > 1.0e-8_dp) then
          status = 2
          allocate(occupancy(0,0))
          return
        end if
      end do
    end do
    allocate(occupancy(nt,k))
    occupancy(1,:) = initial
    do t = 2, nt
      occupancy(t,:) = matmul(transpose(transition(:,:,t-1)), occupancy(t-1,:))
    end do
    status = 0
  end subroutine soprob_markov_ordm_core


  pure subroutine areg_boot_linear(x, y, bootstrap_indices, coefficients, rsquared_app, &
                                   rsquared_validated, coef_boot, nfail, status)
    real(dp), intent(in) :: x(:,:) !! Numeric predictor matrix for the all-linear areg path, observations by predictors.
    real(dp), intent(in) :: y(:) !! Numeric response vector; this specialized path requires complete data.
    integer, intent(in) :: bootstrap_indices(:,:) !! One-based bootstrap row indices, with observations by bootstrap replicate.
    real(dp), allocatable, intent(out) :: coefficients(:) !! Original-sample intercept and linear coefficients.
    real(dp), intent(out) :: rsquared_app !! Apparent original-sample R-squared.
    real(dp), intent(out) :: rsquared_validated !! Bootstrap optimism-corrected R-squared.
    real(dp), allocatable, intent(out) :: coef_boot(:,:) !! Bootstrap coefficients by replicate; failed columns are NaN.
    integer, intent(out) :: nfail !! Number of bootstrap fits that were singular or otherwise unsuccessful.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid input or original fit failure.
    real(dp), allocatable :: res(:), fit(:), cb(:), rb(:), fb(:), xb(:,:), yb(:), pred(:)
    real(dp) :: r2b, r2test, sse, sst, ybar, optimism
    integer :: n, p, b, i, st, good

    n = size(x,1)
    p = size(x,2)
    rsquared_app = nan_dp()
    rsquared_validated = nan_dp()
    nfail = 0
    if (size(y) /= n .or. size(bootstrap_indices,1) /= n .or. n <= p + 1 .or. &
        any(ieee_is_nan(x)) .or. any(ieee_is_nan(y))) then
      status = 1
      allocate(coefficients(0), coef_boot(0,0))
      return
    end if
    if (any(bootstrap_indices < 1) .or. any(bootstrap_indices > n)) then
      status = 2
      allocate(coefficients(0), coef_boot(0,0))
      return
    end if
    call lm_fit_qr_bare(x, y, .true., coefficients, res, fit, rsquared_app, st)
    if (st /= 0) then
      status = 3
      if (.not. allocated(coefficients)) allocate(coefficients(0))
      allocate(coef_boot(0,0))
      return
    end if
    allocate(coef_boot(p+1,size(bootstrap_indices,2)))
    coef_boot = nan_dp()
    allocate(xb(n,p), yb(n), pred(n))
    ybar = sum(y) / real(n,dp)
    sst = sum((y-ybar)**2)
    optimism = 0.0_dp
    good = 0
    do b = 1, size(bootstrap_indices,2)
      do i = 1, n
        xb(i,:) = x(bootstrap_indices(i,b),:)
        yb(i) = y(bootstrap_indices(i,b))
      end do
      call lm_fit_qr_bare(xb, yb, .true., cb, rb, fb, r2b, st)
      if (st /= 0) then
        nfail = nfail + 1
        cycle
      end if
      coef_boot(:,b) = cb
      pred = cb(1)
      if (p > 0) pred = pred + matmul(x, cb(2:p+1))
      sse = sum((y-pred)**2)
      if (sst > tiny(1.0_dp)) then
        r2test = 1.0_dp - sse / sst
        optimism = optimism + r2b - r2test
        good = good + 1
      end if
    end do
    if (good > 0) rsquared_validated = rsquared_app - optimism / real(good,dp)
    status = 0
  end subroutine areg_boot_linear


  pure subroutine fit_mult_impute_combine(coefficients, variances, coef_mean, covariance, &
                                          variance_inflation, missing_info, dfmi, status)
    real(dp), intent(in) :: coefficients(:,:) !! Coefficient estimates, coefficient by imputation.
    real(dp), intent(in) :: variances(:,:,:) !! Within-imputation covariance matrices, coefficient by coefficient by imputation.
    real(dp), allocatable, intent(out) :: coef_mean(:) !! Rubin-combined mean coefficient estimates.
    real(dp), allocatable, intent(out) :: covariance(:,:) !! Rubin total covariance matrix.
    real(dp), allocatable, intent(out) :: variance_inflation(:) !! Total variance divided by mean within-imputation variance.
    real(dp), allocatable, intent(out) :: missing_info(:) !! Fraction of missing information for each coefficient.
    real(dp), allocatable, intent(out) :: dfmi(:) !! Approximate t degrees of freedom for each coefficient.
    integer, intent(out) :: status !! Zero on success; nonzero for inconsistent dimensions or fewer than two imputations.
    real(dp), allocatable :: ubar(:,:), bmat(:,:)
    real(dp) :: tau, u, bb
    integer :: p, m, i, j, k

    p = size(coefficients,1)
    m = size(coefficients,2)
    if (m < 2 .or. size(variances,1) /= p .or. size(variances,2) /= p .or. &
        size(variances,3) /= m .or. p < 1) then
      status = 1
      allocate(coef_mean(0), covariance(0,0), variance_inflation(0), missing_info(0), dfmi(0))
      return
    end if
    allocate(coef_mean(p), covariance(p,p), variance_inflation(p), missing_info(p), dfmi(p))
    allocate(ubar(p,p), bmat(p,p))
    coef_mean = sum(coefficients, dim=2) / real(m,dp)
    ubar = sum(variances, dim=3) / real(m,dp)
    bmat = 0.0_dp
    do k = 1, m
      do j = 1, p
        do i = 1, p
          bmat(i,j) = bmat(i,j) + (coefficients(i,k)-coef_mean(i)) * &
                                      (coefficients(j,k)-coef_mean(j))
        end do
      end do
    end do
    bmat = bmat / real(m-1,dp)
    covariance = ubar + (1.0_dp + 1.0_dp/real(m,dp)) * bmat
    do i = 1, p
      u = ubar(i,i)
      bb = bmat(i,i)
      if (u > 0.0_dp) then
        variance_inflation(i) = covariance(i,i) / u
        tau = (1.0_dp + 1.0_dp/real(m,dp)) * bb / u
        if (tau > 0.0_dp) then
          missing_info(i) = tau / (1.0_dp + tau)
          dfmi(i) = real(m-1,dp) * (1.0_dp + 1.0_dp/tau)**2
        else
          missing_info(i) = 0.0_dp
          dfmi(i) = huge(1.0_dp)
        end if
      else
        variance_inflation(i) = nan_dp()
        missing_info(i) = nan_dp()
        dfmi(i) = nan_dp()
      end if
    end do
    status = 0
  end subroutine fit_mult_impute_combine


  pure subroutine impute_transcan_numeric(variable, positions, imputed_values, imputation, result, status)
    real(dp), intent(in) :: variable(:) !! Original numeric variable, usually containing NaNs at imputed positions.
    integer, intent(in) :: positions(:) !! One-based positions corresponding to rows of imputed_values.
    real(dp), intent(in) :: imputed_values(:,:) !! Candidate imputations, missing observation by imputation replicate.
    integer, intent(in) :: imputation !! One-based imputation column to apply.
    real(dp), allocatable, intent(out) :: result(:) !! Copy of variable with selected imputed values inserted.
    integer, intent(out) :: status !! Zero on success; nonzero for inconsistent dimensions, positions, or imputation number.
    integer :: i

    if (size(imputed_values,1) /= size(positions) .or. imputation < 1 .or. &
        imputation > size(imputed_values,2) .or. any(positions < 1) .or. any(positions > size(variable))) then
      status = 1
      allocate(result(0))
      return
    end if
    result = variable
    do i = 1, size(positions)
      result(positions(i)) = imputed_values(i,imputation)
    end do
    status = 0
  end subroutine impute_transcan_numeric


  pure subroutine est_seq_sim_mean_difference(group, response, looks, estimates, variances, status)
    integer, intent(in) :: group(:,:,:) !! Treatment indicators 0/1, observation by simulation by parameter setting.
    real(dp), intent(in) :: response(:,:,:) !! Simulated responses aligned with group by observation, simulation, and setting.
    integer, intent(in) :: looks(:) !! Increasing cumulative observation counts at which estimates are requested.
    real(dp), allocatable, intent(out) :: estimates(:,:,:) !! Mean(group=1)-mean(group=0), look by simulation by parameter setting.
    real(dp), allocatable, intent(out) :: variances(:,:,:) !! Estimated variance of each mean difference in the same shape.
    integer, intent(out) :: status !! Zero on success; nonzero for dimension, look, group, or insufficient-group-size errors.
    integer :: n, ns, np, il, isim, ip, l, i, n0, n1
    real(dp) :: m0, m1, ss0, ss1, z

    n = size(group,1)
    ns = size(group,2)
    np = size(group,3)
    if (any(shape(response) /= shape(group)) .or. size(looks) < 1 .or. any(looks < 2) .or. &
        any(looks > n) .or. any(group < 0) .or. any(group > 1)) then
      status = 1
      allocate(estimates(0,0,0), variances(0,0,0))
      return
    end if
    allocate(estimates(size(looks),ns,np), variances(size(looks),ns,np))
    estimates = nan_dp()
    variances = nan_dp()
    do ip = 1, np
      do isim = 1, ns
        do il = 1, size(looks)
          l = looks(il)
          n0 = 0
          n1 = 0
          m0 = 0.0_dp
          m1 = 0.0_dp
          do i = 1, l
            if (ieee_is_nan(response(i,isim,ip))) cycle
            if (group(i,isim,ip) == 0) then
              n0 = n0 + 1
              m0 = m0 + response(i,isim,ip)
            else
              n1 = n1 + 1
              m1 = m1 + response(i,isim,ip)
            end if
          end do
          if (n0 < 2 .or. n1 < 2) cycle
          m0 = m0 / real(n0,dp)
          m1 = m1 / real(n1,dp)
          ss0 = 0.0_dp
          ss1 = 0.0_dp
          do i = 1, l
            z = response(i,isim,ip)
            if (ieee_is_nan(z)) cycle
            if (group(i,isim,ip) == 0) then
              ss0 = ss0 + (z-m0)**2
            else
              ss1 = ss1 + (z-m1)**2
            end if
          end do
          estimates(il,isim,ip) = m1 - m0
          variances(il,isim,ip) = ss0 / real(n0*(n0-1),dp) + ss1 / real(n1*(n1-1),dp)
        end do
      end do
    end do
    status = 0
  end subroutine est_seq_sim_mean_difference


  pure subroutine rm_boot_linear(time, y, eval_time, residual_indices, coefficients, &
                                 coef_boot, fitted_eval, status)
    real(dp), intent(in) :: time(:) !! Measurement times for the no-subject linear-time specialization of rm.boot.
    real(dp), intent(in) :: y(:) !! Response values aligned with time; this specialization requires complete data.
    real(dp), intent(in) :: eval_time(:) !! Times where the original and bootstrap fitted curves are evaluated.
    integer, intent(in) :: residual_indices(:,:) !! One-based residual resampling indices, observation by bootstrap replicate.
    real(dp), allocatable, intent(out) :: coefficients(:) !! Original intercept and linear time slope.
    real(dp), allocatable, intent(out) :: coef_boot(:,:) !! Bootstrap intercept/slope, coefficient by replicate.
    real(dp), allocatable, intent(out) :: fitted_eval(:,:) !! Fitted curves at eval_time: original in column 1, then bootstraps.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or a singular original/bootstrap fit.
    real(dp), allocatable :: x(:,:), res(:), fit(:), yb(:), cb(:), rb(:), fb(:)
    real(dp) :: r2
    integer :: n, b, i, st

    n = size(time)
    if (size(y) /= n .or. size(residual_indices,1) /= n .or. n < 3 .or. &
        any(ieee_is_nan(time)) .or. any(ieee_is_nan(y)) .or. any(residual_indices < 1) .or. &
        any(residual_indices > n)) then
      status = 1
      allocate(coefficients(0), coef_boot(0,0), fitted_eval(0,0))
      return
    end if
    allocate(x(n,1))
    x(:,1) = time
    call lm_fit_qr_bare(x, y, .true., coefficients, res, fit, r2, st)
    if (st /= 0) then
      status = 2
      if (.not. allocated(coefficients)) allocate(coefficients(0))
      allocate(coef_boot(0,0), fitted_eval(0,0))
      return
    end if
    allocate(coef_boot(2,size(residual_indices,2)))
    allocate(fitted_eval(size(eval_time),size(residual_indices,2)+1))
    fitted_eval(:,1) = coefficients(1) + coefficients(2) * eval_time
    allocate(yb(n))
    status = 0
    do b = 1, size(residual_indices,2)
      do i = 1, n
        yb(i) = fit(i) + res(residual_indices(i,b))
      end do
      call lm_fit_qr_bare(x, yb, .true., cb, rb, fb, r2, st)
      if (st /= 0) then
        status = 3
        coef_boot(:,b) = nan_dp()
        fitted_eval(:,b+1) = nan_dp()
      else
        coef_boot(:,b) = cb
        fitted_eval(:,b+1) = cb(1) + cb(2) * eval_time
      end if
    end do
    if (status /= 3) status = 0
  end subroutine rm_boot_linear


  pure real(dp) function ordinal_loglik_group(y, group, intercepts, beta) result(ll)
    integer, intent(in) :: y(:) !! One-based ordered outcome category for each observation.
    real(dp), intent(in) :: group(:) !! Numeric group/predictor value aligned with y.
    real(dp), intent(in) :: intercepts(:) !! Descending cumulative-logit intercepts for P(Y >= category j+1).
    real(dp), intent(in) :: beta !! Common proportional-odds slope multiplying group.
    integer :: i, j, k
    real(dp), allocatable :: exceed(:), cell(:)
    real(dp) :: z, pr

    ll = -huge(1.0_dp)
    if (size(group) /= size(y)) return
    k = size(intercepts) + 1
    if (k < 2 .or. any(y < 1) .or. any(y > k) .or. any(ieee_is_nan(group))) return
    if (size(intercepts) > 1) then
      if (any(intercepts(2:) > intercepts(:size(intercepts)-1))) return
    end if
    allocate(exceed(k-1), cell(k))
    ll = 0.0_dp
    do i = 1, size(y)
      do j = 1, k-1
        z = intercepts(j) + beta * group(i)
        if (z >= 0.0_dp) then
          exceed(j) = 1.0_dp / (1.0_dp + exp(-z))
        else
          exceed(j) = exp(z) / (1.0_dp + exp(z))
        end if
      end do
      cell(1) = 1.0_dp - exceed(1)
      do j = 2, k-1
        cell(j) = exceed(j-1) - exceed(j)
      end do
      cell(k) = exceed(k-1)
      pr = cell(y(i))
      if (pr <= tiny(1.0_dp)) then
        ll = -huge(1.0_dp)
        return
      end if
      ll = ll + log(pr)
    end do
  end function ordinal_loglik_group


  pure subroutine optimize_ordinal_group(y, group, intercepts, beta, loglik, status)
    integer, intent(in) :: y(:) !! One-based ordered outcome category.
    real(dp), intent(in) :: group(:) !! Numeric predictor aligned with y.
    real(dp), allocatable, intent(out) :: intercepts(:) !! Fitted descending cumulative-logit intercepts.
    real(dp), intent(out) :: beta !! Fitted common proportional-odds slope.
    real(dp), intent(out) :: loglik !! Maximized log likelihood.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid data or failed initialization.
    integer :: k, j, i, n, n_ge, iter
    real(dp), allocatable :: par(:), cand(:), step(:)
    real(dp) :: p, cur, trial
    logical :: improved

    n = size(y)
    k = maxval(y)
    status = 0
    if (n < 2 .or. size(group) /= n .or. minval(y) < 1 .or. k < 2 .or. any(ieee_is_nan(group))) then
      allocate(intercepts(0))
      beta = nan_dp()
      loglik = nan_dp()
      status = 1
      return
    end if
    allocate(intercepts(k-1), par(k), cand(k), step(k))
    do j = 1, k-1
      n_ge = count(y >= j+1)
      p = (real(n_ge,dp) + 0.5_dp) / (real(n,dp) + 1.0_dp)
      intercepts(j) = log(p / (1.0_dp-p))
    end do
    par(1:k-1) = intercepts
    par(k) = 0.0_dp
    step = 0.5_dp
    cur = ordinal_loglik_group(y, group, par(1:k-1), par(k))
    if (cur <= -0.5_dp*huge(1.0_dp)) then
      beta = nan_dp()
      loglik = nan_dp()
      status = 2
      return
    end if
    do iter = 1, 400
      improved = .false.
      do i = 1, k
        cand = par
        cand(i) = cand(i) + step(i)
        trial = ordinal_loglik_group(y, group, cand(1:k-1), cand(k))
        if (trial > cur) then
          par = cand
          cur = trial
          improved = .true.
        else
          cand = par
          cand(i) = cand(i) - step(i)
          trial = ordinal_loglik_group(y, group, cand(1:k-1), cand(k))
          if (trial > cur) then
            par = cand
            cur = trial
            improved = .true.
          end if
        end if
      end do
      if (.not. improved) step = 0.5_dp * step
      if (maxval(step) < 1.0e-7_dp) exit
    end do
    intercepts = par(1:k-1)
    beta = par(k)
    loglik = cur
  end subroutine optimize_ordinal_group


  pure subroutine ord_test_po_numeric(group, outcome, statistic, df, p_value, beta, status)
    real(dp), intent(in) :: group(:) !! Numeric group/predictor values used in the proportional-odds model.
    integer, intent(in) :: outcome(:) !! One-based ordered response categories with at least two observed levels.
    real(dp), intent(out) :: statistic !! Likelihood-ratio chi-square comparing group model with intercept-only model.
    integer, intent(out) :: df !! Model degrees of freedom; one for this numeric-group specialization.
    real(dp), intent(out) :: p_value !! Upper-tail chi-square probability for the likelihood-ratio statistic.
    real(dp), intent(out) :: beta !! Fitted common proportional-odds group slope.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid input or optimization failure.
    real(dp), allocatable :: ints_alt(:), ints_null(:), zero_group(:)
    real(dp) :: ll_alt, ll_null, beta_null
    integer :: st1, st0

    statistic = nan_dp()
    p_value = nan_dp()
    beta = nan_dp()
    df = 1
    if (size(group) /= size(outcome) .or. size(group) < 3) then
      status = 1
      return
    end if
    allocate(zero_group(size(group)))
    zero_group = 0.0_dp
    call optimize_ordinal_group(outcome, group, ints_alt, beta, ll_alt, st1)
    call optimize_ordinal_group(outcome, zero_group, ints_null, beta_null, ll_null, st0)
    if (st1 /= 0 .or. st0 /= 0) then
      status = 2
      return
    end if
    statistic = max(0.0_dp, 2.0_dp * (ll_alt - ll_null))
    p_value = erfc(sqrt(0.5_dp * statistic))
    status = 0
  end subroutine ord_test_po_numeric


  pure real(dp) function markov_target_error(intercepts, initial_offset, transition_offset, initial_state, &
                                              absorbing, target_rows, target) result(err)
    real(dp), intent(in) :: intercepts(:) !! Candidate ordered-logit intercepts.
    real(dp), intent(in) :: initial_offset(:) !! Initial-state cutpoint offsets used by exact occupancy propagation.
    real(dp), intent(in) :: transition_offset(:,:,:) !! Later-time offsets by time, previous state, and cutpoint.
    integer, intent(in) :: initial_state !! One-based conditioning state used by the occupancy propagation routine.
    logical, intent(in) :: absorbing(:) !! Absorbing-state flags.
    integer, intent(in) :: target_rows(:) !! One-based time rows whose occupancies are constrained.
    real(dp), intent(in) :: target(:,:) !! Target occupancy probabilities, row by constrained time and state.
    real(dp), allocatable :: probs(:,:)
    integer :: i, j, st

    call soprob_markov_ord(intercepts, initial_offset, transition_offset, initial_state, absorbing, probs, st)
    if (st /= 0) then
      err = huge(1.0_dp)
      return
    end if
    if (size(target,1) /= size(target_rows) .or. size(target,2) /= size(probs,2) .or. &
        any(target_rows < 1) .or. any(target_rows > size(probs,1))) then
      err = huge(1.0_dp)
      return
    end if
    err = 0.0_dp
    do i = 1, size(target_rows)
      do j = 1, size(probs,2)
        err = err + abs(probs(target_rows(i),j) - target(i,j))
      end do
    end do
  end function markov_target_error


  pure subroutine int_markov_ord_intercepts(initial_intercepts, initial_offset, transition_offset, initial_state, &
                                             absorbing, target_rows, target, intercepts, criterion, status)
    real(dp), intent(in) :: initial_intercepts(:) !! Starting descending ordered-logit intercepts.
    real(dp), intent(in) :: initial_offset(:) !! Cutpoint offsets for the first modeled time.
    real(dp), intent(in) :: transition_offset(:,:,:) !! Later-time offsets by time, previous state, and cutpoint.
    integer, intent(in) :: initial_state !! One-based initial conditioning state for occupancy propagation.
    logical, intent(in) :: absorbing(:) !! Absorbing-state flags aligned with outcome states.
    integer, intent(in) :: target_rows(:) !! One-based time rows where target occupancy probabilities are specified.
    real(dp), intent(in) :: target(:,:) !! Target occupancy probabilities by requested time row and state.
    real(dp), allocatable, intent(out) :: intercepts(:) !! Optimized descending intercept vector.
    real(dp), intent(out) :: criterion !! Final sum of absolute occupancy-probability errors.
    integer, intent(out) :: status !! Zero on success; nonzero for invalid inputs or occupancy propagation failure.
    real(dp), allocatable :: cand(:), step(:)
    real(dp) :: trial
    integer :: i, iter
    logical :: improved

    allocate(intercepts(size(initial_intercepts)), cand(size(initial_intercepts)), step(size(initial_intercepts)))
    intercepts = initial_intercepts
    if (size(intercepts) == 0 .or. size(target_rows) == 0 .or. &
        (size(intercepts) > 1 .and. any(intercepts(2:) > intercepts(:size(intercepts)-1)))) then
      criterion = nan_dp()
      status = 1
      return
    end if
    criterion = markov_target_error(intercepts, initial_offset, transition_offset, initial_state, absorbing, &
                                    target_rows, target)
    if (criterion >= 0.5_dp * huge(1.0_dp)) then
      status = 2
      return
    end if
    step = 0.5_dp
    do iter = 1, 500
      improved = .false.
      do i = 1, size(intercepts)
        cand = intercepts
        cand(i) = cand(i) + step(i)
        trial = markov_target_error(cand, initial_offset, transition_offset, initial_state, absorbing, &
                                    target_rows, target)
        if (trial < criterion) then
          intercepts = cand
          criterion = trial
          improved = .true.
        else
          cand = intercepts
          cand(i) = cand(i) - step(i)
          trial = markov_target_error(cand, initial_offset, transition_offset, initial_state, absorbing, &
                                      target_rows, target)
          if (trial < criterion) then
            intercepts = cand
            criterion = trial
            improved = .true.
          end if
        end if
      end do
      if (.not. improved) step = 0.5_dp * step
      if (maxval(step) < 1.0e-6_dp) exit
    end do
    status = 0
  end subroutine int_markov_ord_intercepts


  pure subroutine sim_reg_ord_supplied(treatment, outcomes, alpha, betas, statistics, p_values, power, status)
    real(dp), intent(in) :: treatment(:) !! Numeric treatment indicator/predictor aligned with rows of outcomes.
    integer, intent(in) :: outcomes(:,:) !! One-based simulated ordinal outcomes, observation by simulation replicate.
    real(dp), intent(in) :: alpha !! Two-sided likelihood-ratio significance threshold in the open interval (0,1).
    real(dp), allocatable, intent(out) :: betas(:) !! Fitted treatment proportional-odds slopes by replicate.
    real(dp), allocatable, intent(out) :: statistics(:) !! Likelihood-ratio chi-square statistics by replicate.
    real(dp), allocatable, intent(out) :: p_values(:) !! Likelihood-ratio P values by replicate.
    real(dp), intent(out) :: power !! Fraction of successfully fitted replicates with P < alpha.
    integer, intent(out) :: status !! Zero when at least one replicate fits; 1 invalid inputs, 2 no replicate fits.
    integer :: b, df, st, n_ok, n_sig

    allocate(betas(size(outcomes,2)), statistics(size(outcomes,2)), p_values(size(outcomes,2)))
    betas = nan_dp()
    statistics = nan_dp()
    p_values = nan_dp()
    power = nan_dp()
    if (size(outcomes,1) /= size(treatment) .or. size(treatment) < 3 .or. &
        size(outcomes,2) == 0 .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
      status = 1
      return
    end if
    n_ok = 0
    n_sig = 0
    do b = 1, size(outcomes,2)
      call ord_test_po_numeric(treatment, outcomes(:,b), statistics(b), df, p_values(b), betas(b), st)
      if (st == 0) then
        n_ok = n_ok + 1
        if (p_values(b) < alpha) n_sig = n_sig + 1
      end if
    end do
    if (n_ok == 0) then
      status = 2
      return
    end if
    power = real(n_sig,dp) / real(n_ok,dp)
    status = 0
  end subroutine sim_reg_ord_supplied

end module hmisc
