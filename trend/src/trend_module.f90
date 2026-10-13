! SPDX-License-Identifier: GPL-3.0-only
! Derived from trend 1.1.8, Copyright (C) 2015-2018 Thorsten Pohlert.
! Modern free-form translation: see NOTICE.md for provenance and changes.
module trend_module
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use r_kinds, only : dp
   use r_descriptive, only : r_mean, r_sd, r_correlation
   use r_distributions, only : r_pnorm, r_qnorm, r_pt
   use r_sorting, only : r_average_ranks, r_sort
   use r_quantiles, only : r_median
   use kendall_module, only : kendall_result, mann_kendall
   use randtests, only : randomness_result, bartels_rank_test, rt_success
   use exact_rank_tests, only : rank_test_result, wilcox_exact
   implicit none
   private

   integer, parameter, public :: trend_ok = 0
   integer, parameter, public :: trend_invalid = 1
   integer, parameter, public :: trend_insufficient = 2
   integer, parameter, public :: trend_nonfinite = 3
   integer, parameter, public :: trend_degenerate = 4
   integer, parameter, public :: trend_ambiguous = 5
   integer, parameter, public :: trend_limit = 6
   integer, parameter :: default_seed = 104729
   integer, parameter :: max_slope_n = 10000

   type, public :: trend_result
      !! Numerical test result; unavailable scalar fields are quiet NaNs after a call.
      real(dp) :: statistic = 0.0_dp !! z, t, RVN, U, or other documented test statistic.
      real(dp) :: p_value = 0.0_dp !! Tail probability, or NaN when undefined.
      real(dp) :: score = 0.0_dp !! Mann-Kendall score, adjusted score, or phase count.
      real(dp) :: variance = 0.0_dp !! Null variance when applicable.
      real(dp) :: tau = 0.0_dp !! Kendall tau-b when applicable.
      real(dp) :: estimate = 0.0_dp !! Sen slope or partial correlation when applicable.
      real(dp) :: standard_error = 0.0_dp !! ITAM slope standard deviation.
      real(dp) :: correlation = 0.0_dp !! Score correlation in the partial MK test.
      real(dp) :: conf_int(2) = 0.0_dp !! Lower and upper Sen confidence bounds.
      integer :: n = 0 !! Number of retained observations; rows for multivariate data.
      integer :: status = trend_invalid !! Zero success; public trend_* diagnostic otherwise.
      integer :: simulations = 0 !! Number of null samples used for Monte Carlo significance.
      logical :: exact = .false. !! True when Lanzante used exact Wilcoxon significance.
      real(dp), allocatable :: data(:) !! Change-point process in original observation order.
      integer, allocatable :: change_points(:) !! All maximizing one-based split indices.
      real(dp), allocatable :: season_score(:), season_variance(:), season_z(:), season_p(:), season_tau(:)
      real(dp), allocatable :: sorted_halves(:, :) !! ITAM sorted first and last halves, omitting an odd midpoint.
      real(dp), allocatable :: covariance(:, :) !! Multisite score covariance or partial-correlation matrix.
   end type trend_result

   public :: dp, itam, mk_test, smk_test, csmk_test, mult_mk_test, partial_mk_test
   public :: partial_cor_trend_test, sens_slope, sea_sens_slope
   public :: bartels_test, cs_test, wm_test, ww_test, rrod_test
   public :: pettitt_test, lanzante_test, br_test, bu_test, snh_test

contains

   pure subroutine initialize(result)
      type(trend_result), intent(out) :: result !! Reset result, deallocating prior arrays and marking undefined scalars.
      real(dp) :: missing
      missing = ieee_value(0.0_dp, ieee_quiet_nan)
      result%statistic = missing
      result%p_value = missing
      result%score = missing
      result%variance = missing
      result%tau = missing
      result%estimate = missing
      result%standard_error = missing
      result%correlation = missing
      result%conf_int = missing
   end subroutine initialize

   pure subroutine validate_vector(x, minimum_n, result)
      real(dp), intent(in) :: x(:) !! Observations after any documented missing-value omission.
      integer, intent(in) :: minimum_n !! Minimum required number of observations.
      type(trend_result), intent(inout) :: result !! Receives sample size and input validation status.
      result%n = size(x)
      if (size(x) < minimum_n) then
         result%status = trend_insufficient
      else if (any(.not. ieee_is_finite(x))) then
         result%status = trend_nonfinite
      else
         result%status = trend_ok
      end if
   end subroutine validate_vector

   pure integer function get_side(alternative) result(side)
      character(len=*), intent(in), optional :: alternative !! Exact spelling two.sided, greater, or less; default two.sided.
      side = 0
      if (.not. present(alternative)) return
      select case (alternative)
      case ('two.sided')
         side = 0
      case ('greater')
         side = 1
      case ('less')
         side = -1
      case default
         side = 2
      end select
   end function get_side

   pure elemental real(dp) function normal_probability(z, side) result(p)
      real(dp), intent(in) :: z !! Standard-normal test statistic; NaN propagates.
      integer, intent(in) :: side !! Zero two-sided, positive upper tail, negative lower tail.
      if (ieee_is_nan(z)) then
         p = z
         return
      end if
      if (side == 0) then
         p = min(1.0_dp, 2.0_dp * r_pnorm(abs(z), lower_tail = .false.))
      else
         p = r_pnorm(z, lower_tail = side < 0)
      end if
   end function normal_probability

   pure elemental real(dp) function compare_values(a, b) result(s)
      real(dp), intent(in) :: a !! First finite value to compare without subtractive overflow.
      real(dp), intent(in) :: b !! Second finite value; equality gives zero.
      s = 0.0_dp
      if (a > b) s = 1.0_dp
      if (a < b) s = -1.0_dp
   end function compare_values

   pure elemental real(dp) function standardized(score, variance, continuity) result(z)
      real(dp), intent(in) :: score !! Signed Mann-Kendall or aggregated score.
      real(dp), intent(in) :: variance !! Nonnegative null variance of the score.
      logical, intent(in) :: continuity !! Whether to subtract one from the absolute nonzero score.
      z = ieee_value(0.0_dp, ieee_quiet_nan)
      if (variance <= 0.0_dp) return
      z = score
      if (continuity) z = compare_values(score, 0.0_dp) * (abs(score) - 1.0_dp)
      z = z / sqrt(variance)
   end function standardized

   pure subroutine mk_components(x, score, variance, tau)
      real(dp), intent(in) :: x(:) !! Finite ordered values, at least two entries.
      real(dp), intent(out) :: score !! Signed Mann-Kendall pair count.
      real(dp), intent(out) :: variance !! Tie-corrected score variance.
      real(dp), intent(out) :: tau !! Tau-b against the untied observation indices; NaN for a constant series.
      type(kendall_result) :: shared
      if (size(x) >= 3) then
         call mann_kendall(x, shared)
         score = shared%s
         variance = shared%vars
         tau = shared%tau
      else
         score = compare_values(x(2), x(1))
         variance = abs(score)
         tau = score
         if (variance == 0.0_dp) tau = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
   end subroutine mk_components

   pure subroutine mk_test(x, result, alternative, continuity)
      real(dp), intent(in) :: x(:) !! Finite ordered observations, at least three; missing values are rejected.
      type(trend_result), intent(out) :: result !! z, p-value, score, tie-corrected variance, and tau-b.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater, or less.
      logical, intent(in), optional :: continuity !! Apply the upstream one-score continuity correction, default true.
      logical :: correct
      integer :: side
      call initialize(result)
      call validate_vector(x, 3, result)
      if (result%status /= trend_ok) return
      side = get_side(alternative)
      if (side == 2) then
         result%status = trend_invalid
         return
      end if
      correct = .true.
      if (present(continuity)) correct = continuity
      call mk_components(x, result%score, result%variance, result%tau)
      result%statistic = standardized(result%score, result%variance, correct)
      result%p_value = normal_probability(result%statistic, side)
      if (result%variance <= 0.0_dp) result%status = trend_degenerate
   end subroutine mk_test

   pure subroutine itam(x, result, conf_level)
      real(dp), intent(in) :: x(:) !! At least four finite ordered observations; odd midpoint omitted from the two halves.
      type(trend_result), intent(out) :: result !! ITAM slope, standard error, z, p-value, interval, and sorted halves.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      real(dp), allocatable :: first(:), last(:), scaled(:)
      real(dp) :: n, level, width, sigma, scale
      integer :: half
      call initialize(result)
      call validate_vector(x, 4, result)
      if (result%status /= trend_ok) return
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      if (.not. ieee_is_finite(level)) then
         result%status = trend_invalid
         return
      end if
      if (level <= 0.0_dp .or. level >= 1.0_dp) then
         result%status = trend_invalid
         return
      end if
      n = real(size(x), dp)
      half = size(x) / 2
      call r_sort(x(:half), first)
      call r_sort(x(size(x) - half + 1:), last)
      allocate(result%sorted_halves(half, 2))
      result%sorted_halves(:, 1) = first
      result%sorted_halves(:, 2) = last
      scaled = x
      scale = maxval(abs(x))
      if (scale > huge(1.0_dp) / 2.0_dp) then
         scaled = scaled / scale
         scaled = scaled - scaled(1)
      else
         scaled = scaled - scaled(1)
         scale = maxval(abs(scaled))
         if (scale > 0.0_dp) scaled = scaled / scale
      end if
      result%estimate = (2.0_dp * (r_mean(scaled(size(x) - half + 1:)) - r_mean(scaled(:half))) / n) * scale
      result%correlation = stable_correlation(first, last)
      if (.not. ieee_is_finite(result%correlation)) then
         result%status = trend_degenerate
         return
      end if
      sigma = (2.0_dp * sqrt(2.0_dp) / (n * sqrt(n)) * r_sd(scaled) * &
         sqrt(max(0.0_dp, 1.0_dp - result%correlation))) * scale
      result%standard_error = sigma
      width = r_qnorm((1.0_dp - level) / 2.0_dp, lower_tail = .false.) * sigma
      result%conf_int = [result%estimate - width, result%estimate + width]
      if (sigma == 0.0_dp) then
         if (result%estimate == 0.0_dp) then
            result%status = trend_degenerate
            return
         end if
         result%statistic = sign(ieee_value(0.0_dp, ieee_positive_inf), result%estimate)
      else
         result%statistic = result%estimate / sigma
      end if
      result%p_value = normal_probability(result%statistic, 0)
   end subroutine itam

   pure subroutine sens_slope(x, result, conf_level)
      real(dp), intent(in) :: x(:) !! Finite equally spaced observations, 3 to 10000; slopes are per observation interval.
      type(trend_result), intent(out) :: result !! Median pairwise slope, interval, and continuity-corrected MK significance.
      real(dp), intent(in), optional :: conf_level !! Confidence probability strictly between zero and one; default 0.95.
      real(dp), allocatable :: slopes(:), sorted(:)
      real(dp) :: level, critical
      integer :: i, j, k, lower, upper, n
      call mk_test(x, result)
      if (result%status /= trend_ok .and. result%status /= trend_degenerate) return
      n = size(x)
      if (n > max_slope_n) then
         result%status = trend_limit
         return
      end if
      level = 0.95_dp
      if (present(conf_level)) level = conf_level
      if (.not. ieee_is_finite(level)) then
         result%status = trend_invalid
         return
      end if
      if (level <= 0.0_dp .or. level >= 1.0_dp) then
         result%status = trend_invalid
         return
      end if
      allocate(slopes(n * (n - 1) / 2))
      k = 0
      do i = 1, n - 1
         do j = i + 1, n
            k = k + 1
            slopes(k) = (x(j) - x(i)) / real(j - i, dp)
         end do
      end do
      result%estimate = r_median(slopes)
      call r_sort(slopes, sorted)
      critical = r_qnorm((1.0_dp + level) / 2.0_dp) * sqrt(result%variance)
      lower = round_even((real(k, dp) - critical) / 2.0_dp)
      upper = round_even((real(k, dp) + critical) / 2.0_dp + 1.0_dp)
      if (lower >= 1 .and. lower <= k) result%conf_int(1) = sorted(lower)
      if (upper >= 1 .and. upper <= k) result%conf_int(2) = sorted(upper)
   end subroutine sens_slope

   pure integer function round_even(x) result(i)
      real(dp), intent(in) :: x !! Finite rank index within default integer range; half ties round to the even integer.
      integer :: base
      base = floor(x)
      i = base
      if (x - real(base, dp) > 0.5_dp) i = base + 1
      if (x - real(base, dp) == 0.5_dp .and. mod(base, 2) /= 0) i = base + 1
   end function round_even

   pure subroutine sea_sens_slope(x, period, result)
      real(dp), intent(in) :: x(:) !! Finite consecutive seasonal observations comprising complete cycles.
      integer, intent(in) :: period !! Seasons per cycle, at least two; each season needs at least two observations.
      type(trend_result), intent(out) :: result !! Pooled median within-season slope, in units per cycle.
      real(dp), allocatable :: slopes(:)
      integer :: n, i, j, s, k
      call initialize(result)
      call validate_vector(x, 4, result)
      if (result%status /= trend_ok) return
      if (period < 2 .or. period > size(x) / 2) then
         result%status = trend_invalid
         return
      end if
      if (mod(size(x), period) /= 0) then
         result%status = trend_invalid
         return
      end if
      n = size(x) / period
      if (n > max_slope_n .or. int(n, int64) * (n - 1) * period / 2 > 50000000_int64) then
         result%status = trend_limit
         return
      end if
      allocate(slopes(period * n * (n - 1) / 2))
      k = 0
      do s = 1, period
         do i = 1, n - 1
            do j = i + 1, n
               k = k + 1
               slopes(k) = (x(s + (j - 1) * period) - x(s + (i - 1) * period)) / real(j - i, dp)
            end do
         end do
      end do
      result%estimate = r_median(slopes)
   end subroutine sea_sens_slope

   pure subroutine smk_test(x, period, result, alternative, continuity)
      real(dp), intent(in) :: x(:) !! Finite consecutive seasonal observations; incomplete terminal cycles are supported.
      integer, intent(in) :: period !! At least two seasons, each with at least two observations.
      type(trend_result), intent(out) :: result !! Overall and per-season MK statistics in season order.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater, or less.
      logical, intent(in), optional :: continuity !! Correct both individual and overall scores, default true.
      integer :: s, side
      logical :: correct
      call initialize(result)
      call validate_vector(x, 4, result)
      if (result%status /= trend_ok) return
      side = get_side(alternative)
      if (period < 2 .or. period > size(x) / 2 .or. side == 2) then
         result%status = trend_invalid
         return
      end if
      correct = .true.
      if (present(continuity)) correct = continuity
      allocate(result%season_score(period), result%season_variance(period), result%season_z(period), &
         result%season_p(period), result%season_tau(period))
      do s = 1, period
         call mk_components(x(s::period), result%season_score(s), result%season_variance(s), result%season_tau(s))
         result%season_z(s) = standardized(result%season_score(s), result%season_variance(s), correct)
         result%season_p(s) = normal_probability(result%season_z(s), side)
      end do
      result%score = sum(result%season_score)
      result%variance = sum(result%season_variance)
      result%statistic = standardized(result%score, result%variance, correct)
      result%p_value = normal_probability(result%statistic, side)
      if (result%variance <= 0.0_dp) result%status = trend_degenerate
   end subroutine smk_test

   pure real(dp) function score_covariance(x, y) result(covariance)
      real(dp), intent(in) :: x(:) !! Finite first variate, in observation order.
      real(dp), intent(in) :: y(:) !! Finite second variate, same length as x.
      real(dp), allocatable :: rx(:), ry(:)
      real(dp) :: concordance, n
      integer :: i, j
      call r_average_ranks(x, rx)
      call r_average_ranks(y, ry)
      concordance = 0.0_dp
      do i = 1, size(x) - 1
         do j = i + 1, size(x)
            concordance = concordance + compare_values(x(j), x(i)) * compare_values(y(j), y(i))
         end do
      end do
      n = real(size(x), dp)
      covariance = (concordance + 4.0_dp * sum(rx * ry) - n * (n + 1.0_dp)**2) / 3.0_dp
   end function score_covariance

   pure subroutine mult_mk_test(x, result, alternative)
      real(dp), intent(in) :: x(:, :) !! Finite matrix of observations by sites; at least three rows and two columns.
      type(trend_result), intent(out) :: result !! Aggregate score, variance, significance, and site covariance matrix.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater, or less.
      integer :: i, j, side, n, sites
      real(dp) :: score, variance, tau
      call initialize(result)
      n = size(x, 1)
      sites = size(x, 2)
      result%n = n
      if (n < 3 .or. sites < 2) then
         result%status = trend_insufficient
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         result%status = trend_nonfinite
         return
      end if
      side = get_side(alternative)
      if (side == 2) return
      result%status = trend_ok
      allocate(result%covariance(sites, sites))
      result%score = 0.0_dp
      do i = 1, sites
         call mk_components(x(:, i), score, variance, tau)
         result%score = result%score + score
         do j = 1, i
            result%covariance(i, j) = score_covariance(x(:, i), x(:, j))
            result%covariance(j, i) = result%covariance(i, j)
         end do
      end do
      result%variance = sum(result%covariance)
      result%statistic = standardized(result%score, result%variance, n <= 10)
      result%p_value = normal_probability(result%statistic, side)
      if (result%variance <= 0.0_dp) result%status = trend_degenerate
   end subroutine mult_mk_test

   pure subroutine csmk_test(x, period, result, alternative)
      real(dp), intent(in) :: x(:) !! Finite consecutive seasonal data comprising complete cycles.
      integer, intent(in) :: period !! At least two seasons and at least three complete cycles.
      type(trend_result), intent(out) :: result !! Correlated seasonal score, covariance, and significance.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater, or less.
      real(dp), allocatable :: matrix(:, :)
      integer :: season, n
      call initialize(result)
      call validate_vector(x, 6, result)
      if (result%status /= trend_ok) return
      if (period < 2 .or. period > size(x) / 3) then
         result%status = trend_invalid
         return
      end if
      if (mod(size(x), period) /= 0) then
         result%status = trend_invalid
         return
      end if
      n = size(x) / period
      allocate(matrix(n, period))
      do season = 1, period
         matrix(:, season) = x(season::period)
      end do
      call mult_mk_test(matrix, result, alternative)
   end subroutine csmk_test

   pure subroutine partial_mk_test(x, y, result, alternative)
      real(dp), intent(in) :: x(:) !! Finite response observations, at least three, tested against time order.
      real(dp), intent(in) :: y(:) !! Finite covariate observations, same length as x.
      type(trend_result), intent(out) :: result !! Adjusted score, variance, score correlation, z, and significance.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater, or less.
      real(dp) :: sx, sy, unused, tau, n, base_variance, covariance
      integer :: side
      call initialize(result)
      call validate_vector(x, 3, result)
      if (result%status /= trend_ok) return
      if (size(y) /= size(x)) then
         result%status = trend_invalid
         return
      end if
      if (any(.not. ieee_is_finite(y))) then
         result%status = trend_nonfinite
         return
      end if
      side = get_side(alternative)
      if (side == 2) then
         result%status = trend_invalid
         return
      end if
      call mk_components(x, sx, unused, tau)
      call mk_components(y, sy, unused, tau)
      n = real(size(x), dp)
      base_variance = n * (n - 1.0_dp) * (2.0_dp * n + 5.0_dp) / 18.0_dp
      covariance = score_covariance(x, y)
      result%correlation = covariance / base_variance
      result%score = sx - result%correlation * sy
      result%variance = (1.0_dp - result%correlation**2) * base_variance
      result%statistic = standardized(result%score, result%variance, .false.)
      result%p_value = normal_probability(result%statistic, side)
      if (result%variance <= 0.0_dp) result%status = trend_degenerate
   end subroutine partial_mk_test

   pure real(dp) function stable_correlation(x, y) result(rho)
      real(dp), intent(in) :: x(:) !! Finite first variate; correlation is invariant to centering and positive scaling.
      real(dp), intent(in) :: y(:) !! Finite same-length second variate; constant inputs yield NaN.
      real(dp), allocatable :: a(:), b(:)
      real(dp) :: scale
      a = x
      scale = maxval(abs(a))
      if (scale > huge(1.0_dp) / 2.0_dp) a = a / scale
      a = a - a(1)
      scale = maxval(abs(a))
      if (scale > 0.0_dp) a = a / scale
      b = y
      scale = maxval(abs(b))
      if (scale > huge(1.0_dp) / 2.0_dp) b = b / scale
      b = b - b(1)
      scale = maxval(abs(b))
      if (scale > 0.0_dp) b = b / scale
      rho = r_correlation(a, b)
   end function stable_correlation

   pure subroutine partial_cor_trend_test(x, z, result, method)
      real(dp), intent(in) :: x(:) !! Finite response observations, at least three, tested against their positions.
      real(dp), intent(in) :: z(:) !! Finite covariate, same length as x.
      type(trend_result), intent(out) :: result !! Partial correlation, t statistic with n-2 degrees of freedom, and correlations.
      character(len=*), intent(in), optional :: method !! pearson (default) or spearman; no partial-name matching.
      real(dp), allocatable :: a(:), b(:), time(:)
      real(dp) :: numerator, denominator, rho, degrees
      integer :: i
      call initialize(result)
      call validate_vector(x, 3, result)
      if (result%status /= trend_ok) return
      if (size(z) /= size(x)) then
         result%status = trend_invalid
         return
      end if
      if (any(.not. ieee_is_finite(z))) then
         result%status = trend_nonfinite
         return
      end if
      a = x
      b = z
      if (present(method)) then
         select case (method)
         case ('pearson')
         case ('spearman')
            call r_average_ranks(x, a)
            call r_average_ranks(z, b)
         case default
            result%status = trend_invalid
            return
         end select
      end if
      time = [(real(i, dp), i = 1, size(x))]
      allocate(result%covariance(3, 3))
      result%covariance = 1.0_dp
      result%covariance(1, 2) = stable_correlation(time, a)
      result%covariance(1, 3) = stable_correlation(time, b)
      result%covariance(2, 3) = stable_correlation(a, b)
      result%covariance(2, 1) = result%covariance(1, 2)
      result%covariance(3, 1) = result%covariance(1, 3)
      result%covariance(3, 2) = result%covariance(2, 3)
      if (any(.not. ieee_is_finite(result%covariance))) then
         result%status = trend_degenerate
         return
      end if
      numerator = result%covariance(1, 2) - result%covariance(2, 3) * result%covariance(1, 3)
      denominator = sqrt(max(0.0_dp, 1.0_dp - result%covariance(2, 3)**2)) * &
         sqrt(max(0.0_dp, 1.0_dp - result%covariance(1, 3)**2))
      if (numerator == 0.0_dp .and. denominator == 0.0_dp) then
         rho = 0.0_dp
      else if (denominator <= 0.0_dp) then
         result%status = trend_degenerate
         return
      else
         rho = numerator / denominator
      end if
      if (abs(rho) > 1.0_dp + 1.0e-12_dp) then
         result%status = trend_degenerate
         return
      end if
      rho = max(-1.0_dp, min(1.0_dp, rho))
      result%estimate = rho
      degrees = real(size(x) - 2, dp)
      if (abs(rho) == 1.0_dp) then
         ! A perfect partial correlation has an infinite t statistic and zero p-value.
         result%statistic = sign(ieee_value(0.0_dp, ieee_positive_inf), rho)
         result%p_value = 0.0_dp
      else
         result%statistic = sqrt(degrees) * rho / sqrt(1.0_dp - rho**2)
         result%p_value = 2.0_dp * r_pt(abs(result%statistic), degrees, lower_tail = .false.)
      end if
   end subroutine partial_cor_trend_test

   pure subroutine bartels_test(x, result)
      real(dp), intent(in) :: x(:) !! At least ten finite observations; exact ties use average ranks.
      type(trend_result), intent(out) :: result !! Rank von Neumann ratio and upstream lower-tail significance.
      type(randomness_result) :: shared
      call initialize(result)
      call validate_vector(x, 10, result)
      if (result%status /= trend_ok) return
      shared = bartels_rank_test(x, alternative = 'left.sided', pvalue = 'beta')
      if (shared%status /= rt_success) then
         result%status = trend_degenerate
         return
      end if
      result%statistic = shared%rvn
      result%p_value = shared%p_value
      if (size(x) >= 100) result%p_value = r_pnorm(shared%rvn, location = 2.0_dp, &
         scale = sqrt(20.0_dp / (5.0_dp * real(size(x), dp) + 7.0_dp)))
   end subroutine bartels_test

   pure subroutine cs_test(x, result)
      real(dp), intent(in) :: x(:) !! Ordered data; NaNs omitted, infinities rejected, at least three retained values.
      type(trend_result), intent(out) :: result !! Trend's first/last-third Cox-Stuart z and two-sided significance.
      real(dp), allocatable :: a(:)
      real(dp) :: rn, correction
      integer :: n, length_third, greater, less
      call initialize(result)
      a = pack(x, .not. ieee_is_nan(x))
      call validate_vector(a, 3, result)
      if (result%status /= trend_ok) return
      n = size(a)
      length_third = ceiling(real(n, dp) / 3.0_dp)
      greater = count(a(n - length_third + 1:) > a(:length_third))
      less = count(a(n - length_third + 1:) < a(:length_third))
      if (greater + less == 0) then
         result%status = trend_degenerate
         return
      end if
      rn = real(n, dp)
      correction = 0.0_dp
      if (n <= 30) correction = 0.5_dp
      result%score = real(max(greater, less), dp)
      result%statistic = (abs(result%score - rn / 6.0_dp) - correction) / sqrt(rn / 12.0_dp)
      result%p_value = normal_probability(result%statistic, 0)
   end subroutine cs_test

   pure subroutine wm_test(x, result)
      real(dp), intent(in) :: x(:) !! Ordered data; NaNs omitted, zero adjacent differences skipped, infinities rejected.
      type(trend_result), intent(out) :: result !! Wallis-Moore phase-frequency z, p-value, and phase count in score.
      real(dp), allocatable :: a(:), signs(:)
      real(dp) :: rn, correction
      integer :: n, count_signs, i, phases
      call initialize(result)
      a = pack(x, .not. ieee_is_nan(x))
      call validate_vector(a, 2, result)
      if (result%status /= trend_ok) return
      n = size(a)
      signs = compare_values(a(2:), a(:n - 1))
      signs = pack(signs, signs /= 0.0_dp)
      count_signs = size(signs)
      phases = 0
      if (count_signs > 2) then
         do i = 2, count_signs
            if (signs(i) /= signs(i - 1)) phases = phases + 1
         end do
      end if
      phases = max(0, phases - 1)
      rn = real(n, dp)
      correction = 0.0_dp
      if (n <= 30) correction = 0.5_dp
      result%score = real(phases, dp)
      result%statistic = (abs(result%score - (2.0_dp * rn - 7.0_dp) / 3.0_dp) - correction) / &
         sqrt((16.0_dp * rn - 29.0_dp) / 90.0_dp)
      result%p_value = normal_probability(result%statistic, 0)
   end subroutine wm_test

   pure subroutine ww_test(x, result)
      real(dp), intent(in) :: x(:) !! Ordered data; NaNs omitted, infinities rejected; at least three retained values.
      type(trend_result), intent(out) :: result !! Wald-Wolfowitz cyclic-product z and two-sided significance.
      real(dp), allocatable :: a(:)
      real(dp) :: rn, s1, s2, s3, s4, expected, observed, scale
      integer :: n
      call initialize(result)
      a = pack(x, .not. ieee_is_nan(x))
      call validate_vector(a, 3, result)
      if (result%status /= trend_ok) return
      n = size(a)
      ! Centering and scaling preserve the statistic and stabilize the raw moments.
      scale = maxval(abs(a))
      if (scale == 0.0_dp .or. maxval(a) == minval(a)) then
         result%status = trend_degenerate
         return
      end if
      if (scale > huge(1.0_dp) / 2.0_dp) a = a / scale
      a = a - a(1)
      a = a / maxval(abs(a))
      if (n == 3) then
         result%status = trend_degenerate
         return
      end if
      rn = real(n, dp)
      observed = sum(a(:n - 1) * a(2:)) + a(1) * a(n)
      s1 = sum(a)
      s2 = sum(a**2)
      s3 = sum(a**3)
      s4 = sum(a**4)
      expected = (s1**2 - s2) / (rn - 1.0_dp)
      result%variance = (s2**2 - s4) / (rn - 1.0_dp) - expected**2 + &
         (s1**4 - 4.0_dp * s1**2 * s2 + 4.0_dp * s1 * s3 + s2**2 - 2.0_dp * s4) / &
         ((rn - 1.0_dp) * (rn - 2.0_dp))
      if (result%variance <= 0.0_dp) then
         result%status = trend_degenerate
         return
      end if
      result%statistic = (observed - expected) / sqrt(result%variance)
      result%p_value = normal_probability(result%statistic, 0)
   end subroutine ww_test

   pure subroutine rrod_test(x, y, result, alternative)
      real(dp), intent(in) :: x(:) !! First sample; all nonfinite values are omitted as in upstream rrod.test.default.
      real(dp), intent(in) :: y(:) !! Second sample; all nonfinite values are omitted independently.
      type(trend_result), intent(out) :: result !! Robust rank-order distributional z and selected-tail significance.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), greater, or less.
      real(dp), allocatable :: a(:), b(:), nx(:), ny(:)
      real(dp) :: mean_x, mean_y, denominator
      integer :: i, side
      call initialize(result)
      a = pack(x, ieee_is_finite(x))
      b = pack(y, ieee_is_finite(y))
      result%n = size(a) + size(b)
      if (min(size(a), size(b)) < 1) then
         result%status = trend_insufficient
         return
      end if
      side = get_side(alternative)
      if (side == 2) return
      result%status = trend_ok
      allocate(nx(size(a)), ny(size(b)))
      do i = 1, size(a)
         nx(i) = real(count(a(i) > b), dp)
      end do
      do i = 1, size(b)
         ny(i) = real(count(b(i) > a), dp)
      end do
      mean_x = r_mean(nx)
      mean_y = r_mean(ny)
      denominator = mean_x * mean_y + sum((nx - mean_x)**2) + sum((ny - mean_y)**2)
      if (denominator <= 0.0_dp) then
         if (sum(nx) == sum(ny)) then
            result%status = trend_degenerate
            return
         end if
         result%statistic = sign(ieee_value(0.0_dp, ieee_positive_inf), sum(nx) - sum(ny))
      else
         result%statistic = 0.5_dp * (sum(nx) - sum(ny)) / sqrt(denominator)
      end if
      result%p_value = normal_probability(result%statistic, side)
   end subroutine rrod_test

   pure subroutine pettitt_test(x, result)
      real(dp), intent(in) :: x(:) !! At least two finite observations in time order; ties use average ranks.
      type(trend_result), intent(out) :: result !! Pettitt statistic, approximate p-value, rank process, and maximizing splits.
      real(dp), allocatable :: ranks(:)
      real(dp) :: running, n
      integer :: i
      call initialize(result)
      call validate_vector(x, 2, result)
      if (result%status /= trend_ok) return
      call r_average_ranks(x, ranks)
      allocate(result%data(size(x)))
      n = real(size(x), dp)
      running = 0.0_dp
      do i = 1, size(x)
         running = running + ranks(i)
         result%data(i) = 2.0_dp * running - real(i, dp) * (n + 1.0_dp)
      end do
      result%statistic = maxval(abs(result%data))
      result%change_points = pack([(i, i = 1, size(x))], abs(result%data) == result%statistic)
      result%p_value = min(1.0_dp, 2.0_dp * exp(-6.0_dp * result%statistic**2 / (n**3 + n**2)))
   end subroutine pettitt_test

   pure subroutine lanzante_test(x, result, method)
      real(dp), intent(in) :: x(:) !! Finite ordered observations; a unique interior Pettitt maximizing split is required.
      type(trend_result), intent(out) :: result !! Split process and either Wilcoxon W or robust rank-order z with significance.
      character(len=*), intent(in), optional :: method !! wilcox.test (default) or rrod.test; R formula interfaces are omitted.
      type(trend_result) :: robust
      type(rank_test_result) :: wilcoxon
      real(dp), allocatable :: sorted(:)
      real(dp) :: adjusted
      logical :: untied
      integer :: split
      call pettitt_test(x, result)
      if (result%status /= trend_ok) return
      result%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      result%p_value = result%statistic
      if (size(result%change_points) /= 1) then
         result%status = trend_ambiguous
         return
      end if
      split = result%change_points(1)
      if (split >= size(x)) then
         result%status = trend_insufficient
         return
      end if
      if (present(method)) then
         select case (method)
         case ('wilcox.test')
         case ('rrod.test')
            call rrod_test(x(:split), x(split + 1:), robust)
            result%statistic = robust%statistic
            result%p_value = robust%p_value
            result%status = robust%status
            return
         case default
            result%status = trend_invalid
            return
         end select
      end if
      call r_sort(x, sorted)
      untied = all(sorted(2:) /= sorted(:size(x) - 1))
      call wilcox_exact(x(:split), wilcoxon, y = x(split + 1:), &
         exact = untied .and. split < 50 .and. size(x) - split < 50)
      result%statistic = wilcoxon%statistic
      result%p_value = wilcoxon%p_value
      result%exact = wilcoxon%exact
      if (wilcoxon%status /= 0) then
         result%status = trend_degenerate
         return
      end if
      if (.not. wilcoxon%exact) then
         if (wilcoxon%variance <= 0.0_dp) then
            result%status = trend_degenerate
            return
         end if
         adjusted = result%statistic - wilcoxon%expected
         adjusted = adjusted - 0.5_dp * compare_values(adjusted, 0.0_dp)
         result%p_value = normal_probability(adjusted / sqrt(wilcoxon%variance), 0)
      end if
   end subroutine lanzante_test

   pure subroutine br_test(x, result, m, seed, null_samples)
      real(dp), intent(in) :: x(:) !! At least two finite ordered observations with nonzero variance.
      type(trend_result), intent(out) :: result !! Buishand range statistic, cumulative process, splits, and Monte Carlo p-value.
      integer, intent(in), optional :: m !! Positive number of simulations; default 20000, or null_samples column count.
      integer, intent(in), optional :: seed !! Deterministic local RNG seed; default 104729, independent of global RANDOM_SEED.
      real(dp), intent(in), optional :: null_samples(:, :) !! Optional n-by-m finite normal samples from an external RNG.
      call changepoint_mc(x, 1, result, m, seed, null_samples)
   end subroutine br_test

   pure subroutine bu_test(x, result, m, seed, null_samples)
      real(dp), intent(in) :: x(:) !! At least two finite ordered observations with nonzero variance.
      type(trend_result), intent(out) :: result !! Buishand U statistic, scaled cumulative process, splits, and Monte Carlo p-value.
      integer, intent(in), optional :: m !! Positive simulation count; default 20000, or supplied null_samples column count.
      integer, intent(in), optional :: seed !! Deterministic local RNG seed, default 104729; does not affect global RNG state.
      real(dp), intent(in), optional :: null_samples(:, :) !! Optional finite n-by-m null realizations, one realization per column.
      call changepoint_mc(x, 2, result, m, seed, null_samples)
   end subroutine bu_test

   pure subroutine snh_test(x, result, m, seed, null_samples)
      real(dp), intent(in) :: x(:) !! At least two finite ordered observations with nonzero variance.
      type(trend_result), intent(out) :: result !! Standard normal homogeneity maximum, process, splits, and Monte Carlo p-value.
      integer, intent(in), optional :: m !! Positive simulation count; default 20000, or supplied null_samples column count.
      integer, intent(in), optional :: seed !! Deterministic local RNG seed, default 104729; does not affect global RNG state.
      real(dp), intent(in), optional :: null_samples(:, :) !! Optional finite n-by-m null realizations, one realization per column.
      call changepoint_mc(x, 3, result, m, seed, null_samples)
   end subroutine snh_test

   pure subroutine changepoint_mc(x, method, result, m, seed, null_samples)
      real(dp), intent(in) :: x(:) !! Ordered data, at least two finite observations.
      integer, intent(in) :: method !! One Buishand range, two Buishand U, three SNHT.
      type(trend_result), intent(out) :: result !! Observed test statistic, process, maximizing indices, and Monte Carlo tail.
      integer, intent(in), optional :: m !! Positive requested number of null realizations, default 20000.
      integer, intent(in), optional :: seed !! Any default integer, deterministically mapped to a valid Park-Miller state.
      real(dp), intent(in), optional :: null_samples(:, :) !! Optional n-row null samples; m must match their positive column count.
      real(dp), allocatable :: sample(:), process(:)
      real(dp) :: simulated, target
      integer(int64) :: state
      integer :: draws, i, j, exceed, n
      logical :: valid
      call initialize(result)
      call validate_vector(x, 2, result)
      if (result%status /= trend_ok) return
      n = size(x)
      draws = 20000
      if (present(null_samples)) then
         draws = size(null_samples, 2)
         if (size(null_samples, 1) /= n .or. any(.not. ieee_is_finite(null_samples))) then
            result%status = trend_invalid
            return
         end if
      end if
      if (present(m)) then
         if (present(null_samples)) then
            if (m /= draws) then
               result%status = trend_invalid
               return
            end if
         end if
         draws = m
      end if
      if (draws < 1) then
         result%status = trend_invalid
         return
      end if
      call change_statistic(x, method, .false., result%statistic, result%data, valid)
      if (.not. valid) then
         result%status = trend_degenerate
         return
      end if
      if (method == 3) then
         target = result%statistic
         result%change_points = pack([(i, i = 1, n - 1)], result%data == target)
      else
         target = maxval(abs(result%data))
         result%change_points = pack([(i, i = 1, n)], abs(result%data) == target)
      end if
      state = int(default_seed, int64)
      if (present(seed)) state = int(seed, int64)
      state = modulo(state - 1_int64, 2147483646_int64) + 1_int64
      allocate(sample(n))
      exceed = 0
      do j = 1, draws
         if (present(null_samples)) then
            sample = null_samples(:, j)
         else
            call normal_sample(state, sample)
         end if
         ! Upstream simulation uses divisor n, while observed R sd uses n-1.
         call change_statistic(sample, method, .true., simulated, process, valid)
         if (.not. valid) then
            result%status = trend_degenerate
            result%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (simulated > result%statistic) exceed = exceed + 1
      end do
      result%simulations = draws
      result%p_value = real(exceed, dp) / real(draws, dp)
   end subroutine changepoint_mc

   pure subroutine change_statistic(x, method, population, statistic, process, valid)
      real(dp), intent(in) :: x(:) !! Finite observed or null series, at least two values.
      integer, intent(in) :: method !! One range, two U, or three SNHT.
      logical, intent(in) :: population !! Use population variance divisor n for null simulations, otherwise n-1.
      real(dp), intent(out) :: statistic !! Corresponding scalar change-point statistic, or NaN for a constant series.
      real(dp), allocatable, intent(out) :: process(:) !! Scaled cumulative sums or SNHT process of length n or n-1.
      logical, intent(out) :: valid !! True for positive finite standard deviation and a finite statistic.
      real(dp) :: mean, sigma, running, n, total, scale
      real(dp), allocatable :: scaled(:)
      integer :: i, count_x
      count_x = size(x)
      n = real(count_x, dp)
      ! The three statistics are location/scale invariant. Center before computing
      ! moments to avoid cancellation for large offsets, and scale to avoid overflow.
      scaled = x
      scale = maxval(abs(scaled))
      if (scale > huge(1.0_dp) / 2.0_dp) scaled = scaled / scale
      scaled = scaled - scaled(1)
      scale = maxval(abs(scaled))
      if (scale > 0.0_dp) scaled = scaled / scale
      mean = r_mean(scaled)
      sigma = r_sd(scaled)
      if (population) sigma = sigma * sqrt((n - 1.0_dp) / n)
      valid = ieee_is_finite(sigma) .and. sigma > 0.0_dp
      statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. valid) return
      scaled = (scaled - mean) / sigma
      if (method == 3) then
         allocate(process(count_x - 1))
         running = 0.0_dp
         total = sum(scaled)
         do i = 1, count_x - 1
            running = running + scaled(i)
            process(i) = running**2 / real(i, dp) + (total - running)**2 / real(count_x - i, dp)
         end do
         statistic = maxval(process)
      else
         allocate(process(count_x))
         running = 0.0_dp
         do i = 1, count_x
            running = running + scaled(i)
            process(i) = running
         end do
         if (method == 1) then
            statistic = (maxval(process) - minval(process)) / sqrt(n)
         else
            statistic = sum(process(:count_x - 1)**2) / (n * (n + 1.0_dp))
         end if
      end if
      valid = ieee_is_finite(statistic)
   end subroutine change_statistic

   pure subroutine normal_sample(state, sample)
      integer(int64), intent(inout) :: state !! Park-Miller state in 1:2147483646, advanced by two uniforms per normal pair.
      real(dp), intent(out) :: sample(:) !! Standard-normal values generated with the Box-Muller transform.
      real(dp) :: u1, u2, radius, angle
      real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp
      integer :: i
      do i = 1, size(sample), 2
         state = modulo(16807_int64 * state, 2147483647_int64)
         u1 = real(state, dp) / 2147483647.0_dp
         state = modulo(16807_int64 * state, 2147483647_int64)
         u2 = real(state, dp) / 2147483647.0_dp
         radius = sqrt(-2.0_dp * log(u1))
         angle = two_pi * u2
         sample(i) = radius * cos(angle)
         if (i < size(sample)) sample(i + 1) = radius * sin(angle)
      end do
   end subroutine normal_sample

end module trend_module
