! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Computational translation of perm 1.0-0.4 by Michael P. Fay.
module perm_api
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_ordering, only : r_sort_values_in_place
   use r_distributions, only : r_pnorm, r_pchisq, r_qf
   use r_special, only : r_log_choose, r_regularized_beta
   use r_transforms, only : r_expm1, r_log1p
   implicit none
   private
   integer, parameter, public :: perm_ok = 0, perm_invalid = 1, perm_limit = 2, perm_allocation = 3
   integer, parameter, public :: perm_degenerate = 4, perm_numeric = 5
   integer(int64), parameter :: rng_modulus = 2147483647_int64
   integer(int64), parameter :: default_paths = 1000000_int64, default_cells = 10000000_int64
   integer(int64), parameter :: default_work = 100000000_int64
   integer, parameter :: max_observations = 100000, max_draws = 10000000
   type, public :: perm_rng
      integer(int64) :: state = 1234321_int64
   end type perm_rng
   type, public :: perm_control_type
      integer :: nmc = 999
      integer :: digits = 12
      integer(int64) :: seed = 1234321_int64
      real(dp) :: confidence = 0.99_dp
      logical :: set_seed = .true.
      character(len=7) :: tsmethod = 'central'
      integer(int64) :: max_permutations = default_paths
      integer(int64) :: max_cells = default_cells
      integer(int64) :: max_work = default_work
      integer, allocatable :: cm(:, :)
   end type perm_control_type
   type, public :: perm_result
      real(dp) :: p_twosided = -1.0_dp
      real(dp) :: p_absolute = -1.0_dp
      real(dp) :: p_less = -1.0_dp
      real(dp) :: p_greater = -1.0_dp
      real(dp) :: p_equal = -1.0_dp
      real(dp) :: p_value = -1.0_dp
      real(dp) :: statistic = 0.0_dp
      real(dp) :: estimate = 0.0_dp
      real(dp) :: confidence_interval(2) = [-1.0_dp, -1.0_dp]
      integer :: df = 0
      integer :: nmc = 0
      integer :: status = perm_invalid
      character(len=16) :: method = ''
   end type perm_result
   public :: dp, perm_control, choose_matrix, getcnt, calc_pvals, calc_pvals_mc, p_ci
   public :: method_rule_ts1, method_rule_ks1, method_rule_trend1
   public :: twosample_pclt, ksample_pclt, trend_pclt
   public :: twosample_exact_ce, twosample_exact_network
   public :: twosample_exact_mc, ksample_exact_mc, trend_exact_mc
   public :: perm_ts, perm_ks, perm_trend
contains
   pure subroutine perm_control(control, status, nmc, seed, digits, confidence, set_seed, tsmethod, &
      cm, max_permutations, max_cells, max_work)
      type(perm_control_type), intent(out) :: control !! Validated options, with R defaults unless overridden.
      integer, intent(out) :: status !! perm_ok on success, otherwise perm_invalid or perm_limit.
      integer, intent(in), optional :: nmc !! Number of random permutations, from 1 through 10,000,000; default 999.
      integer(int64), intent(in), optional :: seed !! Local generator seed; any int64 value is deterministically normalized.
      integer, intent(in), optional :: digits !! Significant digits for centered-statistic comparisons, 1 through 15; default 12.
      real(dp), intent(in), optional :: confidence !! Monte Carlo interval confidence in [0,1], default 0.99.
      logical, intent(in), optional :: set_seed !! Reset local RNG for each call; false requires an explicit perm_rng argument.
      character(len=*), intent(in), optional :: tsmethod !! Two-sided convention: central (default) or abs.
      integer, intent(in), optional :: cm(:, :) !! Optional zero/one selection rows for complete enumeration.
      integer(int64), intent(in), optional :: max_permutations !! Exact enumeration path cap, default 1,000,000.
      integer(int64), intent(in), optional :: max_cells !! Maximum cells for allocated tables, default 10,000,000.
      integer(int64), intent(in), optional :: max_work !! Maximum estimated enumeration or Monte Carlo operations, default 100,000,000.
      if (present(nmc)) control%nmc = nmc
      if (present(seed)) control%seed = seed
      if (present(digits)) control%digits = digits
      if (present(confidence)) control%confidence = confidence
      if (present(set_seed)) control%set_seed = set_seed
      if (present(tsmethod)) then
         if (tsmethod /= 'central' .and. tsmethod /= 'abs') then
            status = perm_invalid
            return
         end if
         control%tsmethod = tsmethod
      end if
      if (present(max_permutations)) control%max_permutations = max_permutations
      if (present(max_cells)) control%max_cells = max_cells
      if (present(max_work)) control%max_work = max_work
      status = control_status(control)
      if (status /= perm_ok) return
      if (present(cm)) then
         if (size(cm, kind=int64) > control%max_cells) then
            status = perm_limit
            return
         end if
         if (any(cm /= 0 .and. cm /= 1)) then
            status = perm_invalid
            return
         end if
         control%cm = cm
      end if
   end subroutine perm_control

   pure integer function control_status(control) result(status)
      type(perm_control_type), intent(in) :: control !! Options to validate before any array allocation or random draws.
      status = perm_invalid
      if (control%nmc < 1 .or. control%digits < 1 .or. control%digits > 15) return
      if (.not. ieee_is_finite(control%confidence)) return
      if (control%confidence < 0.0_dp .or. control%confidence > 1.0_dp) return
      if (control%tsmethod /= 'central' .and. control%tsmethod /= 'abs') return
      if (control%max_permutations < 1_int64 .or. control%max_cells < 1_int64) return
      if (control%max_work < 1_int64) return
      status = perm_limit
      if (control%nmc > max_draws) return
      if (control%max_permutations >= int(huge(1), int64)) return
      status = perm_ok
   end function control_status

   pure integer function data_status(scores) result(status)
      real(dp), intent(in) :: scores(:) !! Finite scores, with magnitude bounded to prevent subsequent quadratic overflow.
      status = perm_invalid
      if (size(scores) < 2) return
      if (any(.not. ieee_is_finite(scores))) return
      status = perm_limit
      if (size(scores) > max_observations) return
      status = perm_numeric
      if (maxval(abs(scores)) > sqrt(huge(1.0_dp)) / (8.0_dp * real(size(scores), dp))) return
      status = perm_ok
   end function data_status

   pure elemental real(dp) function significant(x, digits) result(value)
      real(dp), intent(in) :: x !! Finite centered statistic; zeros are returned unchanged.
      integer, intent(in) :: digits !! Decimal significant digits, restricted by callers to 1 through 15.
      real(dp) :: power, scaled, whole, fraction, magnitude, factor
      integer :: exponent10
      value = x
      if (x == 0.0_dp) return
      factor = 1.0_dp
      if (abs(x) < 1.0e-200_dp) factor = 1.0e200_dp
      if (abs(x) > 1.0e200_dp) factor = 1.0e-200_dp
      magnitude = abs(x) * factor
      exponent10 = floor(log10(magnitude)) - digits + 1
      power = 10.0_dp ** exponent10
      scaled = magnitude / power
      if (scaled >= 2.0_dp ** 52) return
      whole = real(floor(scaled, kind=int64), dp)
      fraction = scaled - whole
      if (fraction > 0.5_dp .or. (fraction == 0.5_dp .and. modulo(whole, 2.0_dp) /= 0.0_dp)) whole = whole + 1.0_dp
      ! Canonicalize a rounded carry across a decimal decade; otherwise equal
      ! decimal answers can differ by an ulp because their powers differ.
      if (whole >= 10.0_dp ** digits) then
         whole = whole / 10.0_dp
         power = 10.0_dp ** (exponent10 + 1)
      end if
      value = sign((whole * power) / factor, x)
   end function significant

   pure subroutine calc_pvals(statistics, observed, result, digits, weights, two_sided_statistic)
      real(dp), intent(in) :: statistics(:) !! Enumerated finite statistics; at least one value is required.
      real(dp), intent(in) :: observed !! Observed statistic, compared after centering and significant-digit rounding.
      type(perm_result), intent(out) :: result !! All five tail probabilities, or the upper-tail two-sided statistic probability.
      integer, intent(in), optional :: digits !! Significant digits, 1 through 15; default 12.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative finite enumeration weights matching statistics; sum must be positive.
      logical, intent(in), optional :: two_sided_statistic !! True uses the upper tail for an already two-sided statistic.
      real(dp), allocatable :: probability(:), centered(:)
      real(dp) :: mean_value, t0, magnitude
      integer :: d, astat
      logical :: upper_only
      result%status = perm_invalid
      if (size(statistics) < 1 .or. .not. ieee_is_finite(observed)) return
      if (any(.not. ieee_is_finite(statistics))) return
      d = 12
      if (present(digits)) d = digits
      if (d < 1 .or. d > 15) return
      upper_only = .false.
      if (present(two_sided_statistic)) upper_only = two_sided_statistic
      if (max(maxval(abs(statistics)), abs(observed)) > huge(1.0_dp) / 8.0_dp) then
         result%status = perm_numeric
         return
      end if
      allocate(probability(size(statistics)), centered(size(statistics)), stat=astat)
      if (astat /= 0) then
         result%status = perm_allocation
         return
      end if
      probability = 1.0_dp / real(size(statistics), dp)
      if (present(weights)) then
         if (size(weights) /= size(statistics)) return
         if (any(.not. ieee_is_finite(weights)) .or. any(weights < 0.0_dp)) return
         magnitude = maxval(weights)
         if (magnitude <= 0.0_dp) return
         probability = weights / magnitude
         probability = probability / sum(probability)
      end if
      centered = statistics
      t0 = observed
      magnitude = max(maxval(abs(statistics)), abs(observed))
      if (magnitude > 0.0_dp .and. magnitude < 1.0e-200_dp) then
         centered = scale(statistics, -exponent(magnitude))
         t0 = scale(observed, -exponent(magnitude))
      end if
      mean_value = sum(probability * centered)
      centered = significant(centered - mean_value, d)
      t0 = significant(t0 - mean_value, d)
      result%p_less = min(1.0_dp, sum(probability, mask=centered <= t0))
      result%p_greater = min(1.0_dp, sum(probability, mask=centered >= t0))
      result%p_equal = min(1.0_dp, sum(probability, mask=centered == t0))
      result%p_absolute = min(1.0_dp, sum(probability, mask=abs(centered) >= abs(t0)))
      result%p_twosided = min(1.0_dp, 2.0_dp * min(result%p_less, result%p_greater))
      if (upper_only) result%p_twosided = result%p_greater
      result%p_value = result%p_twosided
      result%statistic = observed
      result%status = perm_ok
   end subroutine calc_pvals

   pure subroutine p_ci(successes, trials, interval, status, confidence)
      integer, intent(in) :: successes !! Number of tail hits, 0 through trials, without the added observed permutation.
      integer, intent(in) :: trials !! Number of independent Monte Carlo permutations, 1 through 10,000,000.
      real(dp), intent(out) :: interval(2) !! Equal-tailed Clopper-Pearson lower and upper probability bounds.
      integer, intent(out) :: status !! perm_ok on success; perm_invalid for an invalid count or confidence.
      real(dp), intent(in), optional :: confidence !! Coverage probability in [0,1], default 0.99.
      real(dp) :: level, alpha, a, b
      interval = ieee_value(0.0_dp, ieee_quiet_nan)
      status = perm_invalid
      if (trials < 1 .or. successes < 0 .or. successes > trials) return
      if (trials > max_draws) then
         status = perm_limit
         return
      end if
      level = 0.99_dp
      if (present(confidence)) level = confidence
      if (.not. ieee_is_finite(level)) return
      if (level < 0.0_dp .or. level > 1.0_dp) return
      interval = [0.0_dp, 1.0_dp]
      alpha = (1.0_dp - level) / 2.0_dp
      if (alpha > 0.0_dp) then
         ! Beta(a,b) = a*F(2a,2b)/(b+a*F). Both inversions use the small tail.
         if (successes > 0) then
            a = real(successes, dp)
            b = real(trials - successes + 1, dp)
            interval(1) = beta_lower_quantile(alpha, a, b)
         end if
         if (successes < trials) then
            a = real(trials - successes, dp)
            b = real(successes + 1, dp)
            interval(2) = 1.0_dp - beta_lower_quantile(alpha, a, b)
         end if
      end if
      status = perm_ok
      if (any(.not. ieee_is_finite(interval))) status = perm_numeric
   end subroutine p_ci

   pure elemental real(dp) function beta_lower_quantile(probability, a, b) result(value)
      real(dp), intent(in) :: probability !! Lower-tail probability in (0,0.5], as used by a p-value interval.
      real(dp), intent(in) :: a !! Positive integer-valued first beta shape.
      real(dp), intent(in) :: b !! Positive integer-valued second beta shape.
      real(dp) :: f, lo, hi, mid, cdf
      integer :: iteration
      if (a == 1.0_dp) then
         value = -r_expm1(r_log1p(-probability) / b)
         return
      end if
      if (b == 1.0_dp) then
         value = exp(log(probability) / a)
         return
      end if
      f = r_qf(probability, 2.0_dp * a, 2.0_dp * b)
      value = f / (b / a + f)
      lo = 0.0_dp
      hi = 1.0_dp
      cdf = r_regularized_beta(value, a, b)
      if (cdf < probability) then
         lo = value
      else
         hi = value
      end if
      ! The shared F inverse uses an absolute tolerance. Relative beta refinement
      ! preserves tiny nonzero Clopper-Pearson bounds at confidence near one.
      do iteration = 1, 200
         mid = lo + (hi - lo) / 2.0_dp
         cdf = r_regularized_beta(mid, a, b)
         if (cdf < probability) then
            lo = mid
         else
            hi = mid
         end if
         if (hi - lo <= 8.0_dp * epsilon(1.0_dp) * max(mid, tiny(1.0_dp))) exit
      end do
      value = lo + (hi - lo) / 2.0_dp
   end function beta_lower_quantile

   pure subroutine calc_pvals_mc(statistics, observed, result, digits, alternative, two_sided_statistic, confidence)
      real(dp), intent(in) :: statistics(:) !! Fixed-count simulated statistics, excluding the observed permutation.
      real(dp), intent(in) :: observed !! Observed finite test statistic.
      type(perm_result), intent(out) :: result !! Add-one tail probabilities and the selected tail's Monte Carlo confidence interval.
      integer, intent(in), optional :: digits !! Significant digits for centered comparisons, default 12.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), less, greater, or two.sidedAbs.
      logical, intent(in), optional :: two_sided_statistic !! True selects the upper tail of an already two-sided statistic.
      real(dp), intent(in), optional :: confidence !! Clopper-Pearson confidence in [0,1], default 0.99.
      real(dp), allocatable :: centered(:)
      real(dp) :: mean_value, t0, denom, magnitude
      integer :: d, n, lower, upper, absolute, equal, hits, status
      logical :: upper_only
      character(len=16) :: alt
      call calc_pvals(statistics, observed, result, digits=digits, two_sided_statistic=two_sided_statistic)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      n = size(statistics)
      if (n > max_draws) then
         result%status = perm_limit
         return
      end if
      d = 12
      if (present(digits)) d = digits
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) then
            result%status = perm_invalid
            return
         end if
         alt = alternative
      end if
      upper_only = .false.
      if (present(two_sided_statistic)) upper_only = two_sided_statistic
      if (.not. valid_alternative(alt) .and. .not. upper_only) return
      centered = statistics
      t0 = observed
      magnitude = max(maxval(abs(statistics)), abs(observed))
      if (magnitude > 0.0_dp .and. magnitude < 1.0e-200_dp) then
         centered = scale(statistics, -exponent(magnitude))
         t0 = scale(observed, -exponent(magnitude))
      end if
      mean_value = sum(centered / real(n, dp))
      centered = significant(centered - mean_value, d)
      t0 = significant(t0 - mean_value, d)
      lower = count(centered <= t0)
      upper = count(centered >= t0)
      absolute = count(abs(centered) >= abs(t0))
      equal = count(centered == t0)
      denom = real(n, dp) + 1.0_dp
      result%p_less = real(lower + 1, dp) / denom
      result%p_greater = real(upper + 1, dp) / denom
      result%p_equal = real(equal + 1, dp) / denom
      result%p_absolute = real(absolute + 1, dp) / denom
      result%p_twosided = min(1.0_dp, 2.0_dp * min(result%p_less, result%p_greater))
      hits = min(lower, upper)
      if (alt == 'less') hits = lower
      if (alt == 'greater') hits = upper
      if (alt == 'two.sidedAbs') hits = absolute
      if (upper_only) then
         hits = upper
         result%p_twosided = result%p_greater
      end if
      call p_ci(hits, n, result%confidence_interval, status, confidence)
      if (status /= perm_ok) then
         result%status = status
         return
      end if
      if (alt == 'two.sided' .and. .not. upper_only) then
         result%confidence_interval = 2.0_dp * result%confidence_interval
         result%confidence_interval = min(1.0_dp, result%confidence_interval)
      end if
      call select_probability(result, alt)
      if (upper_only) result%p_value = result%p_twosided
      result%nmc = n
      result%status = perm_ok
   end subroutine calc_pvals_mc

   pure logical function valid_alternative(alternative) result(valid)
      character(len=*), intent(in) :: alternative !! One of the four recognized R tail names.
      valid = alternative == 'two.sided' .or. alternative == 'less' .or. &
         alternative == 'greater' .or. alternative == 'two.sidedAbs'
   end function valid_alternative

   pure subroutine select_probability(result, alternative)
      type(perm_result), intent(inout) :: result !! Tail probabilities from a completed permutation calculation.
      character(len=*), intent(in) :: alternative !! Requested tail name.
      select case (alternative)
      case ('less')
         result%p_value = result%p_less
      case ('greater')
         result%p_value = result%p_greater
      case ('two.sidedAbs')
         result%p_value = result%p_absolute
      case default
         result%p_value = result%p_twosided
      end select
   end subroutine select_probability

   pure integer(int64) function bounded_choose(n, m, cap) result(value)
      integer, intent(in) :: n !! Nonnegative population size, at most 100,000.
      integer, intent(in) :: m !! Selection size between zero and n.
      integer(int64), intent(in) :: cap !! Saturation bound below the default-integer maximum.
      integer :: i, k
      value = 0_int64
      if (n < 0 .or. m < 0 .or. m > n) return
      value = 1_int64
      k = min(m, n - m)
      do i = 1, k
         value = value * int(n - k + i, int64) / int(i, int64)
         if (value > cap) then
            value = cap + 1_int64
            return
         end if
      end do
   end function bounded_choose

   pure subroutine choose_matrix(n, m, matrix, status, max_permutations, max_cells)
      integer, intent(in) :: n !! Population size from zero through 100,000.
      integer, intent(in) :: m !! Number selected, from zero through n.
      integer, allocatable, intent(out) :: matrix(:, :) !! All binary selection rows, ordered as R chooseMatrix.
      integer, intent(out) :: status !! perm_ok, perm_invalid, perm_limit, or perm_allocation.
      integer(int64), intent(in), optional :: max_permutations !! Row count cap, default 1,000,000.
      integer(int64), intent(in), optional :: max_cells !! Total matrix element cap, default 10,000,000.
      integer(int64) :: cap, cells, rows
      integer, allocatable :: index(:)
      integer :: k, i, row, j, astat
      status = perm_invalid
      if (n < 0 .or. m < 0 .or. m > n) return
      cap = default_paths
      cells = default_cells
      if (present(max_permutations)) cap = max_permutations
      if (present(max_cells)) cells = max_cells
      if (cap < 1_int64 .or. cap >= int(huge(1), int64) .or. cells < 1_int64) return
      status = perm_limit
      if (n > max_observations) return
      rows = bounded_choose(n, m, cap)
      if (rows > cap) return
      if (rows > cells / max(1, n)) return
      k = min(m, n - m)
      allocate(matrix(int(rows), n), index(k), stat=astat)
      if (astat /= 0) then
         status = perm_allocation
         return
      end if
      matrix = 0
      index = [(i, i=1,k)]
      do row = 1, int(rows)
         matrix(row, index) = 1
         if (row == rows) exit
         j = k
         do while (j > 0)
            if (index(j) < n - k + j) exit
            j = j - 1
         end do
         index(j) = index(j) + 1
         do i = j + 1, k
            index(i) = index(i - 1) + 1
         end do
      end do
      if (m > n - m) matrix = 1 - matrix
      status = perm_ok
   end subroutine choose_matrix

   pure subroutine getcnt(nodes, edge_start, edge_size, indices, status, max_cells)
      integer, intent(in) :: nodes(:) !! One-based network node indices, repetitions permitted.
      integer, intent(in) :: edge_start(:) !! One-based first edge of each node and a final sentinel.
      integer, intent(in) :: edge_size(:) !! Number of outgoing edges per node, consistent with adjacent edge starts.
      integer, allocatable, intent(out) :: indices(:) !! Concatenated outgoing-edge indices in node order.
      integer, intent(out) :: status !! perm_ok on success or an invalid/network allocation limit status.
      integer(int64), intent(in), optional :: max_cells !! Maximum returned indices, default 10,000,000.
      integer(int64) :: cap, total, width
      integer :: i, node, position, j, astat
      cap = default_cells
      if (present(max_cells)) cap = max_cells
      status = perm_invalid
      if (cap < 1_int64 .or. size(edge_start) /= size(edge_size) + 1) return
      if (any(nodes < 1) .or. any(nodes > size(edge_size))) return
      if (any(edge_start < 1) .or. any(edge_size < 0)) return
      if (any(edge_start(2:) - edge_start(:size(edge_size)) /= edge_size)) return
      total = 0_int64
      do i = 1, size(nodes)
         width = int(edge_size(nodes(i)), int64)
         if (width > cap - total) then
            status = perm_limit
            return
         end if
         total = total + width
      end do
      if (total >= int(huge(1), int64)) then
         status = perm_limit
         return
      end if
      allocate(indices(int(total)), stat=astat)
      if (astat /= 0) then
         status = perm_allocation
         return
      end if
      position = 0
      do i = 1, size(nodes)
         node = nodes(i)
         do j = edge_start(node), edge_start(node + 1) - 1
            position = position + 1
            indices(position) = j
         end do
      end do
      status = perm_ok
   end subroutine getcnt

   pure subroutine label_groups(group, labels, counts, status)
      integer, intent(in) :: group(:) !! Arbitrary integer group labels, interpreted in increasing numeric order.
      integer, allocatable, intent(out) :: labels(:) !! Dense one-based group indices corresponding to observations.
      integer, allocatable, intent(out) :: counts(:) !! Positive sizes of the distinct groups.
      integer, intent(out) :: status !! perm_ok or a size/allocation error.
      real(dp), allocatable :: values(:)
      integer, allocatable :: sizes(:)
      integer :: i, j, ng, n, astat
      n = size(group)
      status = perm_invalid
      if (n < 2) return
      status = perm_limit
      if (n > max_observations) return
      allocate(values(n), labels(n), sizes(n), stat=astat)
      if (astat /= 0) then
         status = perm_allocation
         return
      end if
      values = real(group, dp)
      call r_sort_values_in_place(values)
      ng = 1
      do i = 2, n
         if (values(i) /= values(ng)) then
            ng = ng + 1
            values(ng) = values(i)
         end if
      end do
      sizes = 0
      do i = 1, n
         ! Binary search preserves O(n log(number of groups)) grouping cost.
         block
            integer :: lo, hi
            lo = 1
            hi = ng
            do while (lo < hi)
               j = lo + (hi - lo) / 2
               if (real(group(i), dp) > values(j)) then
                  lo = j + 1
               else
                  hi = j
               end if
            end do
            labels(i) = lo
            sizes(lo) = sizes(lo) + 1
         end block
      end do
      counts = sizes(:ng)
      status = perm_ok
   end subroutine label_groups

   pure subroutine trend_pclt(scores, group, result)
      real(dp), intent(in) :: scores(:) !! Finite response scores, length at least two.
      real(dp), intent(in) :: group(:) !! Finite numerical trend covariate of equal length.
      type(perm_result), intent(out) :: result !! Standard-normal approximation, correlation estimate and both tails.
      real(dp), allocatable :: x(:), y(:)
      real(dp) :: sx, sy, ssx, ssy, correlation
      integer :: status, n
      status = data_status(scores)
      result%status = status
      if (status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      status = data_status(group)
      result%status = status
      if (status /= perm_ok) return
      n = size(scores)
      sx = maxval(abs(scores))
      sy = maxval(abs(group))
      result%status = perm_degenerate
      if (sx == 0.0_dp .or. sy == 0.0_dp) return
      x = scores / sx
      y = group / sy
      x = x - sum(x) / real(n, dp)
      y = y - sum(y) / real(n, dp)
      ssx = sum(x * x)
      ssy = sum(y * y)
      if (ssx == 0.0_dp .or. ssy == 0.0_dp) return
      correlation = sum(x * y) / sqrt(ssx) / sqrt(ssy)
      result%estimate = max(-1.0_dp, min(1.0_dp, correlation))
      result%statistic = sqrt(real(n - 1, dp)) * result%estimate
      result%p_less = r_pnorm(result%statistic)
      result%p_greater = r_pnorm(result%statistic, lower_tail=.false.)
      result%p_twosided = min(1.0_dp, 2.0_dp * min(result%p_less, result%p_greater))
      result%p_absolute = result%p_twosided
      result%p_value = result%p_twosided
      result%p_equal = 0.0_dp
      result%method = 'pclt'
      result%status = perm_ok
   end subroutine trend_pclt

   pure subroutine twosample_pclt(scores, group, result)
      real(dp), intent(in) :: scores(:) !! Finite response scores, one per group label.
      integer, intent(in) :: group(:) !! Exactly two integer labels; the larger label defines the upper-tail sample.
      type(perm_result), intent(out) :: result !! Permutation CLT normal statistic and tail probabilities.
      integer, allocatable :: labels(:), counts(:)
      integer :: status
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      call label_groups(group, labels, counts, status)
      result%status = status
      if (status /= perm_ok) return
      result%status = perm_invalid
      if (size(counts) /= 2) return
      call trend_pclt(scores, real(labels - 1, dp), result)
   end subroutine twosample_pclt

   pure real(dp) function k_statistic(centered, labels, counts) result(statistic)
      real(dp), intent(in) :: centered(:) !! Centered scores whose sum is numerically zero.
      integer, intent(in) :: labels(:) !! Dense group indices, one per centered score.
      integer, intent(in) :: counts(:) !! Positive group sizes indexed by labels.
      real(dp) :: totals(size(counts))
      integer :: i
      totals = 0.0_dp
      do i = 1, size(centered)
         totals(labels(i)) = totals(labels(i)) + centered(i)
      end do
      statistic = sum(totals * totals / real(counts, dp))
   end function k_statistic

   pure subroutine ksample_pclt(scores, group, result)
      real(dp), intent(in) :: scores(:) !! Finite response scores; at least two groups are needed.
      integer, intent(in) :: group(:) !! Integer categorical labels matching scores, with two or more distinct values.
      type(perm_result), intent(out) :: result !! Chi-square approximation and group-count-minus-one degrees of freedom.
      integer, allocatable :: labels(:), counts(:)
      real(dp), allocatable :: centered(:)
      real(dp) :: scale, variance
      integer :: status, n
      result%status = data_status(scores)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      call label_groups(group, labels, counts, status)
      result%status = status
      if (status /= perm_ok) return
      result%status = perm_invalid
      if (size(counts) < 2) return
      result%status = perm_degenerate
      scale = maxval(abs(scores))
      if (scale == 0.0_dp) return
      n = size(scores)
      centered = scores / scale
      centered = centered - sum(centered) / real(n, dp)
      variance = sum(centered * centered)
      if (variance == 0.0_dp) return
      result%statistic = real(n - 1, dp) * k_statistic(centered, labels, counts) / variance
      result%df = size(counts) - 1
      result%p_twosided = r_pchisq(result%statistic, real(result%df, dp), lower_tail=.false.)
      result%p_value = result%p_twosided
      result%method = 'pclt'
      result%status = perm_ok
   end subroutine ksample_pclt

   pure subroutine twosample_exact_ce(scores, group, result, control)
      real(dp), intent(in) :: scores(:) !! Finite scores; arbitrary real values and ties are supported.
      integer, intent(in) :: group(:) !! Exactly two labels, with the larger label defining the sampled group.
      type(perm_result), intent(out) :: result !! Complete-enumeration tails, or a limit/error status without partial p-values.
      type(perm_control_type), intent(in), optional :: control !! Rounding and table caps; optional caller-supplied selection matrix.
      type(perm_control_type) :: options
      integer, allocatable :: labels(:), counts(:), cm(:, :)
      real(dp), allocatable :: statistics(:)
      integer :: status, i, n, m
      options = perm_control_type()
      if (present(control)) options = control
      result%status = control_status(options)
      if (result%status /= perm_ok) return
      result%status = data_status(scores)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      call label_groups(group, labels, counts, status)
      result%status = status
      if (status /= perm_ok) return
      result%status = perm_invalid
      if (size(counts) /= 2) return
      n = size(scores)
      m = counts(2)
      if (allocated(options%cm)) then
         if (size(options%cm, 2) /= n .or. size(options%cm, 1) < 1) return
         if (any(options%cm /= 0 .and. options%cm /= 1)) return
         if (any(sum(options%cm, dim=2) /= m)) return
         result%status = perm_limit
         if (size(options%cm, kind=int64) > options%max_cells) return
         if (int(size(options%cm, 1), int64) > options%max_permutations) return
         cm = options%cm
      else
         call choose_matrix(n, m, cm, status, options%max_permutations, options%max_cells)
         result%status = status
         if (status /= perm_ok) return
      end if
      result%status = perm_limit
      if (int(size(cm, 1), int64) * int(n, int64) > options%max_work) return
      allocate(statistics(size(cm, 1)), stat=status)
      if (status /= 0) then
         result%status = perm_allocation
         return
      end if
      do i = 1, size(statistics)
         statistics(i) = sum(scores, mask=cm(i, :) == 1)
      end do
      call calc_pvals(statistics, sum(scores, mask=labels == 2), result, digits=options%digits)
      result%method = 'exact.ce'
   end subroutine twosample_exact_ce

   pure subroutine twosample_exact_network(scores, group, result, control)
      real(dp), intent(in) :: scores(:) !! Finite arbitrary real scores; tied scores are grouped without lattice rounding.
      integer, intent(in) :: group(:) !! Exactly two labels; larger label defines the sampled group.
      type(perm_result), intent(out) :: result !! Hypergeometric-weighted exact tails, or a resource/error status.
      type(perm_control_type), intent(in), optional :: control !! Significant digits and enumeration, allocation and work caps.
      type(perm_control_type) :: options
      integer, allocatable :: labels(:), counts(:), multiplicity(:)
      integer(int64), allocatable :: paths(:), next_paths(:)
      real(dp), allocatable :: values(:), statistics(:), logweights(:)
      integer(int64) :: operations, cap
      real(dp) :: total, t0
      integer :: status, n, m, ng, i, j, q, position, remaining
      logical :: complement
      options = perm_control_type()
      if (present(control)) options = control
      result%status = control_status(options)
      if (result%status /= perm_ok) return
      result%status = data_status(scores)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      call label_groups(group, labels, counts, status)
      result%status = status
      if (status /= perm_ok) return
      result%status = perm_invalid
      if (size(counts) /= 2) return
      n = size(scores)
      m = minval(counts)
      complement = counts(2) > counts(1)
      t0 = sum(scores, mask=labels == 2)
      total = sum(scores)
      if (m == 1) then
         result%status = perm_limit
         if (int(n, int64) > options%max_permutations) return
         if (3_int64 * n > options%max_cells .or. int(n, int64) > options%max_work) return
         values = scores
         if (complement) values = total - scores
         call calc_pvals(values, t0, result, digits=options%digits)
         result%method = 'exact.network'
         return
      end if
      result%status = perm_limit
      if (int(2 * n + 2 * (m + 1), int64) > options%max_cells) return
      allocate(values(n), multiplicity(n), paths(0:m), next_paths(0:m), stat=status)
      if (status /= 0) then
         result%status = perm_allocation
         return
      end if
      values = scores
      call r_sort_values_in_place(values)
      ng = 1
      multiplicity = 0
      multiplicity(1) = 1
      do i = 2, n
         if (values(i) /= values(ng)) then
            ng = ng + 1
            values(ng) = values(i)
         end if
         multiplicity(ng) = multiplicity(ng) + 1
      end do
      result%status = perm_limit
      if (ng > 512) return
      if (int(2 * n + 2 * (m + 1), int64) > options%max_cells) return
      cap = options%max_permutations
      operations = 0_int64
      paths = 0_int64
      paths(0) = 1_int64
      remaining = n
      do i = 1, ng
         next_paths = 0_int64
         remaining = remaining - multiplicity(i)
         do j = 0, m
            if (paths(j) == 0_int64) cycle
            do q = max(0, m - remaining - j), min(multiplicity(i), m - j)
               operations = operations + 1_int64
               if (operations > options%max_work) return
               next_paths(j + q) = min(cap + 1_int64, next_paths(j + q) + paths(j))
            end do
         end do
         paths = next_paths
      end do
      if (paths(m) > cap) return
      if (paths(m) > (options%max_cells - 2_int64 * n - 2_int64 * (m + 1)) / 2_int64) return
      if (paths(m) * int(ng, int64) > options%max_work - operations) return
      allocate(statistics(int(paths(m))), logweights(int(paths(m))), stat=status)
      if (status /= 0) then
         result%status = perm_allocation
         return
      end if
      position = 0
      call visit(1, m, n, 0.0_dp, 0.0_dp, position, statistics, logweights)
      if (complement) statistics = total - statistics
      logweights = exp(logweights - maxval(logweights))
      call calc_pvals(statistics, t0, result, digits=options%digits, weights=logweights)
      result%method = 'exact.network'
   contains
      pure recursive subroutine visit(stage, needed, available, accumulated, logweight, cursor, stats, weights)
         integer, intent(in) :: stage !! Current distinct-score group, from one through ng.
         integer, intent(in) :: needed !! Number of observations still needed in the selected group.
         integer, intent(in) :: available !! Number of observations in this and later score groups.
         real(dp), intent(in) :: accumulated !! Current selected-score sum.
         real(dp), intent(in) :: logweight !! Logarithm of the current product of binomial multiplicities.
         integer, intent(inout) :: cursor !! Number of enumerated occupancy paths already written.
         real(dp), intent(inout) :: stats(:) !! Preallocated statistics, one entry per feasible occupancy path.
         real(dp), intent(inout) :: weights(:) !! Preallocated logarithmic path multiplicities.
         integer :: take, later
         real(dp) :: value, weight
         later = available - multiplicity(stage)
         do take = max(0, needed - later), min(multiplicity(stage), needed)
            value = accumulated + real(take, dp) * values(stage)
            weight = logweight + r_log_choose(multiplicity(stage), take)
            if (stage == ng) then
               cursor = cursor + 1
               stats(cursor) = value
               weights(cursor) = weight
            else
               call visit(stage + 1, needed - take, later, value, weight, cursor, stats, weights)
            end if
         end do
      end subroutine visit
   end subroutine twosample_exact_network

   pure subroutine random_index(rng, n, value)
      type(perm_rng), intent(inout) :: rng !! Park-Miller state in 1 through 2,147,483,646.
      integer, intent(in) :: n !! Uniform index upper bound, from one through 100,000.
      integer, intent(out) :: value !! Unbiased discrete uniform draw from 1 through n.
      integer(int64) :: limit, raw
      limit = ((rng_modulus - 1_int64) / int(n, int64)) * int(n, int64)
      do
         rng%state = modulo(48271_int64 * rng%state, rng_modulus)
         raw = rng%state - 1_int64
         if (raw < limit) exit
      end do
      value = int(modulo(raw, int(n, int64))) + 1
   end subroutine random_index

   pure subroutine monte_carlo(scores, covariate, labels, counts, mode, result, control, alternative, rng)
      real(dp), intent(in) :: scores(:) !! Validated finite scores.
      real(dp), intent(in) :: covariate(:) !! Numerical covariate for two-sample or trend statistics.
      integer, intent(in) :: labels(:) !! Dense categorical labels for the k-sample statistic.
      integer, intent(in) :: counts(:) !! Positive categorical group counts, unused for linear statistics.
      integer, intent(in) :: mode !! Statistic selector: 1 for a linear sum, 2 for k-sample between-group variation.
      type(perm_result), intent(out) :: result !! Fixed-count Monte Carlo result with add-one correction and p-value interval.
      type(perm_control_type), intent(in) :: control !! Validated numerical controls and explicit resource limits.
      character(len=*), intent(in) :: alternative !! Recognized tail name; ignored for a k-sample statistic.
      type(perm_rng), intent(inout), optional :: rng !! Optional stream, updated on success; required when set_seed is false.
      type(perm_rng) :: stream
      real(dp), allocatable :: statistics(:), centered(:), covariate_work(:)
      integer, allocatable :: permutation(:)
      real(dp) :: observed
      integer :: n, i, j, draw, swap, astat
      result%status = control_status(control)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (.not. valid_alternative(alternative)) return
      n = size(scores)
      result%status = perm_limit
      if (int(n, int64) * int(control%nmc, int64) > control%max_work) return
      if (int(control%nmc, int64) + 3_int64 * n > control%max_cells) return
      result%status = perm_invalid
      if (.not. control%set_seed .and. .not. present(rng)) return
      stream%state = 1_int64 + modulo(control%seed, rng_modulus - 1_int64)
      if (.not. control%set_seed) then
         stream = rng
         if (stream%state < 1_int64 .or. stream%state >= rng_modulus) return
      end if
      allocate(statistics(control%nmc), centered(n), covariate_work(n), permutation(n), stat=astat)
      if (astat /= 0) then
         result%status = perm_allocation
         return
      end if
      ! Power-of-two normalization is exact except when disparate input magnitudes underflow.
      ! It prevents products and quadratic k-sample statistics from spuriously becoming zero.
      centered = scores
      if (maxval(abs(scores)) > 0.0_dp) centered = scale(scores, -exponent(maxval(abs(scores))))
      covariate_work = covariate
      if (maxval(abs(covariate)) > 0.0_dp) then
         covariate_work = scale(covariate, -exponent(maxval(abs(covariate))))
      end if
      if (mode == 2) centered = centered - sum(centered) / real(n, dp)
      if (mode == 2) then
         observed = k_statistic(centered, labels, counts)
      else
         observed = sum(centered * covariate_work)
      end if
      do draw = 1, control%nmc
         permutation = [(i, i=1,n)]
         do i = n, 2, -1
            call random_index(stream, i, j)
            swap = permutation(i)
            permutation(i) = permutation(j)
            permutation(j) = swap
         end do
         if (mode == 2) then
            statistics(draw) = k_statistic(centered, labels(permutation), counts)
         else
            statistics(draw) = sum(centered * covariate_work(permutation))
         end if
      end do
      call calc_pvals_mc(statistics, observed, result, control%digits, alternative, mode == 2, control%confidence)
      result%method = 'exact.mc'
      if (result%status == perm_ok .and. present(rng)) rng = stream
   end subroutine monte_carlo

   pure subroutine twosample_exact_mc(scores, group, result, control, alternative, rng)
      real(dp), intent(in) :: scores(:) !! Finite scores, one per binary group label.
      integer, intent(in) :: group(:) !! Binary labels zero/one, matching the low-level R Monte Carlo routine.
      type(perm_result), intent(out) :: result !! Monte Carlo tails with an added observed permutation and binomial interval.
      type(perm_control_type), intent(in), optional :: control !! Sampling controls; low-level default is 9,999 random permutations.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), less, greater, or two.sidedAbs.
      type(perm_rng), intent(inout), optional :: rng !! Explicit local stream; needed for set_seed=false and updated on success.
      type(perm_control_type) :: options
      character(len=16) :: alt
      integer :: counts(2)
      options = perm_control_type(nmc=9999)
      if (present(control)) options = control
      result%status = data_status(scores)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      if (any(group /= 0 .and. group /= 1)) return
      counts = [count(group == 0), count(group == 1)]
      if (any(counts == 0)) return
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) then
            result%status = perm_invalid
            return
         end if
         alt = alternative
      end if
      call monte_carlo(scores, real(group, dp), group + 1, counts, 1, result, options, alt, rng)
   end subroutine twosample_exact_mc

   pure subroutine ksample_exact_mc(scores, group, result, control, rng)
      real(dp), intent(in) :: scores(:) !! Finite scores, one per categorical group label.
      integer, intent(in) :: group(:) !! Integer categorical labels with at least two distinct groups.
      type(perm_result), intent(out) :: result !! Upper-tail k-sample Monte Carlo p-value and Clopper-Pearson interval.
      type(perm_control_type), intent(in), optional :: control !! Sampling controls; low-level default is 9,999 random permutations.
      type(perm_rng), intent(inout), optional :: rng !! Explicit local stream; needed for set_seed=false and updated on success.
      type(perm_control_type) :: options
      integer, allocatable :: labels(:), counts(:)
      integer :: status
      options = perm_control_type(nmc=9999)
      if (present(control)) options = control
      result%status = data_status(scores)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      call label_groups(group, labels, counts, status)
      result%status = status
      if (status /= perm_ok) return
      result%status = perm_invalid
      if (size(counts) < 2) return
      call monte_carlo(scores, scores, labels, counts, 2, result, options, 'two.sided', rng)
   end subroutine ksample_exact_mc

   pure subroutine trend_exact_mc(scores, group, result, control, alternative, rng)
      real(dp), intent(in) :: scores(:) !! Finite response scores.
      real(dp), intent(in) :: group(:) !! Finite numerical trend covariate, of the same length as scores.
      type(perm_result), intent(out) :: result !! Linear-statistic Monte Carlo tails and uncertainty interval.
      type(perm_control_type), intent(in), optional :: control !! Sampling controls; low-level default is 999 permutations.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), less, greater, or two.sidedAbs.
      type(perm_rng), intent(inout), optional :: rng !! Explicit local stream; needed for set_seed=false and updated on success.
      type(perm_control_type) :: options
      character(len=16) :: alt
      options = perm_control_type()
      if (present(control)) options = control
      result%status = data_status(scores)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      if (size(group) /= size(scores)) return
      result%status = data_status(group)
      if (result%status /= perm_ok) return
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) then
            result%status = perm_invalid
            return
         end if
         alt = alternative
      end if
      call monte_carlo(scores, group, [1], [1], 1, result, options, alt, rng)
   end subroutine trend_exact_mc

   pure integer function unique_count(values) result(number)
      real(dp), intent(in) :: values(:) !! Finite values; an empty vector has zero distinct values.
      real(dp), allocatable :: sorted(:)
      integer :: i
      number = 0
      if (size(values) == 0) return
      sorted = values
      call r_sort_values_in_place(sorted)
      number = 1
      do i = 2, size(sorted)
         if (sorted(i) /= sorted(i - 1)) number = number + 1
      end do
   end function unique_count

   pure subroutine method_rule_ts1(scores, group, method, status, exact, bounds)
      real(dp), intent(in) :: scores(:) !! Finite scores used to count distinct observations within groups.
      integer, intent(in) :: group(:) !! Exactly two integer labels, one per score.
      character(len=*), intent(out) :: method !! exact.network, exact.mc, or pclt; provide at least 13 characters.
      integer, intent(out) :: status !! perm_ok on success or perm_invalid for malformed rule inputs.
      integer, intent(in), optional :: exact !! -1 automatic (default), 0 asymptotic, 1 exact or Monte Carlo.
      integer, intent(in), optional :: bounds(:) !! Positive sample-size bounds; default [1000,200,100,50,16].
      integer, allocatable :: labels(:), counts(:), thresholds(:)
      integer :: choice, kmin, index, threshold
      logical :: below
      method = ''
      status = data_status(scores)
      if (status /= perm_ok) return
      status = perm_invalid
      if (size(group) /= size(scores) .or. len(method) < 13) return
      call label_groups(group, labels, counts, status)
      if (status /= perm_ok) return
      status = perm_invalid
      if (size(counts) /= 2) return
      choice = -1
      if (present(exact)) choice = exact
      if (choice < -1 .or. choice > 1) return
      thresholds = [1000, 200, 100, 50, 16]
      if (present(bounds)) thresholds = bounds
      if (size(thresholds) < 1) return
      if (any(thresholds < 1)) return
      kmin = min(unique_count(pack(scores, labels == 1)), unique_count(pack(scores, labels == 2)))
      index = max(1, kmin - 1)
      threshold = minval(thresholds)
      if (index <= size(thresholds)) threshold = thresholds(index)
      below = size(scores) <= threshold
      method = 'pclt'
      if (choice /= 0) then
         if (below) then
            method = 'exact.network'
         else if (choice == 1) then
            method = 'exact.mc'
         end if
      end if
      status = perm_ok
   end subroutine method_rule_ts1

   pure subroutine method_rule_ks1(scores, group, method, status, exact, bound)
      real(dp), intent(in) :: scores(:) !! Finite scores, used to validate the group sample sizes.
      integer, intent(in) :: group(:) !! Integer categorical group labels matching scores.
      character(len=*), intent(out) :: method !! exact.mc or pclt; provide at least eight characters.
      integer, intent(out) :: status !! perm_ok on success or a validation status.
      integer, intent(in), optional :: exact !! -1 automatic (default), 0 asymptotic, 1 Monte Carlo.
      integer, intent(in), optional :: bound !! Automatic Monte Carlo threshold for the smallest group, default five.
      integer, allocatable :: labels(:), counts(:)
      integer :: choice, threshold
      method = ''
      status = data_status(scores)
      if (status /= perm_ok) return
      status = perm_invalid
      if (size(group) /= size(scores) .or. len(method) < 8) return
      call label_groups(group, labels, counts, status)
      if (status /= perm_ok) return
      status = perm_invalid
      if (size(counts) < 2) return
      choice = -1
      threshold = 5
      if (present(exact)) choice = exact
      if (present(bound)) threshold = bound
      if (choice < -1 .or. choice > 1 .or. threshold < 1) return
      method = 'pclt'
      if (choice == 1 .or. (choice == -1 .and. minval(counts) <= threshold)) method = 'exact.mc'
      status = perm_ok
   end subroutine method_rule_ks1

   pure subroutine method_rule_trend1(scores, covariate, method, status, exact, bound)
      real(dp), intent(in) :: scores(:) !! Finite response scores; sample size determines automatic selection.
      real(dp), intent(in) :: covariate(:) !! Finite numerical covariate, matching scores.
      character(len=*), intent(out) :: method !! exact.mc or pclt; provide at least eight characters.
      integer, intent(out) :: status !! perm_ok on success or a validation status.
      integer, intent(in), optional :: exact !! -1 automatic (default), 0 asymptotic, 1 Monte Carlo.
      integer, intent(in), optional :: bound !! Automatic Monte Carlo threshold for total sample size, default 20.
      integer :: choice, threshold
      method = ''
      status = data_status(scores)
      if (status /= perm_ok) return
      status = perm_invalid
      if (size(covariate) /= size(scores) .or. len(method) < 8) return
      status = data_status(covariate)
      if (status /= perm_ok) return
      status = perm_invalid
      choice = -1
      threshold = 20
      if (present(exact)) choice = exact
      if (present(bound)) threshold = bound
      if (choice < -1 .or. choice > 1 .or. threshold < 1) return
      method = 'pclt'
      if (choice == 1 .or. (choice == -1 .and. size(scores) <= threshold)) method = 'exact.mc'
      status = perm_ok
   end subroutine method_rule_trend1

   pure subroutine perm_ts(x, y, result, alternative, exact, method, control, rng)
      real(dp), intent(in) :: x(:) !! First nonempty finite sample; greater tests whether this sample has larger scores.
      real(dp), intent(in) :: y(:) !! Second nonempty finite sample.
      type(perm_result), intent(out) :: result !! Selected test, mean difference, tail probabilities and optional MC uncertainty.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), less, greater, or two.sidedAbs.
      integer, intent(in), optional :: exact !! -1 automatic (default), 0 asymptotic, 1 exact or Monte Carlo.
      character(len=*), intent(in), optional :: method !! Override: exact.network, exact.ce, exact.mc, or pclt.
      type(perm_control_type), intent(in), optional :: control !! Validated rounding, simulation and resource options.
      type(perm_rng), intent(inout), optional :: rng !! Explicit local stream for reproducible continuation between Monte Carlo calls.
      type(perm_control_type) :: options
      real(dp), allocatable :: scores(:)
      integer, allocatable :: group(:)
      character(len=16) :: selected, alt
      integer :: status, nx, ny
      options = perm_control_type()
      if (present(control)) options = control
      result%status = control_status(options)
      if (result%status /= perm_ok) return
      result%status = perm_invalid
      nx = size(x)
      ny = size(y)
      if (nx < 1 .or. ny < 1) return
      result%status = perm_limit
      if (int(nx, int64) + int(ny, int64) > max_observations) return
      scores = [x, y]
      group = [spread(1, 1, nx), spread(0, 1, ny)]
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) then
            result%status = perm_invalid
            return
         end if
         alt = alternative
      end if
      if (alt == 'two.sided' .and. options%tsmethod == 'abs') alt = 'two.sidedAbs'
      result%status = perm_invalid
      if (.not. valid_alternative(alt)) return
      if (present(method)) then
         if (len_trim(method) > len(selected)) then
            result%status = perm_invalid
            return
         end if
         selected = method
      else
         call method_rule_ts1(scores, group, selected, status, exact)
         result%status = status
         if (status /= perm_ok) return
      end if
      select case (selected)
      case ('exact.network')
         call twosample_exact_network(scores, group, result, options)
      case ('exact.ce')
         call twosample_exact_ce(scores, group, result, options)
      case ('exact.mc')
         call twosample_exact_mc(scores, group, result, options, alt, rng)
      case ('pclt')
         call twosample_pclt(scores, group, result)
      case default
         result%status = perm_invalid
         return
      end select
      if (result%status /= perm_ok) return
      result%estimate = sum(x / real(nx, dp)) - sum(y / real(ny, dp))
      call select_probability(result, alt)
   end subroutine perm_ts

   pure subroutine perm_ks(scores, group, result, exact, method, control, rng)
      real(dp), intent(in) :: scores(:) !! Finite response scores for a k-sample test, not a Kolmogorov-Smirnov test.
      integer, intent(in) :: group(:) !! Integer categorical labels, with at least two distinct groups.
      type(perm_result), intent(out) :: result !! K-sample p-value and asymptotic statistic or MC uncertainty.
      integer, intent(in), optional :: exact !! -1 automatic (default), 0 asymptotic, 1 Monte Carlo.
      character(len=*), intent(in), optional :: method !! Optional explicit pclt or exact.mc selection.
      type(perm_control_type), intent(in), optional :: control !! Rounding, fixed-count simulation and resource options.
      type(perm_rng), intent(inout), optional :: rng !! Explicit local stream for reproducible Monte Carlo continuation.
      type(perm_control_type) :: options
      character(len=16) :: selected
      integer :: status
      options = perm_control_type()
      if (present(control)) options = control
      result%status = control_status(options)
      if (result%status /= perm_ok) return
      if (present(method)) then
         if (len_trim(method) > len(selected)) then
            result%status = perm_invalid
            return
         end if
         selected = method
      else
         call method_rule_ks1(scores, group, selected, status, exact)
         result%status = status
         if (status /= perm_ok) return
      end if
      select case (selected)
      case ('pclt')
         call ksample_pclt(scores, group, result)
      case ('exact.mc')
         call ksample_exact_mc(scores, group, result, options, rng)
      case default
         result%status = perm_invalid
      end select
   end subroutine perm_ks

   pure subroutine perm_trend(scores, covariate, result, alternative, exact, method, control, rng)
      real(dp), intent(in) :: scores(:) !! Finite response scores.
      real(dp), intent(in) :: covariate(:) !! Finite numerical covariate, matching the scores.
      type(perm_result), intent(out) :: result !! Trend-test probability, correlation estimate and any MC uncertainty.
      character(len=*), intent(in), optional :: alternative !! two.sided (default), less, greater, or two.sidedAbs.
      integer, intent(in), optional :: exact !! -1 automatic (default), 0 asymptotic, 1 Monte Carlo.
      character(len=*), intent(in), optional :: method !! Optional explicit pclt or exact.mc selection.
      type(perm_control_type), intent(in), optional :: control !! Rounding, fixed-count simulation and resource options.
      type(perm_rng), intent(inout), optional :: rng !! Explicit local stream for reproducible Monte Carlo continuation.
      type(perm_control_type) :: options
      type(perm_result) :: asymptotic
      character(len=16) :: selected, alt
      integer :: status
      options = perm_control_type()
      if (present(control)) options = control
      result%status = control_status(options)
      if (result%status /= perm_ok) return
      alt = 'two.sided'
      if (present(alternative)) then
         if (len_trim(alternative) > len(alt)) then
            result%status = perm_invalid
            return
         end if
         alt = alternative
      end if
      if (alt == 'two.sided' .and. options%tsmethod == 'abs') alt = 'two.sidedAbs'
      result%status = perm_invalid
      if (.not. valid_alternative(alt)) return
      if (present(method)) then
         if (len_trim(method) > len(selected)) then
            result%status = perm_invalid
            return
         end if
         selected = method
      else
         call method_rule_trend1(scores, covariate, selected, status, exact)
         result%status = status
         if (status /= perm_ok) return
      end if
      select case (selected)
      case ('pclt')
         call trend_pclt(scores, covariate, result)
      case ('exact.mc')
         call trend_exact_mc(scores, covariate, result, options, alt, rng)
      case default
         result%status = perm_invalid
         return
      end select
      if (result%status /= perm_ok) return
      call trend_pclt(scores, covariate, asymptotic)
      result%estimate = ieee_value(0.0_dp, ieee_quiet_nan)
      if (asymptotic%status == perm_ok) result%estimate = asymptotic%estimate
      call select_probability(result, alt)
   end subroutine perm_trend
end module perm_api
