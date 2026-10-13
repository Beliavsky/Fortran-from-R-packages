! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from monobin 0.2.4, Andrija Djurovic (2022).
module monobin
   use r_kinds, only : dp
   use hmisc, only : cut2_explicit
   use r_ordering, only : r_order
   use r_sorting, only : r_average_ranks
   use r_distributions, only : r_pt, r_pnorm
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   public :: dp, bin_options, bin_summary, bin_result, fit_bins, apply_bins
   public :: pct_bin, iso_bin, woe_bin, cum_bin, sts_bin, ndr_bin, mdt_bin
   public :: bin_ok, bin_invalid, bin_no_data, bin_constant_x, bin_constant_y
   integer, parameter :: bin_ok = 0, bin_invalid = 1, bin_no_data = 2
   integer, parameter :: bin_constant_x = 3, bin_constant_y = 4

   type :: bin_options
      !! Controls shared across algorithms; g = 0 selects the algorithm default.
      integer :: g = 0
      integer :: direction = 0 !! 0 = automatic, 1 = increasing, -1 = decreasing.
      integer :: target = 0 !! 0 = automatic, 1 = binary, 2 = continuous.
      logical :: separate_special = .false.
      logical :: woe_trend = .true.
      real(dp) :: min_pct_obs = 0.05_dp
      real(dp) :: min_avg_rate = 0.01_dp
      real(dp) :: woe_gap = 0.1_dp
      real(dp) :: p_value = 0.05_dp
   end type

   type :: bin_summary
      !! Sufficient statistics and WoE/IV for one complete or special bin.
      integer :: count = 0
      logical :: special = .false.
      real(dp) :: special_value = 0.0_dp
      real(dp) :: y_sum = 0.0_dp
      real(dp) :: y_mean = 0.0_dp
      real(dp) :: y_sumsq = 0.0_dp
      real(dp) :: y_m2 = 0.0_dp !! Centered response sum of squares for stable inference.
      real(dp) :: x_mean = 0.0_dp
      real(dp) :: x_min = 0.0_dp
      real(dp) :: x_max = 0.0_dp
      real(dp) :: dist_good = 0.0_dp
      real(dp) :: dist_bad = 0.0_dp
      real(dp) :: woe = 0.0_dp
      real(dp) :: iv = 0.0_dp
      real(dp) :: p_value = 0.0_dp
   end type

   type :: bin_result
      integer :: status = bin_invalid
      integer :: direction = 0
      integer :: n_special = 0
      logical :: binary = .false.
      logical :: separate_special = .false.
      type(bin_summary), allocatable :: bins(:)
      integer, allocatable :: membership(:) !! Original input length; 0 means missing target.
      integer, allocatable :: retained(:) !! Original indices retained after missing-target removal.
      real(dp), allocatable :: special_values(:) !! User-defined finite special values.
      real(dp), allocatable :: breaks(:) !! Lower boundaries of complete bins 2 through K.
   end type
contains
   pure subroutine pct_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'pct', result, options, special)
   end subroutine pct_bin

   pure subroutine iso_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'iso', result, options, special)
   end subroutine iso_bin

   pure subroutine woe_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'woe', result, options, special)
   end subroutine woe_bin

   pure subroutine cum_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'cum', result, options, special)
   end subroutine cum_bin

   pure subroutine sts_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'sts', result, options, special)
   end subroutine sts_bin

   pure subroutine ndr_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'ndr', result, options, special)
   end subroutine ndr_bin

   pure subroutine mdt_bin(x, y, result, options, special)
      real(dp), intent(in) :: x(:) !! Numeric predictor; nonfinite values are always special cases.
      real(dp), intent(in) :: y(:) !! Same-length response; NaNs are omitted and infinities rejected.
      type(bin_result), intent(out) :: result !! Fitted summaries, cutpoints, membership, and status.
      type(bin_options), intent(in), optional :: options !! Algorithm controls; omitted fields use upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Additional predictor values to keep out of complete bins.
      call fit_bins(x, y, 'mdt', result, options, special)
   end subroutine mdt_bin

   pure subroutine fit_bins(x, y, method, result, options, special)
      real(dp), intent(in) :: x(:) !! Predictor observations; NaN and infinities are special cases.
      real(dp), intent(in) :: y(:) !! Responses aligned with x; missing responses are discarded.
      character(len=*), intent(in) :: method !! One of pct, iso, woe, cum, sts, ndr, or mdt.
      type(bin_result), intent(out) :: result !! Fit output; a nonzero status leaves no fitted bins.
      type(bin_options), intent(in), optional :: options !! Shared algorithm controls with upstream defaults.
      real(dp), intent(in), optional :: special(:) !! Extra special values in addition to nonfinite predictors.
      type(bin_options) :: opt
      type(bin_summary), allocatable :: cc(:), sc(:), tmp(:)
      real(dp), allocatable :: xc(:), yc(:), ranks_x(:), ranks_y(:), sv(:)
      integer, allocatable :: keep(:), ci(:), order(:), lab(:)
      real(dp) :: total_y, min_rate, cov
      integer :: n, nc, ns, i, j, k, g, d, min_obs, nvalid
      logical :: binary
      opt = bin_options()
      if (present(options)) opt = options
      allocate(result%bins(0), result%membership(size(x)), result%retained(0), &
         result%breaks(0), result%special_values(0))
      result%membership = 0
      n = size(x)
      if (n /= size(y) .or. n == 0) return
      if (abs(opt%direction) > 1 .or. opt%target < 0 .or. opt%target > 2) return
      if (opt%g < 0) return
      if (.not. all(ieee_is_finite([opt%min_pct_obs, opt%min_avg_rate, opt%woe_gap, opt%p_value]))) return
      if (opt%min_pct_obs < 0.0_dp .or. opt%min_pct_obs > 1.0_dp) return
      if (opt%min_avg_rate < 0.0_dp .or. opt%woe_gap < 0.0_dp) return
      if (opt%p_value <= 0.0_dp .or. opt%p_value > 1.0_dp) return
      select case (method)
      case ('pct', 'iso', 'woe', 'cum', 'sts', 'ndr', 'mdt')
      case default
         return
      end select
      if (present(special)) result%special_values = pack(special, ieee_is_finite(special))
      allocate(keep(n), ci(n), sv(n), sc(n))
      nc = 0
      ns = 0
      nvalid = 0
      total_y = 0.0_dp
      binary = .true.
      do i = 1, n
         if (ieee_is_nan(y(i))) cycle
         if (.not. ieee_is_finite(y(i))) return
         nvalid = nvalid + 1
         keep(nvalid) = i
         total_y = total_y + y(i)
         binary = binary .and. (y(i) == 0.0_dp .or. y(i) == 1.0_dp)
         if (is_special(x(i), result%special_values)) then
            j = 0
            if (opt%separate_special) then
               do k = 1, ns
                  if (same_value(x(i), sv(k))) then
                     j = k
                     exit
                  end if
               end do
            else if (ns > 0) then
               j = 1
            end if
            if (j == 0) then
               ns = ns + 1
               j = ns
               sv(j) = x(i)
               sc(j)%special = .true.
               sc(j)%special_value = x(i)
            end if
            call add_observation(sc(j), x(i), y(i))
            result%membership(i) = j
         else
            nc = nc + 1
            ci(nc) = i
         end if
      end do
      if (opt%target == 1 .and. .not. binary) return
      if (opt%target == 2) binary = .false.
      if (nc == 0) then
         result%status = bin_no_data
         return
      end if
      xc = x(ci(:nc))
      yc = y(ci(:nc))
      if (all(yc == yc(1))) then
         result%status = bin_constant_y
         return
      end if
      if (all(xc == xc(1))) then
         result%status = bin_constant_x
         return
      end if
      call r_order(xc, order)
      xc = xc(order)
      yc = yc(order)
      ci(:nc) = ci(order)
      d = opt%direction
      if (d == 0) then
         call r_average_ranks(xc, ranks_x)
         call r_average_ranks(yc, ranks_y)
         cov = sum((ranks_x - sum(ranks_x) / real(nc, dp)) * &
            (ranks_y - sum(ranks_y) / real(nc, dp)))
         if (cov > 0.0_dp) d = 1
         if (cov < 0.0_dp) d = -1
      end if
      g = opt%g
      if (g == 0) then
         g = 15
         if (method == 'mdt') g = 50
      end if
      min_obs = ceiling(max(30.0_dp, real(nvalid, dp) * opt%min_pct_obs))
      min_rate = opt%min_avg_rate
      if (binary) min_rate = real(ceiling(max(1.0_dp, total_y * opt%min_avg_rate)), dp)
      select case (method)
      case ('pct')
         do
            call percentile_groups(xc, yc, g, cc)
            if (size(cc) <= 1) exit
            if (.not. binary .and. opt%woe_trend) then
               tmp = [sc(:ns), cc]
               call evidence(tmp, binary)
               if (strict_monotone(tmp(ns + 1:)%woe, opt%direction)) exit
            else
               if (strict_monotone(cc%y_mean, opt%direction)) exit
            end if
            g = g - 1
            if (g < 1) exit
         end do
      case ('cum')
         call percentile_groups(xc, yc, g, cc)
         call cumulative_merge(cc, d)
      case ('mdt')
         call tree_bins(xc, yc, g, min_obs, min_rate, binary, d, cc)
      case default
         call isotonic_groups(xc, yc, d, cc)
         if (d == -1) cc = cc(size(cc):1:-1)
         call minimum_correction(cc, min_obs, min_rate, binary)
         if (d == -1) cc = cc(size(cc):1:-1)
         select case (method)
         case ('woe')
            call gap_merge(cc, sc(:ns), binary, opt%woe_gap)
         case ('sts')
            call significance_merge(cc, binary, opt%p_value, .false.)
         case ('ndr')
            call significance_merge(cc, binary, opt%p_value, .true.)
            call nested_reassignment(xc, yc, cc)
         end select
      end select
      result%bins = [sc(:ns), cc]
      if (method == 'sts') result%bins(:ns)%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      call evidence(result%bins, binary)
      result%status = bin_ok
      result%binary = binary
      result%direction = d
      result%n_special = ns
      result%separate_special = opt%separate_special
      result%retained = keep(:nvalid)
      if (size(cc) > 1) result%breaks = cc(2:)%x_min
      call apply_bins(result, x, lab)
      result%membership = lab
      do i = 1, n
         if (ieee_is_nan(y(i))) result%membership(i) = 0
      end do
   end subroutine fit_bins

   pure subroutine apply_bins(model, x, membership, response, woe)
      type(bin_result), intent(in) :: model !! Successfully fitted model with ordered complete-bin boundaries.
      real(dp), intent(in) :: x(:) !! New predictor values; ordinary tails extend to infinity.
      integer, allocatable, intent(out) :: membership(:) !! One-based bin indices; 0 for unseen special categories.
      real(dp), allocatable, intent(out), optional :: response(:) !! Fitted bin response means, NaN for unmapped cases.
      real(dp), allocatable, intent(out), optional :: woe(:) !! Fitted bin weights of evidence, NaN for unmapped cases.
      integer :: i, j, k
      allocate(membership(size(x)))
      membership = 0
      if (present(response)) then
         allocate(response(size(x)))
         response = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      if (present(woe)) then
         allocate(woe(size(x)))
         woe = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
      if (model%status /= bin_ok) return
      do i = 1, size(x)
         j = 0
         if (is_special(x(i), model%special_values)) then
            if (.not. model%separate_special .and. model%n_special > 0) then
               j = 1
            else
               do k = 1, model%n_special
                  if (same_value(x(i), model%bins(k)%special_value)) then
                     j = k
                     exit
                  end if
               end do
            end if
         else
            j = model%n_special + 1 + count(x(i) >= model%breaks)
         end if
         membership(i) = j
         if (j > 0) then
            if (present(response)) response(i) = model%bins(j)%y_mean
            if (present(woe)) woe(i) = model%bins(j)%woe
         end if
      end do
   end subroutine apply_bins

   pure elemental logical function same_value(a, b) result(equal)
      real(dp), intent(in) :: a !! First predictor value, possibly a NaN or infinity.
      real(dp), intent(in) :: b !! Second predictor value; all NaNs share one missing category.
      equal = .false.
      if (ieee_is_nan(a) .and. ieee_is_nan(b)) then
         equal = .true.
      else if (.not. ieee_is_nan(a) .and. .not. ieee_is_nan(b)) then
         equal = a == b
      end if
   end function same_value

   pure logical function is_special(x, special) result(found)
      real(dp), intent(in) :: x !! Predictor value to classify, including IEEE nonfinite values.
      real(dp), intent(in) :: special(:) !! Additional finite special-case codes.
      found = .not. ieee_is_finite(x)
      if (ieee_is_finite(x)) found = any(x == special)
   end function is_special

   pure subroutine add_observation(bin, x, y)
      type(bin_summary), intent(inout) :: bin !! Accumulated summary updated by one observation.
      real(dp), intent(in) :: x !! Predictor value; special bins may contain nonfinite values.
      real(dp), intent(in) :: y !! Finite response belonging to this bin.
      real(dp) :: w, oldmean
      oldmean = bin%y_mean
      if (bin%count == 0) then
         bin%x_min = x
         bin%x_max = x
         bin%x_mean = x
      else
         if (ieee_is_nan(bin%x_min) .or. ieee_is_nan(x)) then
            bin%x_min = ieee_value(0.0_dp, ieee_quiet_nan)
            bin%x_max = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            bin%x_min = min(bin%x_min, x)
            bin%x_max = max(bin%x_max, x)
         end if
         w = 1.0_dp / real(bin%count + 1, dp)
         if (ieee_is_finite(bin%x_mean) .and. ieee_is_finite(x)) then
            bin%x_mean = (1.0_dp - w) * bin%x_mean + w * x
         else if (ieee_is_nan(bin%x_mean) .or. ieee_is_nan(x)) then
            bin%x_mean = ieee_value(0.0_dp, ieee_quiet_nan)
         else if (.not. ieee_is_finite(bin%x_mean) .and. .not. ieee_is_finite(x)) then
            if (bin%x_mean /= x) bin%x_mean = ieee_value(0.0_dp, ieee_quiet_nan)
         else if (.not. ieee_is_finite(x)) then
            bin%x_mean = x
         end if
      end if
      bin%count = bin%count + 1
      bin%y_sum = bin%y_sum + y
      bin%y_sumsq = bin%y_sumsq + y * y
      bin%y_mean = bin%y_sum / real(bin%count, dp)
      bin%y_m2 = bin%y_m2 + (y - oldmean) * (y - bin%y_mean)
   end subroutine add_observation

   pure function combined(a, b) result(c)
      type(bin_summary), intent(in) :: a !! Left adjacent complete-bin summary.
      type(bin_summary), intent(in) :: b !! Right adjacent complete-bin summary.
      type(bin_summary) :: c
      c%count = a%count + b%count
      c%y_sum = a%y_sum + b%y_sum
      c%y_sumsq = a%y_sumsq + b%y_sumsq
      c%y_m2 = a%y_m2 + b%y_m2 + (a%y_mean - b%y_mean) ** 2 * &
         real(a%count, dp) * real(b%count, dp) / real(c%count, dp)
      c%y_mean = c%y_sum / real(c%count, dp)
      c%x_mean = (a%x_mean * real(a%count, dp) + b%x_mean * real(b%count, dp)) / real(c%count, dp)
      c%x_min = min(a%x_min, b%x_min)
      c%x_max = max(a%x_max, b%x_max)
   end function combined

   pure subroutine merge_at(bins, j)
      type(bin_summary), allocatable, intent(inout) :: bins(:) !! Ordered complete bins, shortened by one.
      integer, intent(in) :: j !! Left member of the adjacent pair; 1 <= j < size(bins).
      bins(j) = combined(bins(j), bins(j + 1))
      bins = [bins(:j), bins(j + 2:)]
   end subroutine merge_at

   pure subroutine isotonic_groups(x, y, direction, bins)
      real(dp), intent(in) :: x(:) !! Finite predictor sorted in nondecreasing order.
      real(dp), intent(in) :: y(:) !! Finite aligned response observations.
      integer, intent(in) :: direction !! Increasing 1, decreasing -1, or zero-correlation collapse 0.
      type(bin_summary), allocatable, intent(out) :: bins(:) !! Pooled-adjacent-violator summaries in x order.
      type(bin_summary) :: stack(size(x))
      integer :: i, k
      k = 0
      do i = 1, size(x)
         if (i > 1) then
            if (x(i) == x(i - 1)) then
               call add_observation(stack(k), x(i), y(i))
               cycle
            end if
         end if
         k = k + 1
         call add_observation(stack(k), x(i), y(i))
      end do
      bins = stack(:k)
      i = 1
      do while (i < size(bins))
         if (real(direction, dp) * (bins(i + 1)%y_mean - bins(i)%y_mean) <= 0.0_dp) then
            call merge_at(bins, i)
            i = max(1, i - 1)
         else
            i = i + 1
         end if
      end do
   end subroutine isotonic_groups

   pure subroutine minimum_correction(bins, min_obs, min_rate, binary)
      type(bin_summary), allocatable, intent(inout) :: bins(:) !! Isotonic bins to merge while retaining their x order.
      integer, intent(in) :: min_obs !! Minimum count, including the upstream floor of 30.
      real(dp), intent(in) :: min_rate !! Minimum bad count for binary or response mean for continuous targets.
      logical, intent(in) :: binary !! Whether the target is binary rather than continuous.
      integer :: i, stage, j
      logical :: bad
      do stage = 1, 2
         do
            j = 0
            do i = 1, size(bins)
               if (stage == 1) then
                  bad = bins(i)%count < min_obs
               else if (binary) then
                  bad = bins(i)%y_sum < min_rate
               else
                  bad = bins(i)%y_mean < min_rate
               end if
               if (bad) then
                  j = i
                  exit
               end if
            end do
            if (j == 0 .or. size(bins) == 1) exit
            call merge_at(bins, min(j, size(bins) - 1))
         end do
      end do
   end subroutine minimum_correction

   pure logical function strict_monotone(values, direction) result(ok)
      real(dp), intent(in) :: values(:) !! Ordered bin response means or evidence values.
      integer, intent(in) :: direction !! 0 accepts either strict direction; +/-1 forces the stated sign.
      logical :: increasing, decreasing
      ok = .false.
      if (any(ieee_is_nan(values))) return
      increasing = all(values(2:) > values(:size(values) - 1))
      decreasing = all(values(2:) < values(:size(values) - 1))
      ok = (direction >= 0 .and. increasing) .or. (direction <= 0 .and. decreasing)
   end function strict_monotone

   pure subroutine evidence(bins, binary)
      type(bin_summary), intent(inout) :: bins(:) !! Complete and special summaries receiving normalized WoE and IV.
      logical, intent(in) :: binary !! Binary uses good/bad ratios; continuous uses response/count shares.
      real(dp) :: n, s, a, b
      integer :: i
      n = real(sum(bins%count), dp)
      s = sum(bins%y_sum)
      do i = 1, size(bins)
         if (binary) then
            a = safe_ratio(real(bins(i)%count, dp) - bins(i)%y_sum, n - s)
            b = safe_ratio(bins(i)%y_sum, s)
         else
            a = safe_ratio(bins(i)%y_sum, s)
            b = real(bins(i)%count, dp) / n
         end if
         bins(i)%dist_good = a
         bins(i)%dist_bad = b
         bins(i)%woe = log_ratio(a, b)
         bins(i)%iv = (a - b) * bins(i)%woe
      end do
   end subroutine evidence

   pure elemental real(dp) function safe_ratio(a, b) result(value)
      real(dp), intent(in) :: a !! Numerator of a possibly undefined distribution share.
      real(dp), intent(in) :: b !! Denominator; zero returns the corresponding IEEE result without trapping.
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (b /= 0.0_dp) then
         value = a / b
      else if (a > 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else if (a < 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      end if
   end function safe_ratio

   pure elemental real(dp) function log_ratio(a, b) result(value)
      real(dp), intent(in) :: a !! First distribution share, nonnegative for an ordinary WoE statistic.
      real(dp), intent(in) :: b !! Second distribution share; zero cells retain infinite WoE.
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      if (a > 0.0_dp .and. b > 0.0_dp) then
         value = log(a / b)
      else if (a == 0.0_dp .and. b > 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else if (a > 0.0_dp .and. b == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      end if
   end function log_ratio

   pure subroutine gap_merge(bins, special, binary, gap)
      type(bin_summary), allocatable, intent(inout) :: bins(:) !! Corrected isotonic complete bins, merged in place.
      type(bin_summary), intent(in) :: special(:) !! Special bins included in the evidence normalization.
      logical, intent(in) :: binary !! Chooses binary or continuous evidence formula.
      real(dp), intent(in) :: gap !! Required strict minimum separation of adjacent WoE values.
      type(bin_summary), allocatable :: allbins(:)
      real(dp), allocatable :: differences(:)
      integer :: j, n, ns
      logical :: reverse
      ns = size(special)
      reverse = bins(1)%y_mean > bins(size(bins))%y_mean
      if (reverse) bins = bins(size(bins):1:-1)
      do while (size(bins) > 1)
         allbins = [special, bins]
         call evidence(allbins, binary)
         n = size(bins)
         differences = abs(allbins(ns + 2:)%woe - allbins(ns + 1:ns + n - 1)%woe)
         if (all(differences > gap)) exit
         j = minloc(differences, dim = 1)
         call merge_at(bins, j)
      end do
      if (reverse) bins = bins(size(bins):1:-1)
   end subroutine gap_merge
   pure subroutine percentile_groups(x, y, g, bins)
      !! Narrow adaptation of Hmisc::cut2(g): upper ECDF order-statistic endpoints,
      !! whole tied runs, and repeated cut suppression; not interpolated quantiles.
      real(dp), intent(in) :: x(:) !! Sorted finite predictors, retaining repeated values.
      real(dp), intent(in) :: y(:) !! Responses in the same observation order as x.
      integer, intent(in) :: g !! Requested initial equal-count groups; must be positive.
      type(bin_summary), allocatable, intent(out) :: bins(:) !! Nonempty groups at distinct cut2 upper endpoints.
      type(bin_summary) :: work(size(x))
      integer :: j, first, last, previous, k, i, n
      n = size(x)
      first = 1
      previous = 0
      k = 0
      do j = 1, g
         last = max(1, min(n, ceiling(real(j, dp) * (real(n, dp) / real(g, dp)))))
         do while (last < n)
            if (x(last + 1) /= x(last)) exit
            last = last + 1
         end do
         if (last == previous) cycle
         k = k + 1
         do i = first, last
            call add_observation(work(k), x(i), y(i))
         end do
         previous = last
         first = last + 1
      end do
      bins = work(:k)
   end subroutine percentile_groups

   pure subroutine cumulative_merge(bins, direction)
      type(bin_summary), allocatable, intent(inout) :: bins(:) !! Percentile groups replaced by cumulative-mean pools.
      integer, intent(in) :: direction !! Positive means increasing x trend; zero follows upstream decreasing choice.
      type(bin_summary) :: work(size(bins)), pool
      integer :: first, j, best, k, n, total
      real(dp) :: sy, rate, highest
      logical :: reverse
      reverse = direction == 1
      if (reverse) bins = bins(size(bins):1:-1)
      first = 1
      k = 0
      n = size(bins)
      do while (first <= n)
         sy = 0.0_dp
         total = 0
         highest = -huge(1.0_dp)
         best = first
         do j = first, n
            total = total + bins(j)%count
            sy = sy + bins(j)%y_sum
            rate = sy / real(total, dp)
            if (rate > highest) then
               highest = rate
               best = j
            end if
         end do
         pool = bins(first)
         do j = first + 1, best
            pool = combined(pool, bins(j))
         end do
         k = k + 1
         work(k) = pool
         first = best + 1
      end do
      bins = work(:k)
      if (reverse) bins = bins(size(bins):1:-1)
   end subroutine cumulative_merge

   pure subroutine nested_reassignment(x, y, bins)
      !! Preserve Hmisc explicit-cut endpoint semantics after stepwise regression.
      real(dp), intent(in) :: x(:) !! Ordered complete-case predictor observations.
      real(dp), intent(in) :: y(:) !! Responses aligned with x after dropping missing targets.
      type(bin_summary), allocatable, intent(inout) :: bins(:) !! Selected regression groups, rebuilt using cut2 intervals.
      integer, allocatable :: groups(:)
      real(dp), allocatable :: cuts(:)
      type(bin_summary) :: work(size(bins))
      integer :: status, i
      if (size(bins) < 2) return
      call cut2_explicit(x, [x(1), bins(2:)%x_min], .true., groups, cuts, status)
      if (status /= 0) return
      do i = 1, size(x)
         call add_observation(work(groups(i)), x(i), y(i))
      end do
      bins = work(:maxval(groups))
   end subroutine nested_reassignment

   pure subroutine significance_merge(bins, binary, significance, nested)
      type(bin_summary), allocatable, intent(inout) :: bins(:) !! Isotonic bins merged by weakest adjacent contrast.
      logical, intent(in) :: binary !! Binary responses use proportions or logistic nested-dummy Wald tests.
      real(dp), intent(in) :: significance !! Retain only contrasts with p < significance, in (0,1].
      logical, intent(in) :: nested !! True selects nested-dummy regression; false selects two-sample tests.
      real(dp), allocatable :: p(:), eta(:), weights(:)
      real(dp) :: n1, n2, v1, v2, se2, df, z, pooled, mse, residual
      integer :: j, k, n, nobs
      if (size(bins) == 1) then
         bins(1)%p_value = 99.0_dp
         return
      end if
      do while (size(bins) > 1)
         n = size(bins)
         nobs = sum(bins%count)
         allocate(p(n - 1))
         mse = 0.0_dp
         if (nested .and. .not. binary) then
            residual = sum(bins%y_m2)
            mse = max(0.0_dp, residual) / real(max(1, nobs - n), dp)
         else if (nested) then
            call logistic_groups(bins, eta, weights)
         end if
         do j = 1, n - 1
            n1 = real(bins(j)%count, dp)
            n2 = real(bins(j + 1)%count, dp)
            if (nested .and. binary) then
               se2 = 1.0_dp / weights(j) + 1.0_dp / weights(j + 1)
               z = abs(eta(j + 1) - eta(j)) / sqrt(se2)
               p(j) = 2.0_dp * r_pnorm(-z)
            else if (nested) then
               se2 = mse * (1.0_dp / n1 + 1.0_dp / n2)
               z = standardized_difference(bins(j)%y_mean, bins(j + 1)%y_mean, se2)
               p(j) = 2.0_dp * r_pt(-z, real(max(1, nobs - n), dp))
            else if (binary) then
               pooled = (bins(j)%y_sum + bins(j + 1)%y_sum) / (n1 + n2)
               se2 = pooled * (1.0_dp - pooled) * (1.0_dp / n1 + 1.0_dp / n2)
               z = standardized_difference(bins(j)%y_mean, bins(j + 1)%y_mean, se2)
               p(j) = r_pnorm(-z)
            else
               v1 = max(0.0_dp, bins(j)%y_m2) / (n1 - 1.0_dp)
               v2 = max(0.0_dp, bins(j + 1)%y_m2) / (n2 - 1.0_dp)
               se2 = v1 / n1 + v2 / n2
               z = standardized_difference(bins(j)%y_mean, bins(j + 1)%y_mean, se2)
               df = safe_ratio(se2 ** 2, (v1 / n1) ** 2 / (n1 - 1.0_dp) + (v2 / n2) ** 2 / (n2 - 1.0_dp))
               if (se2 == 0.0_dp) then
                  p(j) = r_pnorm(-z)
               else
                  p(j) = r_pt(-z, df)
               end if
            end if
         end do
         bins(1)%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
         bins(2:)%p_value = p
         if (all(p < significance)) exit
         k = maxloc(p, dim = 1)
         call merge_at(bins, k)
         deallocate(p)
      end do
      if (size(bins) == 1) bins(1)%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
   end subroutine significance_merge

   pure elemental real(dp) function standardized_difference(a, b, variance) result(z)
      real(dp), intent(in) :: a !! First finite response mean or model coefficient.
      real(dp), intent(in) :: b !! Second finite response mean or model coefficient.
      real(dp), intent(in) :: variance !! Nonnegative sampling variance of their difference.
      z = 0.0_dp
      if (variance > 0.0_dp) then
         z = abs(a - b) / sqrt(variance)
      else if (a /= b) then
         z = ieee_value(0.0_dp, ieee_positive_inf)
      end if
   end function standardized_difference

   pure subroutine logistic_groups(bins, eta, weights)
      !! Groupwise form of saturated nested-dummy binomial IRLS; no dense solve
      !! is needed because each group intercept is independently identifiable.
      type(bin_summary), intent(in) :: bins(:) !! Current groups, with binary event counts and denominators.
      real(dp), allocatable, intent(out) :: eta(:) !! Fitted group log odds from R-compatible 25-step IRLS.
      real(dp), allocatable, intent(out) :: weights(:) !! Last working information n*p*(1-p) used by glm covariance.
      real(dp) :: mu(size(bins)), neweta(size(bins)), dev, olddev, p, eps
      integer :: iter, i
      allocate(eta(size(bins)), weights(size(bins)))
      eta = 0.0_dp
      weights = 0.1875_dp * real(bins%count, dp)
      neweta = (log(3.0_dp) + 4.0_dp / 3.0_dp) * (2.0_dp * bins%y_mean - 1.0_dp)
      olddev = 2.0_dp * real(sum(bins%count), dp) * log(4.0_dp / 3.0_dp)
      eps = epsilon(1.0_dp)
      do iter = 1, 25
         eta = neweta
         mu = max(eps, min(1.0_dp - eps, 1.0_dp / (1.0_dp + exp(-eta))))
         dev = 0.0_dp
         do i = 1, size(bins)
            p = mu(i)
            dev = dev - 2.0_dp * (bins(i)%y_sum * log(p) + &
               (real(bins(i)%count, dp) - bins(i)%y_sum) * log(1.0_dp - p))
         end do
         if (abs(dev - olddev) / (0.1_dp + abs(dev)) < 1.0e-8_dp) exit
         olddev = dev
         weights = real(bins%count, dp) * mu * (1.0_dp - mu)
         neweta = eta + (bins%y_mean - mu) / (mu * (1.0_dp - mu))
      end do
   end subroutine logistic_groups

   pure subroutine tree_bins(x, y, g, min_obs, min_rate, binary, direction, bins)
      real(dp), intent(in) :: x(:) !! Ordered finite predictors used by the breadth-first decision tree.
      real(dp), intent(in) :: y(:) !! Aligned finite binary or continuous responses.
      integer, intent(in) :: g !! Number of cut2 candidate prebins recalculated inside every node.
      integer, intent(in) :: min_obs !! Required child count, with upstream floor of 30.
      real(dp), intent(in) :: min_rate !! Required bad count or continuous response mean.
      logical, intent(in) :: binary !! Select Gini gain for binary and residual SSE for continuous responses.
      integer, intent(in) :: direction !! Positive requests increasing, all other values decreasing.
      type(bin_summary), allocatable, intent(out) :: bins(:) !! Final ordered leaf-cut summaries.
      integer :: lower(2 * size(x)), upper(2 * size(x)), side(2 * size(x)), state(2 * size(x))
      integer :: snapshot(2 * size(x)), nnode, active, i, j, node, split, ncut, k
      real(dp) :: threshold(2 * size(x)), nodecut(2 * size(x)), cuts(2 * size(x)), rightmean
      type(bin_summary) :: work(size(x))
      integer, allocatable :: cutorder(:)
      nnode = 1
      lower(1) = 1
      upper(1) = size(x)
      side(1) = 0
      state(1) = 0
      nodecut(1) = ieee_value(0.0_dp, ieee_quiet_nan)
      threshold(1) = huge(1.0_dp)
      if (direction == 1) threshold(1) = -huge(1.0_dp)
      do
         active = 0
         do i = 1, nnode
            if (state(i) == 0) then
               active = active + 1
               snapshot(active) = i
            end if
         end do
         if (active == 0) exit
         do i = 1, active
            node = snapshot(i)
            call node_split(x(lower(node):upper(node)), y(lower(node):upper(node)), g, &
               min_obs, min_rate, binary, direction, side(node), threshold(node), split, rightmean)
            if (split == 0) then
               rightmean = sum(y(lower(node):upper(node))) / real(upper(node) - lower(node) + 1, dp)
               where (state(:nnode) == 0) threshold(:nnode) = rightmean
               state(node) = 2
            else
               split = split + lower(node) - 1
               state(node) = 1
               where (state(:nnode) == 0) threshold(:nnode) = rightmean
               do j = 1, 2
                  k = nnode + j
                  lower(k) = lower(node)
                  upper(k) = upper(node)
                  if (j == 1) upper(k) = split - 1
                  if (j == 2) lower(k) = split
                  side(k) = 2 * j - 3
                  state(k) = 0
                  threshold(k) = rightmean
                  nodecut(k) = x(split)
               end do
               nnode = nnode + 2
            end if
         end do
      end do
      ncut = 0
      do i = 1, nnode
         if (state(i) /= 2 .or. ieee_is_nan(nodecut(i))) cycle
         if (any(cuts(:ncut) == nodecut(i))) cycle
         ncut = ncut + 1
         cuts(ncut) = nodecut(i)
      end do
      call r_order(cuts(:ncut), cutorder)
      cuts(:ncut) = cuts(cutorder)
      k = 1
      do i = 1, size(x)
         j = 1 + count(x(i) >= cuts(:ncut))
         call add_observation(work(j), x(i), y(i))
      end do
      bins = work(:ncut + 1)
   end subroutine tree_bins

   pure subroutine node_split(x, y, g, min_obs, min_rate, binary, direction, side, threshold, split, rightmean)
      real(dp), intent(in) :: x(:) !! Sorted predictor observations in the current tree node.
      real(dp), intent(in) :: y(:) !! Responses belonging to this node.
      integer, intent(in) :: g !! Requested percentile groups for candidate cutpoints.
      integer, intent(in) :: min_obs !! Minimum observations required per child.
      real(dp), intent(in) :: min_rate !! Minimum bad count or continuous target mean per child.
      logical, intent(in) :: binary !! Binary chooses weighted Gini impurity; continuous chooses SSE.
      integer, intent(in) :: direction !! 1 for increasing response, otherwise decreasing.
      integer, intent(in) :: side !! -1 for left child, 1 for right child, 0 for root.
      real(dp), intent(in) :: threshold !! Cross-node monotonic threshold updated breadth-first as in upstream.
      integer, intent(out) :: split !! First right-child observation index, or 0 if no admissible split.
      real(dp), intent(out) :: rightmean !! Right-child response mean for the selected split.
      type(bin_summary), allocatable :: pre(:)
      integer :: n, i, low, high, countleft, countright, sign_dir, sign_thr
      real(dp) :: sleft, sright, mleft, mright, metric, best, total, value
      logical :: ok
      split = 0
      rightmean = 0.0_dp
      call percentile_groups(x, y, g, pre)
      n = size(pre)
      if (n < 4) return
      low = n + 1
      high = 0
      countleft = 0
      sleft = 0.0_dp
      do i = 1, n
         countleft = countleft + pre(i)%count
         sleft = sleft + pre(i)%y_sum
         value = sleft
         if (.not. binary) value = sleft / real(countleft, dp)
         if (countleft >= min_obs .and. value >= min_rate) then
            low = i + 1
            exit
         end if
      end do
      countright = 0
      sright = 0.0_dp
      do i = n, 1, -1
         countright = countright + pre(i)%count
         sright = sright + pre(i)%y_sum
         value = sright
         if (.not. binary) value = sright / real(countright, dp)
         if (countright >= min_obs .and. value >= min_rate) then
            high = i - 2
            exit
         end if
      end do
      if (low > high) return
      total = sum(y)
      sign_dir = -1
      if (direction == 1) sign_dir = 1
      sign_thr = sign_dir
      if (side == -1) sign_thr = -sign_dir
      best = huge(1.0_dp)
      do i = low, high
         countleft = sum(pre(:i - 1)%count)
         countright = size(y) - countleft
         sleft = sum(pre(:i - 1)%y_sum)
         sright = total - sleft
         mleft = sleft / real(countleft, dp)
         mright = sright / real(countright, dp)
         ok = real(sign_dir, dp) * (mright - mleft) > 0.0_dp
         ok = ok .and. real(sign_thr, dp) * (mleft - threshold) > 0.0_dp
         ok = ok .and. real(sign_thr, dp) * (mright - threshold) > 0.0_dp
         if (.not. ok) cycle
         metric = sum((y(:countleft) - mleft) ** 2) + sum((y(countleft + 1:) - mright) ** 2)
         if (metric < best) then
            best = metric
            split = countleft + 1
            rightmean = mright
         end if
      end do
   end subroutine node_split
end module monobin
