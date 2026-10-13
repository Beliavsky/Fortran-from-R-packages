! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of PRROC 1.4 by Jan Grau and Jens Keilwagen.
module prroc_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_ordering, only : r_order
   use r_transforms, only : r_log1p, r_log1pexp, r_logistic
   implicit none
   private

   integer, parameter, public :: prroc_ok = 0
   integer, parameter, public :: prroc_invalid_input = 1
   integer, parameter, public :: prroc_degenerate = 2
   integer, parameter, public :: prroc_curve_limit = 3

   type, public :: prroc_summary
      real(dp) :: auc = 0.0_dp
      real(dp) :: auc_davis_goadrich = 0.0_dp
      real(dp), allocatable :: curve(:,:)
   end type prroc_summary

   type, extends(prroc_summary), public :: prroc_result
      type(prroc_summary), allocatable :: maximum
      type(prroc_summary), allocatable :: minimum
      type(prroc_summary), allocatable :: random
   end type prroc_result

   public :: dp, pr_curve, roc_curve

contains

   pure subroutine pr_curve(scores_class0, result, scores_class1, weights_class0, weights_class1, &
         sorted, curve, min_step_size, max_compute, min_compute, rand_compute, dg_compute, status, max_curve_points)
      real(dp), intent(in) :: scores_class0(:) !! Scores of positive observations, or all observations for soft labels.
      type(prroc_result), intent(out) :: result !! Exact PR AUC, optional Davis-Goadrich AUC, curves and references.
      real(dp), intent(in), optional :: scores_class1(:) !! Negative scores; defaults to scores_class0.
      real(dp), intent(in), optional :: weights_class0(:) !! Nonnegative finite positive weights; defaults to ones.
      real(dp), intent(in), optional :: weights_class1(:) !! Negative weights; defaults to 1 - weights_class0, or ones.
      logical, intent(in), optional :: sorted !! If true, score vectors must already be increasing; default false.
      logical, intent(in), optional :: curve !! Return curve rows (recall, precision, threshold); default false.
      real(dp), intent(in), optional :: min_step_size !! Positive interpolation step in original positive-weight units.
      logical, intent(in), optional :: max_compute !! Compute the upstream maximum-reference construction; default false.
      logical, intent(in), optional :: min_compute !! Compute the upstream minimum-reference construction; default false.
      logical, intent(in), optional :: rand_compute !! Return random-reference prevalence and two curve endpoints.
      logical, intent(in), optional :: dg_compute !! Compute Davis-Goadrich AUC when hard-label eligible; default true.
      integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 degenerate classes, 3 curve-size limit.
      integer, intent(in), optional :: max_curve_points !! Maximum rows per curve; default 1000000, at least two.
      call calculate(scores_class0, result, .true., scores_class1, weights_class0, weights_class1, sorted, &
         curve, min_step_size, max_compute, min_compute, rand_compute, dg_compute, status, max_curve_points)
   end subroutine pr_curve

   pure subroutine roc_curve(scores_class0, result, scores_class1, weights_class0, weights_class1, &
         sorted, curve, max_compute, min_compute, rand_compute, status, max_curve_points)
      real(dp), intent(in) :: scores_class0(:) !! Positive scores, or all scores when positive labels are supplied.
      type(prroc_result), intent(out) :: result !! Trapezoidal ROC AUC, optional curve and reference results.
      real(dp), intent(in), optional :: scores_class1(:) !! Negative-class scores; defaults to scores_class0.
      real(dp), intent(in), optional :: weights_class0(:) !! Finite nonnegative positive weights; defaults to ones.
      real(dp), intent(in), optional :: weights_class1(:) !! Negative weights; defaults to 1 - weights_class0, or ones.
      logical, intent(in), optional :: sorted !! Increasing score inputs may skip sorting; checked when true.
      logical, intent(in), optional :: curve !! Return rows (false-positive rate, sensitivity, threshold); default false.
      logical, intent(in), optional :: max_compute !! Include the maximum-reference construction; default false.
      logical, intent(in), optional :: min_compute !! Include the minimum-reference construction; default false.
      logical, intent(in), optional :: rand_compute !! Include random ROC AUC 0.5 and diagonal endpoints.
      integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 degenerate classes, 3 curve-size limit.
      integer, intent(in), optional :: max_curve_points !! Maximum rows per curve; default 1000000, at least two.
      call calculate(scores_class0, result, .false., scores_class1, weights_class0, weights_class1, sorted, &
         curve, max_compute = max_compute, min_compute = min_compute, rand_compute = rand_compute, &
         status = status, max_curve_points = max_curve_points)
   end subroutine roc_curve

   pure subroutine calculate(scores0, result, is_pr, scores1, weights0, weights1, sorted, curve, &
         min_step_size, max_compute, min_compute, rand_compute, dg_compute, status, max_curve_points)
      real(dp), intent(in) :: scores0(:) !! Supplied first score vector, interpreted as the positive class.
      type(prroc_result), intent(out) :: result !! Main numerical result and requested reference summaries.
      logical, intent(in) :: is_pr !! True selects precision-recall; false selects ROC.
      real(dp), intent(in), optional :: scores1(:) !! Second score vector, defaulting to scores0.
      real(dp), intent(in), optional :: weights0(:) !! Optional first vector of finite nonnegative weights.
      real(dp), intent(in), optional :: weights1(:) !! Optional second vector of finite nonnegative weights.
      logical, intent(in), optional :: sorted !! Whether supplied scores are in increasing order.
      logical, intent(in), optional :: curve !! Whether observed and nonrandom reference points are requested.
      real(dp), intent(in), optional :: min_step_size !! Positive weight increment for PR interpolation.
      logical, intent(in), optional :: max_compute !! Whether to calculate the best-reference ordering.
      logical, intent(in), optional :: min_compute !! Whether to calculate the worst-reference ordering.
      logical, intent(in), optional :: rand_compute !! Whether to include a random-classifier reference.
      logical, intent(in), optional :: dg_compute !! Whether eligible unweighted Davis-Goadrich AUC is requested.
      integer, intent(out), optional :: status !! Summary of input, class-degeneracy or curve-limit errors.
      integer, intent(in), optional :: max_curve_points !! Per-curve row allocation limit, default one million.
      real(dp), allocatable :: scores(:), pos(:), neg(:), best(:), worst(:)
      real(dp) :: scale, step, total_pos, total_neg, log_weight_ratio
      logical :: want_curve, want_dg, eligible
      integer :: code, reference_code, limit

      call clear_summary(result%prroc_summary)
      if (present(status)) status = prroc_invalid_input
      want_curve = .false.
      if (present(curve)) want_curve = curve
      want_dg = .true.
      if (present(dg_compute)) want_dg = dg_compute
      limit = 1000000
      if (present(max_curve_points)) limit = max_curve_points
      if (limit < 2) return
      if (present(min_step_size)) then
         if (.not. ieee_is_finite(min_step_size)) return
         if (min_step_size <= 0.0_dp) return
      end if
      call prepare(scores0, scores1, weights0, weights1, sorted, scores, pos, neg, best, worst, &
         scale, log_weight_ratio, eligible, code)
      if (code /= prroc_ok) then
         if (present(status)) status = code
         return
      end if
      total_pos = sum(pos)
      total_neg = sum(neg)
      if (total_pos <= 0.0_dp .or. (.not. is_pr .and. total_neg <= 0.0_dp)) then
         if (present(status)) status = prroc_degenerate
         return
      end if
      if (scale < 1.0_dp) then
         step = huge(1.0_dp)
         if (scale >= 1.0_dp / huge(1.0_dp)) step = 1.0_dp / scale
      else
         step = 1.0_dp / scale
      end if
      if (present(weights0)) then
         step = total_pos / 100.0_dp
         if (scale >= 1.0_dp) then
            step = min(step, 1.0_dp / scale)
         else if (total_pos * scale > 100.0_dp) then
            step = 1.0_dp / scale
         end if
      end if
      if (present(min_step_size)) then
         if (scale < 1.0_dp) then
            if (min_step_size > huge(step) * scale) then
               step = huge(step)
            else
               step = min_step_size / scale
            end if
         else
            step = min_step_size / scale
         end if
      end if
      call evaluate(scores, pos, neg, is_pr, eligible .and. want_dg, want_curve, step, scale, log_weight_ratio, limit, &
         result%prroc_summary, code)
      if (present(max_compute)) then
         if (max_compute) then
            allocate(result%maximum)
            call evaluate(best, pos, neg, is_pr, eligible .and. want_dg, want_curve, step, scale, log_weight_ratio, limit, &
               result%maximum, reference_code)
            code = max(code, reference_code)
         end if
      end if
      if (present(min_compute)) then
         if (min_compute) then
            allocate(result%minimum)
            call evaluate(worst, pos, neg, is_pr, eligible .and. want_dg, want_curve, step, scale, log_weight_ratio, limit, &
               result%minimum, reference_code)
            code = max(code, reference_code)
         end if
      end if
      if (present(rand_compute)) then
         if (rand_compute) then
            allocate(result%random)
            call clear_summary(result%random)
            deallocate(result%random%curve)
            allocate(result%random%curve(2, 3))
            result%random%curve = 0.0_dp
            result%random%curve(2, 1) = 1.0_dp
            if (is_pr) then
               result%random%auc = weighted_precision(total_pos, total_neg, log_weight_ratio)
               result%random%auc_davis_goadrich = result%random%auc
               result%random%curve(:, 2) = result%random%auc
            else
               result%random%auc = 0.5_dp
               result%random%curve(2, 2) = 1.0_dp
            end if
         end if
      end if
      if (present(status)) status = code
   end subroutine calculate

   pure subroutine clear_summary(result)
      type(prroc_summary), intent(out) :: result !! Reset AUCs to quiet NaN and create a zero-row curve.
      result%auc = ieee_value(0.0_dp, ieee_quiet_nan)
      result%auc_davis_goadrich = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(result%curve(0, 3))
   end subroutine clear_summary

   pure subroutine prepare(scores0, scores1, weights0, weights1, sorted, scores, pos, neg, best, worst, &
         scale, log_weight_ratio, eligible, status)
      real(dp), intent(in) :: scores0(:) !! First score vector; NaNs are invalid, infinities are permitted.
      real(dp), intent(in), optional :: scores1(:) !! Second scores, otherwise the same observations as scores0.
      real(dp), intent(in), optional :: weights0(:) !! First weights, otherwise all ones.
      real(dp), intent(in), optional :: weights1(:) !! Second weights, otherwise complementary or all ones.
      logical, intent(in), optional :: sorted !! True requests an ascending-order validation instead of sorting.
      real(dp), allocatable, intent(out) :: scores(:) !! Combined observed scores with zero-weight rows removed.
      real(dp), allocatable, intent(out) :: pos(:) !! Positive masses divided by their class maximum.
      real(dp), allocatable, intent(out) :: neg(:) !! Negative masses divided by their class maximum.
      real(dp), allocatable, intent(out) :: best(:) !! Upstream maximum-reference scores for each combined row.
      real(dp), allocatable, intent(out) :: worst(:) !! Upstream minimum-reference scores for each combined row.
      real(dp), intent(out) :: scale !! Maximum original positive weight, used for interpolation steps.
      real(dp), intent(out) :: log_weight_ratio !! Logarithm of negative-to-positive weight scale ratio.
      logical, intent(out) :: eligible !! Whether R's complementary binary-label Davis-Goadrich test passes.
      integer, intent(out) :: status !! Input validation or all-zero-weight status.
      real(dp), allocatable :: s0(:), s1(:), w0(:), w1(:), all_pos(:), all_neg(:), all_best(:), all_worst(:)
      integer, allocatable :: order(:)
      logical, allocatable :: keep(:)
      logical :: same, distinct_unweighted, already_sorted
      integer :: n0, n1
      real(dp) :: negative_scale

      status = prroc_invalid_input
      eligible = .false.
      scale = 1.0_dp
      log_weight_ratio = 0.0_dp
      s0 = scores0
      n0 = size(s0)
      n1 = n0
      if (present(scores1)) n1 = size(scores1)
      if (n0 + n1 == 0) return
      if (any(ieee_is_nan(s0))) return
      allocate(w0(n0), w1(n1))
      w0 = 1.0_dp
      if (present(weights0)) then
         if (size(weights0) /= n0) return
         w0 = weights0
      end if
      if (any(.not. ieee_is_finite(w0))) return
      if (any(w0 < 0.0_dp)) return
      already_sorted = .false.
      if (present(sorted)) already_sorted = sorted
      if (already_sorted) then
         if (any(s0(2:) < s0(:n0 - 1))) return
      else
         call r_order(s0, order)
         s0 = s0(order)
         w0 = w0(order)
      end if
      ! Match R's lazy defaults: absent second scores/weights resolve after the first sort.
      s1 = s0
      if (present(scores1)) s1 = scores1
      if (any(ieee_is_nan(s1))) return
      w1 = 1.0_dp
      if (present(weights1)) then
         if (size(weights1) /= n1) return
         w1 = weights1
      else if (present(weights0)) then
         if (n0 /= n1) return
         w1 = 1.0_dp - w0
      end if
      if (any(.not. ieee_is_finite(w1))) return
      if (any(w1 < 0.0_dp)) return
      if (already_sorted) then
         if (any(s1(2:) < s1(:n1 - 1))) return
      else
         call r_order(s1, order)
         s1 = s1(order)
         w1 = w1(order)
      end if
      same = .false.
      if (n0 == n1) same = all(s0 == s1)
      distinct_unweighted = .not. same .and. .not. present(weights0) .and. .not. present(weights1)
      eligible = distinct_unweighted
      if (same) eligible = all(w0 == 1.0_dp - w1) .and. all(w0 == 0.0_dp .or. w0 == 1.0_dp)
      scale = maxval(w0)
      negative_scale = maxval(w1)
      if (scale <= 0.0_dp) then
         status = prroc_degenerate
         return
      end if
      if (negative_scale <= 0.0_dp) negative_scale = scale
      log_weight_ratio = log(negative_scale) - log(scale)
      allocate(all_pos(n0 + n1), all_neg(n0 + n1), all_best(n0 + n1), all_worst(n0 + n1))
      all_pos = 0.0_dp
      all_neg = 0.0_dp
      all_pos(:n0) = w0 / scale
      all_neg(n0 + 1:) = w1 / negative_scale
      if (same) then
         all_best = [w0, w0]
         all_worst = -all_best
      else
         all_best(:n0) = 1.0_dp
         all_best(n0 + 1:) = 0.0_dp
         if (distinct_unweighted) then
            all_worst = -all_best
         else
            all_worst(:n0) = 0.0_dp
            all_worst(n0 + 1:) = 1.0_dp
         end if
      end if
      keep = all_pos > 0.0_dp .or. all_neg > 0.0_dp
      scores = pack([s0, s1], keep)
      pos = pack(all_pos, keep)
      neg = pack(all_neg, keep)
      best = pack(all_best, keep)
      worst = pack(all_worst, keep)
      status = prroc_ok
   end subroutine prepare

   pure subroutine evaluate(scores, pos, neg, is_pr, do_dg, want_curve, step, scale, log_weight_ratio, limit, result, status)
      real(dp), intent(in) :: scores(:) !! Combined score vector whose equal thresholds are grouped exactly.
      real(dp), intent(in) :: pos(:) !! Positive masses normalized within their class, conformable with scores.
      real(dp), intent(in) :: neg(:) !! Negative masses normalized within their class, conformable with scores.
      logical, intent(in) :: is_pr !! Select nonlinear PR integration instead of ROC trapezoids.
      logical, intent(in) :: do_dg !! Binary-label data permit Davis-Goadrich unit-positive interpolation.
      logical, intent(in) :: want_curve !! Whether to allocate the optional curve matrix.
      real(dp), intent(in) :: step !! PR interpolation step in normalized positive-mass units.
      real(dp), intent(in) :: scale !! Original positive-class weight scale; one for eligible binary labels.
      real(dp), intent(in) :: log_weight_ratio !! Log ratio between negative and positive original weight scales.
      integer, intent(in) :: limit !! Maximum permitted number of returned curve rows.
      type(prroc_summary), intent(out) :: result !! Calculated area, eligible secondary area, and curve rows.
      integer, intent(out) :: status !! Success or a curve-size limit; AUCs survive a curve-size limit.
      integer, allocatable :: order(:), steps(:)
      real(dp), allocatable :: threshold(:), p(:), n(:), dp_group(:), dn_group(:), points(:,:)
      integer :: i, j, g, groups, rows, k, count_dg
      real(dp) :: total_pos, total_neg, pp, nn, increment, fraction, previous_precision, precision, ratio
      logical :: interpolate

      call clear_summary(result)
      status = prroc_ok
      call r_order(scores, order)
      allocate(threshold(size(scores)), p(size(scores)), n(size(scores)), &
         dp_group(size(scores)), dn_group(size(scores)))
      groups = 0
      do i = size(scores), 1, -1
         j = order(i)
         if (groups == 0) then
            groups = 1
            threshold(groups) = scores(j)
            dp_group(groups) = 0.0_dp
            dn_group(groups) = 0.0_dp
         else if (scores(j) /= threshold(groups)) then
            groups = groups + 1
            threshold(groups) = scores(j)
            dp_group(groups) = 0.0_dp
            dn_group(groups) = 0.0_dp
         end if
         dp_group(groups) = dp_group(groups) + pos(j)
         dn_group(groups) = dn_group(groups) + neg(j)
      end do
      pp = 0.0_dp
      nn = 0.0_dp
      do g = 1, groups
         pp = pp + dp_group(g)
         nn = nn + dn_group(g)
         p(g) = pp
         n(g) = nn
      end do
      total_pos = pp
      total_neg = nn
      result%auc = 0.0_dp
      if (do_dg .and. is_pr) result%auc_davis_goadrich = 0.0_dp
      pp = 0.0_dp
      nn = 0.0_dp
      do g = 1, groups
         if (is_pr) then
            result%auc = result%auc + dp_group(g) / total_pos * &
               segment_precision(pp, nn, dp_group(g), dn_group(g), log_weight_ratio)
            if (do_dg .and. dp_group(g) > 0.0_dp) then
               count_dg = max(1, nint(dp_group(g) * scale))
               previous_precision = weighted_precision(p(g), n(g), log_weight_ratio)
               if (pp + nn > 0.0_dp) previous_precision = weighted_precision(pp, nn, log_weight_ratio)
               do k = 1, count_dg
                  fraction = real(k, dp) / real(count_dg, dp)
                  precision = weighted_precision(pp + fraction * dp_group(g), nn + fraction * dn_group(g), &
                     log_weight_ratio)
                  result%auc_davis_goadrich = result%auc_davis_goadrich + &
                     (previous_precision + precision) * 0.5_dp * dp_group(g) / total_pos / real(count_dg, dp)
                  previous_precision = precision
               end do
            end if
         else
            result%auc = result%auc + dn_group(g) / total_neg * (pp / total_pos + p(g) / total_pos) * 0.5_dp
         end if
         pp = p(g)
         nn = n(g)
      end do
      if (.not. want_curve) return
      if (.not. is_pr) then
         if (groups >= limit) then
            status = prroc_curve_limit
            return
         end if
         deallocate(result%curve)
         allocate(result%curve(groups + 1, 3))
         result%curve(1, :) = [0.0_dp, 0.0_dp, threshold(1)]
         do g = 1, groups
            result%curve(g + 1, :) = [n(g) / total_neg, p(g) / total_pos, threshold(min(g + 1, groups))]
         end do
         return
      end if
      allocate(steps(groups))
      steps = 0
      rows = groups + 2
      if (rows > limit) then
         status = prroc_curve_limit
         return
      end if
      do g = 2, groups
         interpolate = dp_group(g) > step .and. weighted_precision(p(g), n(g), log_weight_ratio) /= &
            weighted_precision(p(g - 1), n(g - 1), log_weight_ratio)
         if (.not. interpolate) cycle
         if (step <= 0.0_dp) then
            status = prroc_curve_limit
            return
         end if
         if (dp_group(g) / real(limit, dp) > step) then
            status = prroc_curve_limit
            return
         end if
         ratio = dp_group(g) / step
         if (ratio + 1.0e-12_dp >= real(limit - rows, dp)) then
            status = prroc_curve_limit
            return
         end if
         steps(g) = floor(ratio + 1.0e-12_dp) + 1
         if (steps(g) > limit - rows) then
            status = prroc_curve_limit
            return
         end if
         rows = rows + steps(g)
      end do
      allocate(points(rows - 2, 3))
      j = 0
      do g = 1, groups
         do k = 0, steps(g) - 1
            increment = min(real(k, dp) * step, dp_group(g))
            fraction = increment / dp_group(g)
            j = j + 1
            points(j, :) = [(p(g - 1) + increment) / total_pos, &
               weighted_precision(p(g - 1) + increment, n(g - 1) + fraction * dn_group(g), &
               log_weight_ratio), threshold(g)]
         end do
         j = j + 1
         points(j, :) = [p(g) / total_pos, weighted_precision(p(g), n(g), log_weight_ratio), threshold(g)]
      end do
      deallocate(result%curve)
      allocate(result%curve(rows, 3))
      result%curve(2:rows - 1, :) = points(rows - 2:1:-1, :)
      result%curve(1, :) = [1.0_dp, points(rows - 2, 2), points(rows - 2, 3)]
      result%curve(rows, :) = [0.0_dp, points(1, 2), points(1, 3)]
   end subroutine evaluate

   pure elemental real(dp) function weighted_precision(p, n, log_weight_ratio) result(value)
      real(dp), intent(in) :: p !! Nonnegative positive mass in positive-class scale units.
      real(dp), intent(in) :: n !! Nonnegative negative mass in negative-class scale units.
      real(dp), intent(in) :: log_weight_ratio !! Logarithm of negative-to-positive class scale ratio.
      if (p <= 0.0_dp) then
         value = 0.0_dp
      else if (n <= 0.0_dp) then
         value = 1.0_dp
      else
         value = r_logistic(log(p) - log(n) - log_weight_ratio)
      end if
   end function weighted_precision

   pure elemental real(dp) function log_mass(p, n, log_weight_ratio) result(value)
      real(dp), intent(in) :: p !! Nonnegative positive mass; at least one mass must be positive.
      real(dp), intent(in) :: n !! Nonnegative negative mass in its separately normalized units.
      real(dp), intent(in) :: log_weight_ratio !! Logarithm of negative-to-positive mass unit ratio.
      real(dp) :: lp, ln
      if (p <= 0.0_dp) then
         value = log(n) + log_weight_ratio
      else if (n <= 0.0_dp) then
         value = log(p)
      else
         lp = log(p)
         ln = log(n) + log_weight_ratio
         value = max(lp, ln) + r_log1p(exp(-abs(lp - ln)))
      end if
   end function log_mass

   pure elemental real(dp) function segment_precision(p0, n0, dp_segment, dn_segment, log_weight_ratio) result(value)
      real(dp), intent(in) :: p0 !! Cumulative positive mass before this threshold segment.
      real(dp), intent(in) :: n0 !! Cumulative negative mass before this threshold segment.
      real(dp), intent(in) :: dp_segment !! Nonnegative positive mass added at this threshold.
      real(dp), intent(in) :: dn_segment !! Nonnegative negative mass added at this threshold.
      real(dp), intent(in) :: log_weight_ratio !! Logarithm of negative-to-positive normalized mass units.
      real(dp) :: slope_precision, log_ratio, complement, fraction, log_fraction
      if (dp_segment + dn_segment <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      slope_precision = weighted_precision(dp_segment, dn_segment, log_weight_ratio)
      if (p0 + n0 <= 0.0_dp) then
         value = slope_precision
      else
         log_fraction = log_mass(dp_segment, dn_segment, log_weight_ratio) - log_mass(p0, n0, log_weight_ratio)
         if (log_fraction < 0.0_dp) then
            fraction = exp(log_fraction)
            if (fraction < 1.0e-8_dp) then
               complement = fraction * (0.5_dp - fraction / 3.0_dp)
               log_ratio = 1.0_dp - complement
            else
               log_ratio = r_log1p(fraction) / fraction
               complement = 1.0_dp - log_ratio
            end if
         else
            log_ratio = exp(-log_fraction) * r_log1pexp(log_fraction)
            complement = 1.0_dp - log_ratio
         end if
         value = slope_precision * complement + weighted_precision(p0, n0, log_weight_ratio) * log_ratio
      end if
   end function segment_precision

end module prroc_api
