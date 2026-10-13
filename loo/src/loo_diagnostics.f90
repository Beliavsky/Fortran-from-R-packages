! SPDX-License-Identifier: GPL-3.0-or-later
module loo_diagnostics
   use loo_types
   use loo_importance, only: importance_weights
   use r_ordering, only: r_order
   use r_quantiles, only: r_quantile_type7
   use posterior_pareto, only: pareto_khat
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_negative_inf
   implicit none
   private
   public :: E_loo, expectation_result, pareto_k_ids, pareto_k_table, mcse_loo
   public :: loo_compare, comparison_result, compare, pareto_threshold
   type :: expectation_result
      real(dp), allocatable :: value(:, :)
      real(dp), allocatable :: pareto_k(:)
      integer :: status = 1
   end type expectation_result
   type :: comparison_result
      integer, allocatable :: order(:)
      real(dp), allocatable :: elpd_diff(:), se_diff(:), p_worse(:)
      ! Sorted model totals and their standard errors.
      real(dp), allocatable :: elpd(:), se_elpd(:)
      integer :: status = 1
   end type comparison_result
   interface E_loo
      module procedure E_loo_matrix, E_loo_vector
   end interface E_loo
contains
   pure elemental function pareto_threshold(n_draws) result(threshold)
      integer, intent(in) :: n_draws !! Positive number of posterior draws; one yields negative infinity.
      real(dp) :: threshold
      threshold = loo_nan()
      if (n_draws == 1) threshold = ieee_value(0.0_dp, ieee_negative_inf)
      if (n_draws > 1) threshold = min(1.0_dp - 1.0_dp / log10(real(n_draws, dp)), 0.7_dp)
   end function pareto_threshold

   pure function E_loo_matrix(x, object, statistic, probs, log_ratios) result(out)
      real(dp), intent(in) :: x(:, :) !! Predictive quantity, posterior draws by observations, matching the PSIS object.
      type(importance_sampling_result), intent(in) :: object !! Successful PSIS object with normalized weights and Pareto tail lengths.
      character(*), intent(in), optional :: statistic !! 'mean' (default), 'variance', 'sd', or 'quantile'.
      real(dp), intent(in), optional :: probs(:) !! Probabilities strictly between zero and one, required for quantile.
      real(dp), intent(in), optional :: log_ratios(:, :) !! Raw log ratios for improved function-specific k diagnostics; default smoothed ratios.
      type(expectation_result) :: out
      real(dp), allocatable :: w(:, :), lr(:, :), h(:), ratio(:), cumulative(:)
      integer, allocatable :: order(:)
      character(len = 8) :: mode
      real(dp) :: mean_x, value, den, k_r, k_hr, p
      integer :: j, q, i, nq, unique_count
      if (object%status /= 0) return
      if (object%method /= 'psis') return
      if (any(shape(x) /= shape(object%log_weights))) return
      mode = 'mean'
      if (present(statistic)) mode = statistic
      if (mode /= 'mean' .and. mode /= 'variance' .and. mode /= 'sd' .and. mode /= 'quantile') return
      nq = 1
      if (mode == 'quantile') then
         if (.not. present(probs)) return
         if (size(probs) < 1) return
         if (any(.not. ieee_is_finite(probs))) return
         if (any(probs <= 0.0_dp) .or. any(probs >= 1.0_dp)) return
         nq = size(probs)
      end if
      lr = importance_weights(object)
      if (present(log_ratios)) then
         if (any(shape(log_ratios) /= shape(x))) return
         if (any(.not. ieee_is_finite(log_ratios))) return
         lr = log_ratios
      end if
      w = importance_weights(object, logarithm = .false.)
      allocate(out%value(nq, size(x, 2)), out%pareto_k(size(x, 2)))
      out%value = loo_nan()
      do j = 1, size(x, 2)
         if (all(ieee_is_finite(x(:, j)))) then
            ! Centering avoids the cancellation of E[x**2]-E[x]**2.
            mean_x = x(1, j) + sum(w(:, j) * (x(:, j) - x(1, j))) / sum(w(:, j))
            select case (mode)
            case ('mean')
               out%value(1, j) = mean_x
            case ('variance', 'sd')
               den = 1.0_dp - sum(w(:, j) ** 2)
               if (den > 0.0_dp) then
                  value = sum(w(:, j) * (x(:, j) - mean_x) ** 2) / den
                  if (mode == 'sd') value = sqrt(value)
                  out%value(1, j) = value
               end if
            case ('quantile')
               call r_order(x(:, j), order)
               cumulative = w(order, j)
               do i = 2, size(cumulative)
                  cumulative(i) = cumulative(i) + cumulative(i - 1)
               end do
               cumulative = cumulative / cumulative(size(cumulative))
               do q = 1, nq
                  p = probs(q)
                  if (all(w(:, j) == w(1, j))) then
                     out%value(q, j) = r_quantile_type7(x(:, j), p)
                  else
                     i = 1
                     do while (i < size(cumulative))
                        if (cumulative(i) >= p) exit
                        i = i + 1
                     end do
                     value = x(order(i), j)
                     if (i > 1) value = x(order(i - 1), j) + (value - x(order(i - 1), j)) * &
                        (p - cumulative(i - 1)) / (cumulative(i) - cumulative(i - 1))
                     out%value(q, j) = value
                  end if
               end do
            end select
         end if
         ratio = exp(lr(:, j) - maxval(lr(:, j)))
         k_r = pareto_khat(ratio, tail = 'right', ndraws_tail = object%tail_len(j))
         out%pareto_k(j) = k_r
         if (mode == 'quantile' .or. any(.not. ieee_is_finite(x(:, j)))) cycle
         ! The release skips function-specific k for constant or two-valued quantities.
         call r_order(x(:, j), order)
         unique_count = 1
         do i = 2, size(order)
            if (x(order(i), j) /= x(order(i - 1), j)) unique_count = unique_count + 1
         end do
         if (unique_count <= 2) cycle
         h = x(:, j)
         if (mode == 'variance' .or. mode == 'sd') h = h ** 2
         k_hr = pareto_khat(h * ratio, tail = 'both', ndraws_tail = object%tail_len(j))
         if (ieee_is_nan(k_r)) then
            out%pareto_k(j) = k_hr
         else if (.not. ieee_is_nan(k_hr)) then
            out%pareto_k(j) = max(k_r, k_hr)
         end if
      end do
      out%status = 0
   end function E_loo_matrix

   pure function E_loo_vector(x, object, statistic, probs, log_ratios) result(out)
      real(dp), intent(in) :: x(:) !! Predictive quantity for one observation, one value per PSIS draw.
      type(importance_sampling_result), intent(in) :: object !! PSIS result containing exactly one observation.
      character(*), intent(in), optional :: statistic !! 'mean', 'variance', 'sd', or 'quantile'; default mean.
      real(dp), intent(in), optional :: probs(:) !! Strictly interior quantile probabilities when requested.
      real(dp), intent(in), optional :: log_ratios(:) !! Raw log importance ratios with the same length as x.
      type(expectation_result) :: out
      if (present(log_ratios)) then
         out = E_loo_matrix(reshape(x, [size(x), 1]), object, statistic, probs, &
            reshape(log_ratios, [size(log_ratios), 1]))
      else
         out = E_loo_matrix(reshape(x, [size(x), 1]), object, statistic, probs)
      end if
   end function E_loo_vector

   pure function pareto_k_ids(object, threshold) result(ids)
      type(importance_sampling_result), intent(in) :: object !! Successful importance result whose Pareto estimates are to be checked.
      real(dp), intent(in), optional :: threshold !! Maximum reliable k; default depends on the number of draws, capped at 0.7.
      integer, allocatable :: ids(:)
      real(dp) :: limit
      integer :: i
      allocate(ids(0))
      if (object%status /= 0) return
      limit = pareto_threshold(size(object%log_weights, 1))
      if (present(threshold)) limit = threshold
      ids = pack([(i, i = 1,size(object%pareto_k))], object%pareto_k > limit)
   end function pareto_k_ids

   pure function pareto_k_table(object) result(table)
      type(importance_sampling_result), intent(in) :: object !! Importance result with per-observation k and effective sample size estimates.
      real(dp) :: table(3, 3), limit
      integer :: i, group, total
      table = loo_nan()
      if (object%status /= 0) return
      table(:, 1:2) = 0.0_dp
      limit = pareto_threshold(size(object%log_weights, 1))
      do i = 1, size(object%pareto_k)
         if (ieee_is_nan(object%pareto_k(i))) cycle
         if (object%pareto_k(i) <= limit) then
            group = 1
         else if (object%pareto_k(i) <= 1.0_dp) then
            group = 2
         else
            group = 3
         end if
         table(group, 1) = table(group, 1) + 1.0_dp
         if (group == 1) then
            if (ieee_is_nan(table(1, 3))) then
               table(1, 3) = object%n_eff(i)
            else
               table(1, 3) = min(table(1, 3), object%n_eff(i))
            end if
         end if
      end do
      total = nint(sum(table(:, 1)))
      if (total > 0) table(:, 2) = table(:, 1) / real(total, dp)
   end function pareto_k_table

   pure function mcse_loo(object, threshold) result(value)
      type(loo_result), intent(in) :: object !! Successful LOO result including pointwise Monte Carlo errors and k diagnostics.
      real(dp), intent(in), optional :: threshold !! Reliability cutoff; default sample-size-dependent threshold capped at 0.7.
      real(dp) :: value, limit
      value = loo_nan()
      if (object%status /= 0 .or. object%method /= 'loo') return
      limit = pareto_threshold(object%n_draws)
      if (present(threshold)) limit = threshold
      if (any(object%importance%pareto_k > limit)) return
      value = sqrt(sum(object%pointwise(:, 2) ** 2))
   end function mcse_loo

   pure function loo_compare(pointwise) result(out)
      real(dp), intent(in) :: pointwise(:, :) !! Pointwise predictive log scores: observations by at least two models on the same data.
      type(comparison_result) :: out
      real(dp) :: totals(size(pointwise, 2)), delta(size(pointwise, 1)), mean_delta, se
      integer :: m, n, j, model, best
      m = size(pointwise, 2)
      n = size(pointwise, 1)
      if (m < 2 .or. n < 2) return
      if (any(.not. ieee_is_finite(pointwise))) return
      totals = sum(pointwise, dim = 1)
      call r_order(-totals, out%order)
      allocate(out%elpd_diff(m), out%se_diff(m), out%p_worse(m), out%elpd(m), out%se_elpd(m))
      best = out%order(1)
      do j = 1, m
         model = out%order(j)
         delta = pointwise(:, model) - pointwise(:, best)
         mean_delta = sum(delta) / real(n, dp)
         out%elpd_diff(j) = sum(delta)
         se = sqrt(real(n, dp) * sum((delta - mean_delta) ** 2) / real(n - 1, dp))
         out%se_diff(j) = se
         out%p_worse(j) = loo_nan()
         if (out%elpd_diff(j) /= 0.0_dp) then
            if (se > 0.0_dp) then
               out%p_worse(j) = 0.5_dp * erfc(out%elpd_diff(j) / (sqrt(2.0_dp) * se))
            else
               out%p_worse(j) = 1.0_dp
            end if
         end if
         out%elpd(j) = totals(model)
         mean_delta = totals(model) / real(n, dp)
         out%se_elpd(j) = sqrt(real(n, dp) * sum((pointwise(:, model) - mean_delta) ** 2) / real(n - 1, dp))
      end do
      out%status = 0
   end function loo_compare

   pure function compare(first, second) result(difference)
      real(dp), intent(in) :: first(:) !! Pointwise log predictive scores from the first model.
      real(dp), intent(in) :: second(:) !! Scores from the second model on the identical observations.
      real(dp) :: difference(2), delta(size(first)), mean_delta
      integer :: n
      difference = loo_nan()
      n = size(first)
      if (size(second) /= n .or. n < 2) return
      if (any(.not. ieee_is_finite(first)) .or. any(.not. ieee_is_finite(second))) return
      delta = second - first
      mean_delta = sum(delta) / real(n, dp)
      difference(1) = sum(delta)
      difference(2) = sqrt(real(n, dp) * sum((delta - mean_delta) ** 2) / real(n - 1, dp))
   end function compare
end module loo_diagnostics
