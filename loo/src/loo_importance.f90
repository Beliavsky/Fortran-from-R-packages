! SPDX-License-Identifier: GPL-3.0-or-later
module loo_importance
   use loo_types, only: dp, importance_sampling_result, loo_nan
   use matrixstats, only: log_sum_exp
   use r_ordering, only: r_order
   use r_transforms, only: r_expm1
   use posterior_pareto, only: gpdfit, qgeneralized_pareto
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, ieee_positive_inf
   implicit none
   private
   public :: psis, sis, tis, importance_sampling, importance_weights, psislw
   interface psis
      module procedure psis_vector, psis_matrix, psis_array
   end interface psis
   interface sis
      module procedure sis_vector, sis_matrix, sis_array
   end interface sis
   interface tis
      module procedure tis_vector, tis_matrix, tis_array
   end interface tis
contains
   pure function importance_sampling(log_ratios, method, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :) !! Draws by observations; finite values or negative infinity, at least one finite per column.
      character(*), intent(in) :: method !! Importance algorithm: 'psis', 'sis', or 'tis'.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiencies, scalar array or one per observation; NaN becomes one.
      type(importance_sampling_result) :: out
      integer :: n, s, j
      real(dp) :: cutoff, shift, normalizer
      real(dp), allocatable :: centered(:)
      n = size(log_ratios, 2)
      s = size(log_ratios, 1)
      if (n < 1 .or. s < 1) return
      if (method /= 'psis' .and. method /= 'sis' .and. method /= 'tis') return
      if (any(ieee_is_nan(log_ratios))) return
      if (any(log_ratios == ieee_value(0.0_dp, ieee_positive_inf))) return
      do j = 1, n
         if (.not. any(ieee_is_finite(log_ratios(:, j)))) return
      end do
      allocate(out%r_eff(n))
      out%r_eff = 1.0_dp
      if (present(r_eff)) then
         if (size(r_eff) /= 1 .and. size(r_eff) /= n) return
         if (size(r_eff) == 1) then
            out%r_eff = r_eff(1)
         else
            out%r_eff = r_eff
         end if
         where (ieee_is_nan(out%r_eff)) out%r_eff = 1.0_dp
      end if
      if (any(.not. ieee_is_finite(out%r_eff)) .or. any(out%r_eff <= 0.0_dp)) return
      allocate(out%log_weights(s, n), out%normalized_log_weights(s, n), out%norm_const_log(n), out%pareto_k(n), out%n_eff(n), out%tail_len(n))
      out%method = method
      out%log_weights = log_ratios
      out%pareto_k = 0.0_dp
      do j = 1, n
         if (out%r_eff(j) <= 225.0_dp / real(s, dp)) then
            out%tail_len(j) = ceiling(0.2_dp * real(s, dp))
         else
            out%tail_len(j) = ceiling(3.0_dp * sqrt(real(s, dp) / out%r_eff(j)))
         end if
         shift = maxval(log_ratios(:, j))
         centered = log_ratios(:, j) - shift
         if (method == 'psis') then
            block
               real(dp) :: smoothed(s)
               call smooth_column(centered, out%tail_len(j), smoothed, out%pareto_k(j))
               centered = smoothed
            end block
         else if (method == 'tis') then
            cutoff = log_sum_exp(centered) - 0.5_dp * log(real(s, dp))
            centered = min(centered, cutoff)
         end if
         normalizer = log_sum_exp(centered)
         out%log_weights(:, j) = centered + shift
         out%norm_const_log(j) = normalizer + shift
         out%normalized_log_weights(:, j) = centered - normalizer
         out%n_eff(j) = out%r_eff(j) / sum(exp(2.0_dp * out%normalized_log_weights(:, j)))
      end do
      out%status = 0
   end function importance_sampling

   pure subroutine smooth_column(x, tail_length, lw, k)
      real(dp), intent(in) :: x(:) !! Finite or negative-infinite log ratios with at least one finite value.
      integer, intent(in) :: tail_length !! Upper order-statistic tail count, strictly below the sample size.
      real(dp), intent(out) :: lw(size(x)) !! Smoothed, capped log weights in original input order and location.
      real(dp), intent(out) :: k !! Estimated generalized Pareto shape, positive infinity for an unusable tail.
      integer, allocatable :: order(:)
      real(dp), allocatable :: tail(:), excess(:)
      real(dp) :: shift, cutoff, sigma, fitted_k, q
      integer :: n, j, st
      n = size(x)
      shift = maxval(x)
      lw = x - shift
      k = ieee_value(0.0_dp, ieee_positive_inf)
      if (tail_length >= 5 .and. tail_length < n) then
         call r_order(lw, order)
         tail = lw(order(n - tail_length + 1:))
         if (maxval(tail) - minval(tail) >= epsilon(1.0_dp) / 100.0_dp) then
            cutoff = lw(order(n - tail_length))
            ! Stable counterpart of exp(tail)-exp(cutoff), preserving tiny excesses.
            if (ieee_is_finite(cutoff)) then
               excess = exp(tail) * (-r_expm1(cutoff - tail))
            else
               excess = exp(tail)
            end if
            call gpdfit(excess, fitted_k, sigma, sort_x = .false., status = st)
            if (st == 0 .and. ieee_is_finite(fitted_k)) then
               k = fitted_k
               do j = 1, tail_length
                  q = qgeneralized_pareto((real(j, dp) - 0.5_dp) / real(tail_length, dp), &
                     0.0_dp, sigma, k) + exp(cutoff)
                  lw(order(n - tail_length + j)) = log(q)
               end do
            end if
         end if
      end if
      lw = min(lw, 0.0_dp) + shift
   end subroutine smooth_column

   pure function importance_weights(object, logarithm, normalize) result(w)
      type(importance_sampling_result), intent(in) :: object !! Successful importance-sampling result with stored log weights.
      logical, intent(in), optional :: logarithm !! Return log weights when true (default); ordinary weights otherwise.
      logical, intent(in), optional :: normalize !! Normalize each observation to unit total mass; default true.
      real(dp), allocatable :: w(:, :)
      logical :: log_output, normalized
      allocate(w(0, 0))
      if (object%status /= 0) return
      log_output = .true.
      normalized = .true.
      if (present(logarithm)) log_output = logarithm
      if (present(normalize)) normalized = normalize
      w = object%log_weights
      if (normalized) w = object%normalized_log_weights
      if (.not. log_output) w = exp(w)
   end function importance_weights

   pure function psis_vector(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(reshape(log_ratios, [size(log_ratios), 1]), 'psis', r_eff)
   end function psis_vector

   pure function psis_matrix(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(log_ratios, 'psis', r_eff)
   end function psis_matrix

   pure function psis_array(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :, :) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(reshape(log_ratios, [size(log_ratios, 1) * size(log_ratios, 2), size(log_ratios, 3)]), 'psis', r_eff)
   end function psis_array

   pure function sis_vector(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(reshape(log_ratios, [size(log_ratios), 1]), 'sis', r_eff)
   end function sis_vector

   pure function sis_matrix(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(log_ratios, 'sis', r_eff)
   end function sis_matrix

   pure function sis_array(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :, :) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(reshape(log_ratios, [size(log_ratios, 1) * size(log_ratios, 2), size(log_ratios, 3)]), 'sis', r_eff)
   end function sis_array

   pure function tis_vector(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(reshape(log_ratios, [size(log_ratios), 1]), 'tis', r_eff)
   end function tis_vector

   pure function tis_matrix(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(log_ratios, 'tis', r_eff)
   end function tis_matrix

   pure function tis_array(log_ratios, r_eff) result(out)
      real(dp), intent(in) :: log_ratios(:, :, :) !! Log ratios, with draws first; arrays are iterations by chains by observations.
      real(dp), intent(in), optional :: r_eff(:) !! Positive relative efficiency, one value or one per observation; default one.
      type(importance_sampling_result) :: out
      out = importance_sampling(reshape(log_ratios, [size(log_ratios, 1) * size(log_ratios, 2), size(log_ratios, 3)]), 'tis', r_eff)
   end function tis_array

   pure function psislw(log_ratios, wcp, wtrunc) result(out)
      real(dp), intent(in) :: log_ratios(:, :) !! Finite draw-by-observation log ratios for the deprecated PSIS algorithm.
      real(dp), intent(in), optional :: wcp !! Tail proportion in (0,1), default 0.2; uses R type-seven quantile cutoff.
      real(dp), intent(in), optional :: wtrunc !! Nonnegative truncation exponent, default 0.75; zero disables truncation.
      type(importance_sampling_result) :: out
      integer, allocatable :: order(:)
      real(dp), allocatable :: x(:), tail(:), excess(:)
      real(dp) :: prop, truncation, pos, cutoff, k, sigma, q
      integer :: n, j, lo, t, i, st
      prop = 0.2_dp
      truncation = 0.75_dp
      if (present(wcp)) prop = wcp
      if (present(wtrunc)) truncation = wtrunc
      if (.not. ieee_is_finite(prop) .or. .not. ieee_is_finite(truncation)) return
      if (prop <= 0.0_dp .or. prop >= 1.0_dp .or. truncation < 0.0_dp) return
      if (any(.not. ieee_is_finite(log_ratios))) return
      out = sis(log_ratios)
      if (out%status /= 0) return
      n = size(log_ratios, 1)
      out%method = 'psis'
      do j = 1, size(log_ratios, 2)
         x = log_ratios(:, j) - maxval(log_ratios(:, j))
         call r_order(x, order)
         pos = 1.0_dp + real(n - 1, dp) * (1.0_dp - prop)
         lo = floor(pos)
         cutoff = x(order(lo))
         if (lo < n) cutoff = cutoff + (pos - real(lo, dp)) * (x(order(lo + 1)) - cutoff)
         cutoff = max(-700.0_dp, cutoff)
         t = count(x > cutoff)
         out%tail_len(j) = t
         out%pareto_k(j) = ieee_value(0.0_dp, ieee_positive_inf)
         if (t >= 5) then
            tail = x(order(n - t + 1:))
            if (maxval(tail) > minval(tail)) then
               excess = exp(tail) * (-r_expm1(cutoff - tail))
               call gpdfit(excess, k, sigma, wip = .false., min_grid_pts = 80, sort_x = .false., status = st)
               if (st == 0 .and. ieee_is_finite(k)) then
                  out%pareto_k(j) = k
                  do i = 1, t
                     q = qgeneralized_pareto((real(i, dp) - 0.5_dp) / real(t, dp), 0.0_dp, sigma, k)
                     x(order(n - t + i)) = log(q + exp(cutoff))
                  end do
               end if
            end if
         end if
         if (truncation > 0.0_dp) then
            cutoff = (truncation - 1.0_dp) * log(real(n, dp)) + log_sum_exp(x)
            x = min(x, cutoff)
         end if
         out%log_weights(:, j) = x - log_sum_exp(x)
         out%norm_const_log(j) = log_sum_exp(out%log_weights(:, j))
         out%normalized_log_weights(:, j) = out%log_weights(:, j)
         out%n_eff(j) = 1.0_dp / sum(exp(2.0_dp * out%log_weights(:, j)))
      end do
   end function psislw
end module loo_importance
