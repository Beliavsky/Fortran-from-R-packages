! SPDX-License-Identifier: GPL-3.0-or-later
! Numerical translation of loo 2.10.1, R/crps.R and R/loo_predictive_metric.R.
! See the package NOTICE and upstream license for attribution.
module loo_metrics
   use r_kinds, only: dp
   use r_descriptive, only: r_mean, r_sd
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   implicit none
   private
   public :: predictive_metric_result, crps, scrps, loo_crps, loo_scrps, loo_predictive_metric

   type :: predictive_metric_result
      real(dp) :: estimate = 0.0_dp
      real(dp) :: se = 0.0_dp
      real(dp), allocatable :: pointwise(:)
      integer :: status = 0 ! 0 success, 1 invalid argument, 2 undefined scaled score or missing class.
   end type predictive_metric_result

   interface crps
      module procedure crps_matrix, crps_vector
   end interface
   interface scrps
      module procedure scrps_matrix, scrps_vector
   end interface
contains
   ! Scores have the upstream sign convention: larger is better. Unlike R's
   ! randomized default, omitted permutations preserve the supplied pairing.
   pure function crps_matrix(x, x2, y, permutations) result(out)
      real(dp), intent(in) :: x(:, :) !! Predictive draws, shaped draws by observations, all finite.
      real(dp), intent(in) :: x2(:, :) !! Independent predictive draws of the same shape as x.
      real(dp), intent(in) :: y(:) !! Observed responses, one per column of x.
      integer, intent(in), optional :: permutations(:, :) !! Each column is a permutation of 1:S; absent means identity.
      type(predictive_metric_result) :: out
      out = score_core(x, x2, y, .false., permutations)
   end function crps_matrix

   pure function scrps_matrix(x, x2, y, permutations) result(out)
      real(dp), intent(in) :: x(:, :) !! Predictive draws, shaped draws by observations, all finite.
      real(dp), intent(in) :: x2(:, :) !! Independent predictive draws of the same shape as x.
      real(dp), intent(in) :: y(:) !! Observed responses, one per column of x.
      integer, intent(in), optional :: permutations(:, :) !! Each column is a permutation of 1:S; absent means identity.
      type(predictive_metric_result) :: out
      out = score_core(x, x2, y, .true., permutations)
   end function scrps_matrix

   pure function crps_vector(x, x2, y, permutations) result(out)
      real(dp), intent(in) :: x(:) !! Finite predictive draws for a single response.
      real(dp), intent(in) :: x2(:) !! Independent predictive draws, the same length as x.
      real(dp), intent(in) :: y !! Single observed response.
      integer, intent(in), optional :: permutations(:, :) !! Permutation columns for x2; absent means identity.
      type(predictive_metric_result) :: out
      out = score_core(reshape(x, [size(x), 1]), reshape(x2, [size(x2), 1]), [y], .false., permutations)
   end function crps_vector

   pure function scrps_vector(x, x2, y, permutations) result(out)
      real(dp), intent(in) :: x(:) !! Finite predictive draws for a single response.
      real(dp), intent(in) :: x2(:) !! Independent predictive draws, the same length as x.
      real(dp), intent(in) :: y !! Single observed response.
      integer, intent(in), optional :: permutations(:, :) !! Permutation columns for x2; absent means identity.
      type(predictive_metric_result) :: out
      out = score_core(reshape(x, [size(x), 1]), reshape(x2, [size(x2), 1]), [y], .true., permutations)
   end function scrps_vector

   ! Supply PSIS-normalized marginal weights for -log_lik, and separately
   ! normalized joint weights for -log_lik-log_lik(permutation,:). The joint
   ! weights are NOT the product of the marginal PSIS-smoothed weights.
   pure function loo_crps(x, x2, y, weights, joint_weights, permutations) result(out)
      real(dp), intent(in) :: x(:, :) !! Predictive draws, shaped S draws by N observations.
      real(dp), intent(in) :: x2(:, :) !! Independent predictive draws of the same shape as x.
      real(dp), intent(in) :: y(:) !! Observed responses, length N.
      real(dp), intent(in) :: weights(:, :) !! Nonnegative normalized ordinary marginal LOO weights, shape S by N.
      real(dp), intent(in) :: joint_weights(:, :, :) !! Normalized paired LOO weights, shape S by N by permutations.
      integer, intent(in), optional :: permutations(:, :) !! P columns of permutations of 1:S; absent means one identity.
      type(predictive_metric_result) :: out
      out = score_core(x, x2, y, .false., permutations, weights, joint_weights)
   end function loo_crps

   pure function loo_scrps(x, x2, y, weights, joint_weights, permutations) result(out)
      real(dp), intent(in) :: x(:, :) !! Predictive draws, shaped S draws by N observations.
      real(dp), intent(in) :: x2(:, :) !! Independent predictive draws of the same shape as x.
      real(dp), intent(in) :: y(:) !! Observed responses, length N.
      real(dp), intent(in) :: weights(:, :) !! Nonnegative normalized ordinary marginal LOO weights, shape S by N.
      real(dp), intent(in) :: joint_weights(:, :, :) !! Normalized paired LOO weights, shape S by N by permutations.
      integer, intent(in), optional :: permutations(:, :) !! P columns of permutations of 1:S; absent means one identity.
      type(predictive_metric_result) :: out
      out = score_core(x, x2, y, .true., permutations, weights, joint_weights)
   end function loo_scrps

   pure function loo_predictive_metric(x, y, weights, metric) result(out)
      real(dp), intent(in) :: x(:, :) !! Finite predictive draws, shaped S draws by N observations.
      real(dp), intent(in) :: y(:) !! Observed responses, length N; binary 0 or 1 for classification.
      real(dp), intent(in) :: weights(:, :) !! Normalized nonnegative ordinary LOO weights, shaped S by N.
      character(len = *), intent(in), optional :: metric !! mae (default), mse, rmse, acc, or balanced_acc.
      type(predictive_metric_result) :: out
      real(dp) :: prediction(size(y)), tn, tp, n
      integer :: j, nneg, npos
      character(len = :), allocatable :: selected

      call initialize_result(out, size(y))
      if (size(x, 1) == 0 .or. size(y) == 0) return
      if (size(x, 2) /= size(y)) return
      if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) return
      if (any(shape(weights) /= shape(x))) return
      if (.not. normalized_columns(weights)) return
      do j = 1, size(y)
         prediction(j) = sum(weights(:, j) * x(:, j)) / sum(weights(:, j))
      end do
      selected = 'mae'
      if (present(metric)) selected = metric
      select case (selected)
      case ('mae')
         out%pointwise = abs(y - prediction)
      case ('mse', 'rmse')
         out%pointwise = (y - prediction) ** 2
      case ('acc', 'balanced_acc')
         if (any(y /= 0.0_dp .and. y /= 1.0_dp)) return
         if (any(prediction < 0.0_dp .or. prediction > 1.0_dp)) return
         out%pointwise = merge(1.0_dp, 0.0_dp, merge(1.0_dp, 0.0_dp, prediction > 0.5_dp) == y)
      case default
         return
      end select
      call summarize_result(out)
      n = real(size(y), dp)
      select case (selected)
      case ('rmse')
         if (out%estimate > 0.0_dp) then
            out%se = out%se / (2.0_dp * sqrt(out%estimate))
         else
            out%se = ieee_value(0.0_dp, ieee_quiet_nan)
         end if
         out%estimate = sqrt(out%estimate)
      case ('acc')
         out%se = sqrt(out%estimate * (1.0_dp - out%estimate) / n)
      case ('balanced_acc')
         nneg = count(y == 0.0_dp)
         npos = size(y) - nneg
         if (min(nneg, npos) == 0) then
            out%estimate = ieee_value(0.0_dp, ieee_quiet_nan)
            out%se = out%estimate
            out%status = 2
            return
         end if
         tn = sum(out%pointwise, mask = y == 0.0_dp) / real(nneg, dp)
         tp = sum(out%pointwise, mask = y == 1.0_dp) / real(npos, dp)
         out%estimate = 0.5_dp * (tp + tn)
         out%se = sqrt((tp * (1.0_dp - tp) + tn * (1.0_dp - tn)) / (4.0_dp * n))
      end select
   end function loo_predictive_metric

   pure function score_core(x, x2, y, scaled, permutations, weights, joint_weights) result(out)
      real(dp), intent(in) :: x(:, :) !! Finite predictive draws, shaped draws by observations.
      real(dp), intent(in) :: x2(:, :) !! Independent draws, with the same shape as x.
      real(dp), intent(in) :: y(:) !! Responses, one per predictive-draw column.
      logical, intent(in) :: scaled !! True selects the scaled CRPS definition.
      integer, intent(in), optional :: permutations(:, :) !! Explicit x2 row permutations; absent means identity.
      real(dp), intent(in), optional :: weights(:, :) !! Normalized marginal weights, supplied together with joint_weights.
      real(dp), intent(in), optional :: joint_weights(:, :, :) !! Normalized joint weights for every permutation.
      type(predictive_metric_result) :: out
      real(dp) :: exx(size(y)), exy(size(y)), distance(size(x, 1))
      integer :: s, n, p, i, j, repeat, row

      call initialize_result(out, size(y))
      s = size(x, 1)
      n = size(y)
      if (s == 0 .or. n == 0) return
      if (any(shape(x2) /= shape(x)) .or. size(x, 2) /= n) return
      if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(x2))) return
      if (.not. all(ieee_is_finite(y))) return
      p = 1
      if (present(permutations)) then
         if (size(permutations, 1) /= s) return
         p = size(permutations, 2)
         if (p == 0) return
         if (any(permutations < 1) .or. any(permutations > s)) return
         do repeat = 1, p
            do i = 1, s
               if (count(permutations(:, repeat) == i) /= 1) return
            end do
         end do
      end if
      if (present(weights) .neqv. present(joint_weights)) return
      if (present(weights) .and. present(joint_weights)) then
         if (any(shape(weights) /= shape(x))) return
         if (any(shape(joint_weights) /= [s, n, p])) return
         if (.not. normalized_columns(weights)) return
         do repeat = 1, p
            if (.not. normalized_columns(joint_weights(:, :, repeat))) return
         end do
      end if
      exx = 0.0_dp
      do j = 1, n
         distance = abs(x(:, j) - y(j))
         if (present(weights)) then
            exy(j) = sum(weights(:, j) * distance) / sum(weights(:, j))
         else
            exy(j) = r_mean(distance)
         end if
         do repeat = 1, p
            do i = 1, s
               row = i
               if (present(permutations)) row = permutations(i, repeat)
               distance(i) = abs(x(i, j) - x2(row, j))
            end do
            if (present(joint_weights)) then
               exx(j) = exx(j) + sum(joint_weights(:, j, repeat) * distance) &
                  / (sum(joint_weights(:, j, repeat)) * real(p, dp))
            else
               exx(j) = exx(j) + r_mean(distance) / real(p, dp)
            end if
         end do
      end do
      out%pointwise = 0.5_dp * exx - exy
      if (scaled) then
         do j = 1, n
            if (exx(j) > 0.0_dp) then
               out%pointwise(j) = -exy(j) / exx(j) - 0.5_dp * log(exx(j))
            else
               ! Upstream has either -Inf + Inf or 0/0, both yielding NaN.
               out%pointwise(j) = ieee_value(0.0_dp, ieee_quiet_nan)
            end if
         end do
      end if
      call summarize_result(out)
      if (scaled .and. any(exx == 0.0_dp)) out%status = 2
   end function score_core

   pure logical function normalized_columns(weights) result(valid)
      real(dp), intent(in) :: weights(:, :) !! Ordinary weights, each column finite, nonnegative, and summing to one.
      valid = .false.
      if (size(weights, 1) == 0) return
      if (.not. all(ieee_is_finite(weights))) return
      if (any(weights < 0.0_dp)) return
      valid = all(abs(sum(weights, dim = 1) - 1.0_dp) <= 1.0e-8_dp)
   end function normalized_columns

   pure subroutine initialize_result(out, n)
      type(predictive_metric_result), intent(out) :: out !! Invalid result, initialized to quiet NaNs and status 1.
      integer, intent(in) :: n !! Number of pointwise entries to allocate.
      out%estimate = ieee_value(0.0_dp, ieee_quiet_nan)
      out%se = out%estimate
      out%status = 1
      allocate(out%pointwise(n))
      out%pointwise = out%estimate
   end subroutine initialize_result

   pure subroutine summarize_result(out)
      type(predictive_metric_result), intent(inout) :: out !! Pointwise values to summarize by their mean and standard error.
      out%estimate = r_mean(out%pointwise)
      out%se = ieee_value(0.0_dp, ieee_quiet_nan)
      if (size(out%pointwise) > 1 .and. all(ieee_is_finite(out%pointwise))) then
         out%se = r_sd(out%pointwise) / sqrt(real(size(out%pointwise), dp))
      end if
      out%status = 0
   end subroutine summarize_result
end module loo_metrics
