! SPDX-License-Identifier: GPL-3.0-or-later
! Numerical translation of loo 2.10.1, R/loo_model_weights.R.
! See package NOTICE and upstream license for attribution.
module loo_model_weights
   use r_kinds, only: dp
   use r_stability, only: r_log_sum_exp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_quiet_nan, ieee_value
   implicit none
   private
   public :: stacking_weights, pseudobma_weights, model_weights
contains
   ! Maximize the concave log score on the closed probability simplex by
   ! pairwise vertex exchange, using an exact one-dimensional line search.
   ! A duality gap certifies convergence; boundary solutions are permitted.
   pure subroutine stacking_weights(lpd_point, weights, status, objective, iterations, tolerance, max_iterations, kkt_gap)
      real(dp), intent(in) :: lpd_point(:, :) !! N observations by K models of pointwise log predictive densities.
      real(dp), intent(out) :: weights(:) !! K nonnegative optimal model weights, summing to one; NaNs on invalid input.
      integer, intent(out), optional :: status !! 0 converged, 1 invalid input, 2 iteration limit or numerical stall.
      real(dp), intent(out), optional :: objective !! Maximized total log predictive density, also returned if unconverged.
      integer, intent(out), optional :: iterations !! Number of completed vertex-exchange iterations.
      real(dp), intent(in), optional :: tolerance !! Positive per-observation duality-gap threshold; default 1e-10.
      integer, intent(in), optional :: max_iterations !! Positive iteration limit; default 10000.
      real(dp), intent(out), optional :: kkt_gap !! Nonnegative final per-observation primal-dual gap.
      real(dp) :: density(size(lpd_point, 1), size(lpd_point, 2)), row_max(size(lpd_point, 1))
      real(dp) :: mixture(size(lpd_point, 1)), direction(size(lpd_point, 1)), gradient(size(lpd_point, 2))
      real(dp) :: tol, gap, lower, upper, alpha, derivative, mid, previous, nreal
      integer :: n, k, limit, iter, best, worst, j, bisect, code, performed

      weights = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = 1
      if (present(objective)) objective = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(iterations)) iterations = 0
      if (present(kkt_gap)) kkt_gap = ieee_value(0.0_dp, ieee_quiet_nan)
      n = size(lpd_point, 1)
      k = size(lpd_point, 2)
      if (size(weights) /= k .or. n == 0 .or. k < 2) return
      if (.not. valid_lpd(lpd_point)) return
      tol = 1.0e-10_dp
      if (present(tolerance)) tol = tolerance
      if (.not. ieee_is_finite(tol) .or. tol <= 0.0_dp) return
      limit = 10000
      if (present(max_iterations)) limit = max_iterations
      if (limit <= 0) return
      nreal = real(n, dp)
      row_max = maxval(lpd_point, dim = 2)
      do j = 1, k
         density(:, j) = exp(lpd_point(:, j) - row_max)
      end do
      weights = 1.0_dp / real(k, dp)
      code = 2
      performed = 0
      do iter = 1, limit
         mixture = matmul(density, weights)
         do j = 1, k
            gradient(j) = sum(density(:, j) / mixture)
         end do
         gap = max(0.0_dp, (maxval(gradient) - dot_product(weights, gradient)) / nreal)
         if (gap <= tol) then
            code = 0
            exit
         end if
         best = maxloc(gradient, dim = 1)
         worst = minloc(gradient, dim = 1, mask = weights > 0.0_dp)
         if (best == worst) exit
         direction = density(:, best) - density(:, worst)
         upper = weights(worst)
         derivative = exchange_derivative(mixture, direction, upper)
         if (derivative >= 0.0_dp) then
            alpha = upper
         else
            lower = 0.0_dp
            do bisect = 1, 70
               mid = lower + 0.5_dp * (upper - lower)
               derivative = exchange_derivative(mixture, direction, mid)
               if (derivative > 0.0_dp) then
                  lower = mid
               else
                  upper = mid
               end if
            end do
            alpha = lower + 0.5_dp * (upper - lower)
         end if
         previous = weights(best)
         weights(best) = weights(best) + alpha
         weights(worst) = max(0.0_dp, weights(worst) - alpha)
         performed = iter
         if (weights(best) == previous) exit
      end do
      ! Recompute the certificate even when the final permitted step converged.
      weights = weights / sum(weights)
      mixture = matmul(density, weights)
      do j = 1, k
         gradient(j) = sum(density(:, j) / mixture)
      end do
      gap = max(0.0_dp, (maxval(gradient) - dot_product(weights, gradient)) / nreal)
      if (gap <= tol) code = 0
      if (present(status)) status = code
      if (present(iterations)) iterations = performed
      if (present(objective)) objective = sum(row_max + log(mixture))
      if (present(kkt_gap)) kkt_gap = gap
   end subroutine stacking_weights

   ! With no bootstrap matrix this computes plain pseudo-BMA. Supplying B by N
   ! normalized Dirichlet draws computes pseudo-BMA+ without imposing an RNG.
   pure subroutine pseudobma_weights(lpd_point, weights, bootstrap_weights, status)
      real(dp), intent(in) :: lpd_point(:, :) !! N observations by K models of pointwise log predictive densities.
      real(dp), intent(out) :: weights(:) !! K model weights summing to one, or quiet NaNs for invalid arguments.
      real(dp), intent(in), optional :: bootstrap_weights(:, :) !! B by N ordinary Bayesian-bootstrap weights; rows sum to one.
      integer, intent(out), optional :: status !! 0 success, 1 invalid input or no model with finite total score.
      real(dp) :: shifted(size(lpd_point, 1), size(lpd_point, 2)), totals(size(lpd_point, 2))
      real(dp) :: replicate_weights(size(lpd_point, 2)), row_max(size(lpd_point, 1)), norm
      integer :: n, k, b, nb, i, j

      weights = ieee_value(0.0_dp, ieee_quiet_nan)
      if (present(status)) status = 1
      n = size(lpd_point, 1)
      k = size(lpd_point, 2)
      if (size(weights) /= k .or. n == 0 .or. k < 2) return
      if (.not. valid_lpd(lpd_point)) return
      row_max = maxval(lpd_point, dim = 2)
      do j = 1, k
         shifted(:, j) = lpd_point(:, j) - row_max
      end do
      if (.not. present(bootstrap_weights)) then
         totals = sum(shifted, dim = 1)
         if (.not. any(ieee_is_finite(totals))) return
         norm = r_log_sum_exp(totals)
         if (.not. ieee_is_finite(norm)) return
         weights = exp(totals - norm)
      else
         nb = size(bootstrap_weights, 1)
         if (nb == 0 .or. size(bootstrap_weights, 2) /= n) return
         if (.not. all(ieee_is_finite(bootstrap_weights))) return
         if (any(bootstrap_weights < 0.0_dp)) return
         if (any(abs(sum(bootstrap_weights, dim = 2) - 1.0_dp) > 1.0e-8_dp)) return
         weights = 0.0_dp
         do b = 1, nb
            totals = 0.0_dp
            do i = 1, n
               ! Avoid zero times negative infinity for sparse bootstrap rows.
               if (bootstrap_weights(b, i) > 0.0_dp) then
                  totals = totals + bootstrap_weights(b, i) * shifted(i, :) * real(n, dp)
               end if
            end do
            if (.not. any(ieee_is_finite(totals))) then
               weights = ieee_value(0.0_dp, ieee_quiet_nan)
               return
            end if
            norm = r_log_sum_exp(totals)
            if (.not. ieee_is_finite(norm)) then
               weights = ieee_value(0.0_dp, ieee_quiet_nan)
               return
            end if
            replicate_weights = exp(totals - norm)
            weights = weights + replicate_weights / real(nb, dp)
         end do
      end if
      weights = weights / sum(weights)
      if (present(status)) status = 0
   end subroutine pseudobma_weights

   ! Numerical loo_model_weights wrapper, accepting already computed pointwise
   ! LOO/K-fold scores instead of R lists of likelihood matrices or S3 objects.
   pure subroutine model_weights(lpd_point, weights, method, bootstrap_weights, status, tolerance, max_iterations)
      real(dp), intent(in) :: lpd_point(:, :) !! Pointwise LOO/K-fold log scores, N observations by K models.
      real(dp), intent(out) :: weights(:) !! Model probabilities, length K, or quiet NaNs on invalid arguments.
      character(len = *), intent(in), optional :: method !! stacking (default), pseudobma, or pseudobma+.
      real(dp), intent(in), optional :: bootstrap_weights(:, :) !! B by N normalized bootstrap rows, required for pseudobma+.
      integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 stacking failed to converge.
      real(dp), intent(in), optional :: tolerance !! Stacking per-observation duality-gap tolerance; default 1e-10.
      integer, intent(in), optional :: max_iterations !! Stacking maximum iterations; default 10000.
      character(len = :), allocatable :: selected

      selected = 'stacking'
      if (present(method)) selected = method
      select case (selected)
      case ('stacking')
         call stacking_weights(lpd_point, weights, status, tolerance = tolerance, max_iterations = max_iterations)
      case ('pseudobma', 'pseudobma+')
         if (selected == 'pseudobma+' .and. .not. present(bootstrap_weights)) then
            weights = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(status)) status = 1
            return
         end if
         call pseudobma_weights(lpd_point, weights, bootstrap_weights, status)
      case default
         weights = ieee_value(0.0_dp, ieee_quiet_nan)
         if (present(status)) status = 1
      end select
   end subroutine model_weights

   pure real(dp) function exchange_derivative(mixture, direction, step) result(derivative)
      real(dp), intent(in) :: mixture(:) !! Current positive mixture density after rowwise scaling.
      real(dp), intent(in) :: direction(:) !! Per-observation difference of entering and leaving model densities.
      real(dp), intent(in) :: step !! Nonnegative probability transferred from the leaving model.
      real(dp) :: candidate(size(mixture))
      candidate = mixture + step * direction
      if (any(candidate <= 0.0_dp)) then
         derivative = -huge(1.0_dp)
      else
         derivative = sum(direction / candidate)
      end if
   end function exchange_derivative

   pure logical function valid_lpd(lpd_point) result(valid)
      real(dp), intent(in) :: lpd_point(:, :) !! Log densities permitting negative infinity but not NaN or positive infinity.
      valid = .false.
      if (any(ieee_is_nan(lpd_point))) return
      if (any(.not. ieee_is_finite(lpd_point) .and. lpd_point > 0.0_dp)) return
      if (any(count(ieee_is_finite(lpd_point), dim = 2) == 0)) return
      valid = .true.
   end function valid_lpd
end module loo_model_weights
