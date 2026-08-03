! SPDX-License-Identifier: GPL-2.0-only
module glmnet_gaussian
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use glmnet_kinds, only : dp, glmnet_eps, glmnet_huge
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument, &
      glmnet_nonfinite_input, glmnet_max_iterations, glmnet_invalid_response
   use glmnet_types, only : glmnet_control_type, glmnet_path_result, &
      glmnet_family_gaussian, glmnet_family_mgaussian, family_name
   use glmnet_utils, only : normalize_weights, prepare_design, make_lambda_sequence, &
      soft_threshold, group_soft_threshold, clamp_value, count_nonzero_rows, &
      all_finite_vector, all_finite_matrix, weighted_mean
   implicit none
   private
   public :: fit_gaussian_path, fit_mgaussian_path, solve_weighted_elastic_net
contains
   subroutine fit_gaussian_path(x, y, result, control, weights_in, offset_in, &
      lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y(:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      type(glmnet_control_type) :: ctl
      real(dp), allocatable :: weights(:), offset(:), xw(:,:), xm(:), xs(:)
      real(dp), allocatable :: penalty(:), lower(:), upper(:), lambda(:), yw(:)
      real(dp), allocatable :: beta_std(:), residual(:), fitted(:)
      logical, allocatable :: usable(:), excluded(:)
      real(dp) :: alpha_sequence, lambda_max, null_mean, null_rss, rss
      real(dp) :: old_objective, current_objective, fractional
      integer :: n, p, l, status, iter, passes
      logical :: converged

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      p = size(x, 2)
      call initialize_empty_result(result, glmnet_family_gaussian, n, p, 1)
      if (n < 2 .or. p < 1 .or. size(y) /= n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (.not. all_finite_vector(y)) then
         result%status = glmnet_nonfinite_input
         return
      end if
      call normalize_weights(n, weights_in, weights, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      allocate(offset(n))
      if (present(offset_in)) then
         if (size(offset_in) /= n .or. .not. all_finite_vector(offset_in)) then
            result%status = glmnet_invalid_argument
            return
         end if
         offset = offset_in
      else
         offset = 0.0_dp
      end if
      call prepare_design(x, weights, ctl, xw, xm, xs, usable, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call initialize_constraints(p, penalty_factor_in, lower_in, upper_in, &
         excluded_in, xs, usable, penalty, lower, upper, excluded, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      allocate(yw(n), beta_std(p), residual(n), fitted(n))
      yw = y - offset
      if (ctl%intercept) then
         null_mean = weighted_mean(yw, weights)
      else
         null_mean = 0.0_dp
      end if
      null_rss = sum(weights * (yw - null_mean) ** 2)
      if (null_rss <= 100.0_dp * glmnet_eps) then
         result%status = glmnet_invalid_response
         return
      end if
      alpha_sequence = max(ctl%alpha, 1.0e-3_dp)
      lambda_max = 0.0_dp
      do l = 1, p
         if (.not. excluded(l)) then
            lambda_max = max(lambda_max, abs(sum(weights * xw(:, l) * &
               (yw - null_mean))) / max(alpha_sequence * penalty(l), glmnet_eps))
         end if
      end do
      call make_lambda_sequence(lambda_max, n, p, ctl, lambda_in, lambda, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call allocate_path_arrays(result, size(lambda), p, 1)
      result%lambda = lambda
      result%nulldev = null_rss
      result%x_mean = xm
      result%x_scale = xs
      result%standardize = ctl%standardize
      result%intercept_fitted = ctl%intercept
      beta_std = 0.0_dp
      fitted = null_mean
      old_objective = huge(1.0_dp)
      do l = 1, size(lambda)
         call solve_weighted_elastic_net(xw, yw, weights, lambda(l), ctl%alpha, &
            penalty, lower, upper, excluded, ctl%intercept, beta_std, &
            result%intercept(1, l), iter, passes, converged, ctl%threshold, &
            ctl%max_iterations)
         result%npasses = result%npasses + passes
         fitted = result%intercept(1, l) + matmul(xw, beta_std)
         residual = yw - fitted
         rss = sum(weights * residual ** 2)
         current_objective = 0.5_dp * rss + elastic_net_penalty(beta_std, &
            penalty, lambda(l), ctl%alpha)
         result%objective(l) = current_objective
         result%dev_ratio(l) = max(0.0_dp, min(1.0_dp, 1.0_dp - rss / null_rss))
         result%iterations(l) = iter
         result%converged(l) = converged
         result%beta(:, 1, l) = beta_std / xs
         result%intercept(1, l) = result%intercept(1, l) - &
            dot_product(xm, result%beta(:, 1, l))
         result%df(l) = count(abs(result%beta(:, 1, l)) > 1.0e-12_dp)
         if (.not. converged .and. result%status == glmnet_success) &
            result%status = glmnet_max_iterations
         if (.not. present(lambda_in) .and. l >= ctl%minimum_lambda_count) then
            fractional = abs(old_objective - current_objective) / &
               max(abs(old_objective), glmnet_eps)
            if (result%dev_ratio(l) >= ctl%deviance_max) then
               call truncate_path(result, l)
               exit
            end if
            if (fractional < ctl%fractional_deviance .and. &
                result%dev_ratio(l) > 0.0_dp) then
               call truncate_path(result, l)
               exit
            end if
         end if
         old_objective = current_objective
      end do
   end subroutine fit_gaussian_path

   subroutine solve_weighted_elastic_net(x, y, weights, lambda, alpha, penalty, &
      lower, upper, excluded, fit_intercept, beta, intercept, iterations, passes, &
      converged, threshold, max_iterations)
      real(dp), intent(in) :: x(:,:), y(:), weights(:), lambda, alpha
      real(dp), intent(in) :: penalty(:), lower(:), upper(:)
      logical, intent(in) :: excluded(:), fit_intercept
      real(dp), intent(inout) :: beta(:)
      real(dp), intent(out) :: intercept
      integer, intent(out) :: iterations, passes
      logical, intent(out) :: converged
      real(dp), intent(in) :: threshold
      integer, intent(in) :: max_iterations
      real(dp), allocatable :: residual(:)
      real(dp) :: old_beta, new_beta, partial, denominator, change, max_change
      real(dp) :: delta_intercept
      integer :: j, p
      p = size(x, 2)
      allocate(residual(size(y)))
      if (fit_intercept) then
         intercept = sum(weights * (y - matmul(x, beta))) / max(sum(weights), glmnet_eps)
      else
         intercept = 0.0_dp
      end if
      residual = y - intercept - matmul(x, beta)
      converged = .false.
      passes = 0
      do iterations = 1, max_iterations
         max_change = 0.0_dp
         if (fit_intercept) then
            delta_intercept = sum(weights * residual) / max(sum(weights), glmnet_eps)
            intercept = intercept + delta_intercept
            residual = residual - delta_intercept
            max_change = max(max_change, abs(delta_intercept))
         end if
         do j = 1, p
            if (excluded(j)) then
               if (abs(beta(j)) > 0.0_dp) then
                  residual = residual + x(:, j) * beta(j)
                  beta(j) = 0.0_dp
               end if
               cycle
            end if
            old_beta = beta(j)
            partial = sum(weights * x(:, j) * (residual + x(:, j) * old_beta))
            denominator = sum(weights * x(:, j) ** 2) + &
               lambda * (1.0_dp - alpha) * penalty(j)
            new_beta = soft_threshold(partial, lambda * alpha * penalty(j)) / &
               max(denominator, glmnet_eps)
            new_beta = clamp_value(new_beta, lower(j), upper(j))
            change = new_beta - old_beta
            if (abs(change) > 0.0_dp) then
               residual = residual - x(:, j) * change
               beta(j) = new_beta
               max_change = max(max_change, abs(change))
            end if
         end do
         passes = passes + 1
         if (max_change <= threshold * (1.0_dp + maxval(abs(beta)))) then
            converged = .true.
            exit
         end if
      end do
      if (.not. converged) iterations = max_iterations
   end subroutine solve_weighted_elastic_net

   subroutine fit_mgaussian_path(x, y, result, control, weights_in, offset_in, &
      lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y(:,:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:,:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      type(glmnet_control_type) :: ctl
      real(dp), allocatable :: weights(:), offset(:,:), xw(:,:), xm(:), xs(:)
      real(dp), allocatable :: penalty(:), lower(:), upper(:), lambda(:), yw(:,:)
      real(dp), allocatable :: beta_std(:,:), intercept_std(:), residual(:,:)
      logical, allocatable :: usable(:), excluded(:)
      real(dp) :: lambda_max, alpha_sequence, nulldev, rss, norm_gradient
      real(dp) :: current_objective, old_objective, fractional
      integer :: n, p, q, l, j, status, iter, passes
      logical :: converged

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      p = size(x, 2)
      q = size(y, 2)
      call initialize_empty_result(result, glmnet_family_mgaussian, n, p, q)
      if (n < 2 .or. p < 1 .or. q < 1 .or. size(y, 1) /= n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (.not. all_finite_matrix(y)) then
         result%status = glmnet_nonfinite_input
         return
      end if
      call normalize_weights(n, weights_in, weights, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      allocate(offset(n, q))
      if (present(offset_in)) then
         if (size(offset_in, 1) /= n .or. size(offset_in, 2) /= q .or. &
             .not. all_finite_matrix(offset_in)) then
            result%status = glmnet_invalid_argument
            return
         end if
         offset = offset_in
      else
         offset = 0.0_dp
      end if
      call prepare_design(x, weights, ctl, xw, xm, xs, usable, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call initialize_constraints(p, penalty_factor_in, lower_in, upper_in, &
         excluded_in, xs, usable, penalty, lower, upper, excluded, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      allocate(yw(n, q), beta_std(p, q), intercept_std(q), residual(n, q))
      yw = y - offset
      if (ctl%intercept) then
         do j = 1, q
            intercept_std(j) = weighted_mean(yw(:, j), weights)
         end do
      else
         intercept_std = 0.0_dp
      end if
      residual = yw - spread(intercept_std, 1, n)
      nulldev = sum(spread(weights, 2, q) * residual ** 2)
      if (nulldev <= 100.0_dp * glmnet_eps) then
         result%status = glmnet_invalid_response
         return
      end if
      alpha_sequence = max(ctl%alpha, 1.0e-3_dp)
      lambda_max = 0.0_dp
      do j = 1, p
         if (.not. excluded(j)) then
            norm_gradient = sqrt(sum([(sum(weights * xw(:, j) * residual(:, l)) ** 2, l=1,q)]))
            lambda_max = max(lambda_max, norm_gradient / &
               max(alpha_sequence * penalty(j), glmnet_eps))
         end if
      end do
      call make_lambda_sequence(lambda_max, n, p, ctl, lambda_in, lambda, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call allocate_path_arrays(result, size(lambda), p, q)
      result%lambda = lambda
      result%nulldev = nulldev
      result%x_mean = xm
      result%x_scale = xs
      result%standardize = ctl%standardize
      result%intercept_fitted = ctl%intercept
      beta_std = 0.0_dp
      old_objective = huge(1.0_dp)
      do l = 1, size(lambda)
         call solve_group_elastic_net(xw, yw, weights, lambda(l), ctl%alpha, &
            penalty, lower, upper, excluded, ctl%intercept, beta_std, &
            intercept_std, iter, passes, converged, ctl%threshold, ctl%max_iterations)
         result%npasses = result%npasses + passes
         residual = yw - spread(intercept_std, 1, n) - matmul(xw, beta_std)
         rss = sum(spread(weights, 2, q) * residual ** 2)
         current_objective = 0.5_dp * rss + group_elastic_net_penalty(beta_std, &
            penalty, lambda(l), ctl%alpha)
         result%objective(l) = current_objective
         result%dev_ratio(l) = max(0.0_dp, min(1.0_dp, 1.0_dp - rss / nulldev))
         result%iterations(l) = iter
         result%converged(l) = converged
         do j = 1, q
            result%beta(:, j, l) = beta_std(:, j) / xs
            result%intercept(j, l) = intercept_std(j) - dot_product(xm, result%beta(:, j, l))
         end do
         result%df(l) = count_nonzero_rows(result%beta(:, :, l))
         if (.not. converged .and. result%status == glmnet_success) &
            result%status = glmnet_max_iterations
         if (.not. present(lambda_in) .and. l >= ctl%minimum_lambda_count) then
            fractional = abs(old_objective - current_objective) / &
               max(abs(old_objective), glmnet_eps)
            if (result%dev_ratio(l) >= ctl%deviance_max .or. &
                (fractional < ctl%fractional_deviance .and. result%dev_ratio(l) > 0.0_dp)) then
               call truncate_path(result, l)
               exit
            end if
         end if
         old_objective = current_objective
      end do
   end subroutine fit_mgaussian_path

   subroutine solve_group_elastic_net(x, y, weights, lambda, alpha, penalty, lower, &
      upper, excluded, fit_intercept, beta, intercept, iterations, passes, converged, &
      threshold, max_iterations)
      real(dp), intent(in) :: x(:,:), y(:,:), weights(:), lambda, alpha
      real(dp), intent(in) :: penalty(:), lower(:), upper(:)
      logical, intent(in) :: excluded(:), fit_intercept
      real(dp), intent(inout) :: beta(:,:), intercept(:)
      integer, intent(out) :: iterations, passes
      logical, intent(out) :: converged
      real(dp), intent(in) :: threshold
      integer, intent(in) :: max_iterations
      real(dp), allocatable :: residual(:,:), partial(:), candidate(:), delta(:)
      real(dp) :: denominator, max_change
      integer :: j, k, n, p, q
      n = size(x, 1)
      p = size(x, 2)
      q = size(y, 2)
      allocate(residual(n, q), partial(q), candidate(q), delta(q))
      residual = y - spread(intercept, 1, n) - matmul(x, beta)
      converged = .false.
      passes = 0
      do iterations = 1, max_iterations
         max_change = 0.0_dp
         if (fit_intercept) then
            do k = 1, q
               delta(k) = sum(weights * residual(:, k)) / max(sum(weights), glmnet_eps)
            end do
            intercept = intercept + delta
            residual = residual - spread(delta, 1, n)
            max_change = max(max_change, maxval(abs(delta)))
         end if
         do j = 1, p
            if (excluded(j)) then
               if (maxval(abs(beta(j, :))) > 0.0_dp) then
                  residual = residual + matmul(reshape(x(:, j), [n, 1]), reshape(beta(j, :), [1, q]))
                  beta(j, :) = 0.0_dp
               end if
               cycle
            end if
            do k = 1, q
               partial(k) = sum(weights * x(:, j) * &
                  (residual(:, k) + x(:, j) * beta(j, k)))
            end do
            denominator = sum(weights * x(:, j) ** 2) + &
               lambda * (1.0_dp - alpha) * penalty(j)
            call group_soft_threshold(partial, lambda * alpha * penalty(j), candidate)
            candidate = candidate / max(denominator, glmnet_eps)
            candidate = min(max(candidate, lower(j)), upper(j))
            delta = candidate - beta(j, :)
            if (maxval(abs(delta)) > 0.0_dp) then
               residual = residual - matmul(reshape(x(:, j), [n, 1]), reshape(delta, [1, q]))
               beta(j, :) = candidate
               max_change = max(max_change, maxval(abs(delta)))
            end if
         end do
         passes = passes + 1
         if (max_change <= threshold * (1.0_dp + maxval(abs(beta)))) then
            converged = .true.
            exit
         end if
      end do
      if (.not. converged) iterations = max_iterations
   end subroutine solve_group_elastic_net

   subroutine initialize_constraints(p, penalty_in, lower_in, upper_in, excluded_in, &
      scale, usable, penalty, lower, upper, excluded, status)
      integer, intent(in) :: p
      real(dp), intent(in), optional :: penalty_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      real(dp), intent(in) :: scale(:)
      logical, intent(in) :: usable(:)
      real(dp), allocatable, intent(out) :: penalty(:), lower(:), upper(:)
      logical, allocatable, intent(out) :: excluded(:)
      integer, intent(out) :: status
      status = glmnet_success
      allocate(penalty(p), lower(p), upper(p), excluded(p))
      penalty = 1.0_dp
      lower = -glmnet_huge
      upper = glmnet_huge
      excluded = .not. usable
      if (present(penalty_in)) then
         if (size(penalty_in) /= p .or. .not. all_finite_vector(penalty_in) .or. &
             any(penalty_in < 0.0_dp)) then
            status = glmnet_invalid_argument
            return
         end if
         penalty = penalty_in
      end if
      if (present(lower_in)) then
         if (size(lower_in) /= p .or. .not. all_finite_vector(lower_in)) then
            status = glmnet_invalid_argument
            return
         end if
         lower = lower_in * scale
      end if
      if (present(upper_in)) then
         if (size(upper_in) /= p .or. .not. all_finite_vector(upper_in)) then
            status = glmnet_invalid_argument
            return
         end if
         upper = upper_in * scale
      end if
      if (any(lower > upper)) then
         status = glmnet_invalid_argument
         return
      end if
      if (present(excluded_in)) then
         if (size(excluded_in) /= p) then
            status = glmnet_invalid_argument
            return
         end if
         excluded = excluded .or. excluded_in
      end if
   end subroutine initialize_constraints

   pure function elastic_net_penalty(beta, penalty, lambda, alpha) result(value)
      real(dp), intent(in) :: beta(:), penalty(:), lambda, alpha
      real(dp) :: value
      value = lambda * sum(penalty * (alpha * abs(beta) + &
         0.5_dp * (1.0_dp - alpha) * beta ** 2))
   end function elastic_net_penalty

   pure function group_elastic_net_penalty(beta, penalty, lambda, alpha) result(value)
      real(dp), intent(in) :: beta(:,:), penalty(:), lambda, alpha
      real(dp) :: value
      integer :: j
      value = 0.0_dp
      do j = 1, size(beta, 1)
         value = value + lambda * penalty(j) * (alpha * sqrt(sum(beta(j, :) ** 2)) + &
            0.5_dp * (1.0_dp - alpha) * sum(beta(j, :) ** 2))
      end do
   end function group_elastic_net_penalty

   subroutine initialize_empty_result(result, family_code, n, p, q)
      type(glmnet_path_result), intent(out) :: result
      integer, intent(in) :: family_code, n, p, q
      result%family_code = family_code
      result%family = family_name(family_code)
      result%status = glmnet_success
      result%nobs = n
      result%nvars = p
      result%nout = q
      result%nlambda = 0
      result%npasses = 0
   end subroutine initialize_empty_result

   subroutine allocate_path_arrays(result, nlambda, p, q)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: nlambda, p, q
      result%nlambda = nlambda
      allocate(result%lambda(nlambda), result%intercept(q, nlambda), &
         result%beta(p, q, nlambda), result%dev_ratio(nlambda), &
         result%objective(nlambda), result%df(nlambda), result%iterations(nlambda), &
         result%converged(nlambda))
      result%intercept = 0.0_dp
      result%beta = 0.0_dp
      result%dev_ratio = 0.0_dp
      result%objective = 0.0_dp
      result%df = 0
      result%iterations = 0
      result%converged = .false.
   end subroutine allocate_path_arrays

   subroutine truncate_path(result, keep)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: keep
      real(dp), allocatable :: lambda(:), intercept(:,:), beta(:,:,:), dev(:), objective(:)
      integer, allocatable :: df(:), iterations(:)
      logical, allocatable :: converged(:)
      if (keep >= result%nlambda) return
      allocate(lambda(keep), intercept(result%nout, keep), &
         beta(result%nvars, result%nout, keep), dev(keep), objective(keep), &
         df(keep), iterations(keep), converged(keep))
      lambda = result%lambda(:keep)
      intercept = result%intercept(:, :keep)
      beta = result%beta(:, :, :keep)
      dev = result%dev_ratio(:keep)
      objective = result%objective(:keep)
      df = result%df(:keep)
      iterations = result%iterations(:keep)
      converged = result%converged(:keep)
      call move_alloc(lambda, result%lambda)
      call move_alloc(intercept, result%intercept)
      call move_alloc(beta, result%beta)
      call move_alloc(dev, result%dev_ratio)
      call move_alloc(objective, result%objective)
      call move_alloc(df, result%df)
      call move_alloc(iterations, result%iterations)
      call move_alloc(converged, result%converged)
      result%nlambda = keep
   end subroutine truncate_path
end module glmnet_gaussian
