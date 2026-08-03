! SPDX-License-Identifier: GPL-2.0-only
module glmnet_glm
   use glmnet_kinds, only : dp, glmnet_eps, glmnet_huge
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument, &
      glmnet_nonfinite_input, glmnet_max_iterations, glmnet_invalid_response, &
      glmnet_numerical_failure
   use glmnet_types, only : glmnet_control_type, glmnet_path_result, &
      glmnet_family_binomial, glmnet_family_poisson, family_name, &
      glmnet_family_working_interface
   use glmnet_utils, only : normalize_weights, prepare_design, make_lambda_sequence, &
      all_finite_vector, weighted_mean, logistic, safe_log, safe_exp
   use glmnet_gaussian, only : solve_weighted_elastic_net
   implicit none
   private
   public :: fit_binomial_path, fit_poisson_path, fit_custom_family_path
contains
   subroutine fit_binomial_path(x, y, result, control, weights_in, offset_in, &
      lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y(:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      call fit_irls_path(x, y, glmnet_family_binomial, result, control, weights_in, &
         offset_in, lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
   end subroutine fit_binomial_path

   subroutine fit_poisson_path(x, y, result, control, weights_in, offset_in, &
      lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y(:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      call fit_irls_path(x, y, glmnet_family_poisson, result, control, weights_in, &
         offset_in, lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
   end subroutine fit_poisson_path

   subroutine fit_irls_path(x, y, family_code, result, control, weights_in, &
      offset_in, lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y(:)
      integer, intent(in) :: family_code
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      type(glmnet_control_type) :: ctl
      real(dp), allocatable :: weights(:), offset(:), xw(:,:), xm(:), xs(:)
      real(dp), allocatable :: penalty(:), lower(:), upper(:), lambda(:)
      real(dp), allocatable :: beta_std(:), eta(:), mu(:), var(:), working(:), irls_weight(:)
      logical, allocatable :: usable(:), excluded(:)
      real(dp) :: null_mean, null_intercept, intercept_std, nulldev, deviance, objective
      real(dp) :: lambda_max, alpha_sequence, old_intercept
      real(dp) :: old_objective, fractional
      integer :: n, p, j, l, outer, status, inner_iter, inner_passes
      logical :: converged, inner_converged

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      p = size(x, 2)
      call initialize_result(result, family_code, n, p)
      if (n < 2 .or. p < 1 .or. size(y) /= n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (.not. all_finite_vector(y)) then
         result%status = glmnet_nonfinite_input
         return
      end if
      if (family_code == glmnet_family_binomial) then
         if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            result%status = glmnet_invalid_response
            return
         end if
      else
         if (any(y < 0.0_dp)) then
            result%status = glmnet_invalid_response
            return
         end if
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
      allocate(beta_std(p), eta(n), mu(n), var(n), working(n), irls_weight(n))
      beta_std = 0.0_dp
      if (family_code == glmnet_family_binomial) then
         null_mean = min(max(weighted_mean(y, weights), ctl%probability_min), &
            1.0_dp - ctl%probability_min)
         if (ctl%intercept) then
            null_intercept = log(null_mean / (1.0_dp - null_mean)) - weighted_mean(offset, weights)
         else
            null_intercept = 0.0_dp
         end if
         eta = null_intercept + offset
         mu = logistic(eta)
         nulldev = binomial_deviance(y, spread(null_mean, 1, n), weights, ctl%probability_min)
      else
         null_mean = max(weighted_mean(y, weights), ctl%probability_min)
         if (ctl%intercept) then
            null_intercept = log(null_mean) - weighted_mean(offset, weights)
         else
            null_intercept = 0.0_dp
         end if
         eta = null_intercept + offset
         mu = safe_exp(eta)
         nulldev = poisson_deviance(y, spread(null_mean, 1, n), weights)
      end if
      if (nulldev <= 100.0_dp * glmnet_eps) then
         result%status = glmnet_invalid_response
         return
      end if
      alpha_sequence = max(ctl%alpha, 1.0e-3_dp)
      lambda_max = 0.0_dp
      do j = 1, p
         if (.not. excluded(j)) then
            lambda_max = max(lambda_max, abs(sum(weights * xw(:, j) * (y - mu))) / &
               max(alpha_sequence * penalty(j), glmnet_eps))
         end if
      end do
      call make_lambda_sequence(lambda_max, n, p, ctl, lambda_in, lambda, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call allocate_result(result, size(lambda), p)
      result%lambda = lambda
      result%nulldev = nulldev
      result%x_mean = xm
      result%x_scale = xs
      result%standardize = ctl%standardize
      result%intercept_fitted = ctl%intercept
      intercept_std = null_intercept
      deviance = nulldev
      objective = 0.5_dp * nulldev
      old_objective = huge(1.0_dp)
      do l = 1, size(lambda)
         converged = .false.
         do outer = 1, ctl%max_outer_iterations
            old_intercept = intercept_std
            eta = old_intercept + matmul(xw, beta_std) + offset
            if (family_code == glmnet_family_binomial) then
               mu = min(max(logistic(eta), ctl%probability_min), 1.0_dp - ctl%probability_min)
               var = max(mu * (1.0_dp - mu), ctl%probability_min)
            else
               mu = min(safe_exp(min(eta, ctl%eta_max)), exp(ctl%eta_max))
               var = max(mu, ctl%probability_min)
            end if
            irls_weight = weights * var
            working = eta + (y - mu) / var - offset
            call solve_weighted_elastic_net(xw, working, irls_weight, lambda(l), &
               ctl%alpha, penalty, lower, upper, excluded, ctl%intercept, beta_std, &
               intercept_std, inner_iter, inner_passes, inner_converged, &
               ctl%threshold, ctl%max_iterations)
            result%npasses = result%npasses + inner_passes
            eta = intercept_std + matmul(xw, beta_std) + offset
            if (family_code == glmnet_family_binomial) then
               mu = min(max(logistic(eta), ctl%probability_min), 1.0_dp - ctl%probability_min)
               deviance = binomial_deviance(y, mu, weights, ctl%probability_min)
               objective = 0.5_dp * deviance + penalty_value(beta_std, penalty, &
                  lambda(l), ctl%alpha)
            else
               mu = min(safe_exp(min(eta, ctl%eta_max)), exp(ctl%eta_max))
               deviance = poisson_deviance(y, mu, weights)
               objective = 0.5_dp * deviance + penalty_value(beta_std, penalty, &
                  lambda(l), ctl%alpha)
            end if
            if (outer > 1) then
               if (abs(result%objective(l) - objective) <= ctl%threshold * &
                   (1.0_dp + abs(objective))) then
                  converged = .true.
                  exit
               end if
            end if
            result%objective(l) = objective
         end do
         result%objective(l) = objective
         result%dev_ratio(l) = max(0.0_dp, min(1.0_dp, 1.0_dp - deviance / nulldev))
         result%iterations(l) = min(outer, ctl%max_outer_iterations)
         result%converged(l) = converged
         result%beta(:, 1, l) = beta_std / xs
         result%intercept(1, l) = intercept_std - &
            dot_product(xm, result%beta(:, 1, l))
         result%df(l) = count(abs(result%beta(:, 1, l)) > 1.0e-12_dp)
         if (.not. converged .and. result%status == glmnet_success) &
            result%status = glmnet_max_iterations
         if (.not. present(lambda_in) .and. l >= ctl%minimum_lambda_count) then
            fractional = abs(old_objective - objective) / max(abs(old_objective), glmnet_eps)
            if (result%dev_ratio(l) >= ctl%deviance_max .or. &
                (fractional < ctl%fractional_deviance .and. result%dev_ratio(l) > 0.0_dp)) then
               call truncate_result(result, l)
               exit
            end if
         end if
         old_objective = objective
         ! Restore standardized intercept for the next warm start.
         result%intercept(1, l) = result%intercept(1, l) + &
            dot_product(xm, result%beta(:, 1, l))
      end do
      ! Convert all retained intercepts to the original x scale.
      do l = 1, result%nlambda
         result%intercept(1, l) = intercept_std - &
            dot_product(xm, result%beta(:, 1, l))
      end do
   end subroutine fit_irls_path

   subroutine fit_custom_family_path(x, y, family_working, result, control, &
      weights_in, offset_in, lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y(:)
      procedure(glmnet_family_working_interface) :: family_working
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      type(glmnet_control_type) :: ctl
      real(dp), allocatable :: weights(:), offset(:), xw(:,:), xm(:), xs(:)
      real(dp), allocatable :: penalty(:), lower(:), upper(:), lambda(:), beta_std(:)
      real(dp), allocatable :: eta(:), working(:), irls_weight(:)
      logical, allocatable :: usable(:), excluded(:)
      real(dp) :: lambda_max, intercept_std, deviance, nulldev, objective, old_objective
      integer :: n, p, j, l, outer, status, inner_iter, inner_passes
      logical :: converged, inner_converged

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      p = size(x, 2)
      call initialize_result(result, 0, n, p)
      result%family = 'custom'
      if (n < 2 .or. p < 1 .or. size(y) /= n .or. .not. all_finite_vector(y)) then
         result%status = glmnet_invalid_argument
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
      allocate(beta_std(p), eta(n), working(n), irls_weight(n))
      beta_std = 0.0_dp
      if (ctl%intercept) then
         intercept_std = weighted_mean(y, weights) - weighted_mean(offset, weights)
      else
         intercept_std = 0.0_dp
      end if
      eta = intercept_std + offset
      call family_working(y, eta, weights, working, irls_weight, nulldev, status)
      if (status /= glmnet_success .or. nulldev <= 0.0_dp) then
         result%status = glmnet_invalid_response
         return
      end if
      lambda_max = 0.0_dp
      do j = 1, p
         if (.not. excluded(j)) lambda_max = max(lambda_max, &
            abs(sum(irls_weight * xw(:, j) * (working - intercept_std))) / &
            max(max(ctl%alpha, 1.0e-3_dp) * penalty(j), glmnet_eps))
      end do
      call make_lambda_sequence(lambda_max, n, p, ctl, lambda_in, lambda, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call allocate_result(result, size(lambda), p)
      result%lambda = lambda
      result%nulldev = nulldev
      result%x_mean = xm
      result%x_scale = xs
      result%standardize = ctl%standardize
      result%intercept_fitted = ctl%intercept
      deviance = nulldev
      objective = 0.5_dp * nulldev
      old_objective = huge(1.0_dp)
      do l = 1, size(lambda)
         converged = .false.
         do outer = 1, ctl%max_outer_iterations
            eta = intercept_std + matmul(xw, beta_std) + offset
            call family_working(y, eta, weights, working, irls_weight, deviance, status)
            if (status /= glmnet_success) exit
            call solve_weighted_elastic_net(xw, working - offset, irls_weight, &
               lambda(l), ctl%alpha, penalty, lower, upper, excluded, ctl%intercept, &
               beta_std, intercept_std, inner_iter, inner_passes, inner_converged, &
               ctl%threshold, ctl%max_iterations)
            result%npasses = result%npasses + inner_passes
            eta = intercept_std + matmul(xw, beta_std) + offset
            call family_working(y, eta, weights, working, irls_weight, deviance, status)
            if (status /= glmnet_success) exit
            objective = 0.5_dp * deviance + penalty_value(beta_std, penalty, lambda(l), ctl%alpha)
            if (outer > 1 .and. abs(result%objective(l) - objective) <= &
                ctl%threshold * (1.0_dp + abs(objective))) then
               converged = .true.
               exit
            end if
            result%objective(l) = objective
         end do
         if (status /= glmnet_success) then
            result%status = glmnet_numerical_failure
            call truncate_result(result, l)
            return
         end if
         result%objective(l) = objective
         result%dev_ratio(l) = max(0.0_dp, min(1.0_dp, 1.0_dp - deviance / nulldev))
         result%iterations(l) = min(outer, ctl%max_outer_iterations)
         result%converged(l) = converged
         result%beta(:, 1, l) = beta_std / xs
         result%intercept(1, l) = intercept_std - dot_product(xm, result%beta(:, 1, l))
         result%df(l) = count(abs(result%beta(:, 1, l)) > 1.0e-12_dp)
         if (.not. converged .and. result%status == glmnet_success) result%status = glmnet_max_iterations
         old_objective = objective
      end do
   end subroutine fit_custom_family_path

   pure function binomial_deviance(y, mu, weights, pmin) result(value)
      real(dp), intent(in) :: y(:), mu(:), weights(:), pmin
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         if (y(i) > pmin) value = value + 2.0_dp * weights(i) * y(i) * &
            safe_log(y(i) / max(mu(i), pmin))
         if (y(i) < 1.0_dp - pmin) value = value + 2.0_dp * weights(i) * &
            (1.0_dp - y(i)) * safe_log((1.0_dp - y(i)) / max(1.0_dp - mu(i), pmin))
      end do
   end function binomial_deviance

   pure function poisson_deviance(y, mu, weights) result(value)
      real(dp), intent(in) :: y(:), mu(:), weights(:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, size(y)
         if (y(i) > 0.0_dp) then
            value = value + 2.0_dp * weights(i) * &
               (y(i) * safe_log(y(i) / max(mu(i), glmnet_eps)) - (y(i) - mu(i)))
         else
            value = value + 2.0_dp * weights(i) * mu(i)
         end if
      end do
   end function poisson_deviance

   pure function penalty_value(beta, penalty, lambda, alpha) result(value)
      real(dp), intent(in) :: beta(:), penalty(:), lambda, alpha
      real(dp) :: value
      value = lambda * sum(penalty * (alpha * abs(beta) + &
         0.5_dp * (1.0_dp - alpha) * beta ** 2))
   end function penalty_value

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

   subroutine initialize_result(result, family_code, n, p)
      type(glmnet_path_result), intent(out) :: result
      integer, intent(in) :: family_code, n, p
      result%family_code = family_code
      result%family = family_name(family_code)
      result%status = glmnet_success
      result%nobs = n
      result%nvars = p
      result%nout = 1
      result%nlambda = 0
      result%npasses = 0
   end subroutine initialize_result

   subroutine allocate_result(result, nlambda, p)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: nlambda, p
      result%nlambda = nlambda
      allocate(result%lambda(nlambda), result%intercept(1, nlambda), &
         result%beta(p, 1, nlambda), result%dev_ratio(nlambda), &
         result%objective(nlambda), result%df(nlambda), result%iterations(nlambda), &
         result%converged(nlambda))
      result%intercept = 0.0_dp
      result%beta = 0.0_dp
      result%dev_ratio = 0.0_dp
      result%objective = 0.0_dp
      result%df = 0
      result%iterations = 0
      result%converged = .false.
   end subroutine allocate_result

   subroutine truncate_result(result, keep)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: keep
      real(dp), allocatable :: lambda(:), intercept(:,:), beta(:,:,:), dev(:), objective(:)
      integer, allocatable :: df(:), iterations(:)
      logical, allocatable :: converged(:)
      if (keep >= result%nlambda) return
      allocate(lambda(keep), intercept(1, keep), beta(result%nvars, 1, keep), &
         dev(keep), objective(keep), df(keep), iterations(keep), converged(keep))
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
   end subroutine truncate_result
end module glmnet_glm
