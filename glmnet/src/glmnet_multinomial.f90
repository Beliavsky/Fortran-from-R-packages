! SPDX-License-Identifier: GPL-2.0-only
module glmnet_multinomial
   use glmnet_kinds, only : dp, glmnet_eps, glmnet_huge
   use glmnet_status, only : glmnet_success, glmnet_invalid_argument, &
      glmnet_nonfinite_input, glmnet_invalid_response, glmnet_max_iterations
   use glmnet_types, only : glmnet_control_type, glmnet_path_result, &
      glmnet_family_multinomial, family_name
   use glmnet_utils, only : normalize_weights, prepare_design, make_lambda_sequence, &
      all_finite_matrix, all_finite_vector, group_soft_threshold, soft_threshold
   implicit none
   private
   public :: fit_multinomial_path, fit_multinomial_matrix_path
contains
   subroutine fit_multinomial_path(x, class_id, result, control, weights_in, &
      offset_in, lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in, nclass)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: class_id(:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:,:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      integer, intent(in), optional :: nclass
      real(dp), allocatable :: y(:,:)
      integer :: n, k, i
      n = size(x, 1)
      if (present(nclass)) then
         k = nclass
      else if (size(class_id) > 0) then
         k = maxval(class_id)
      else
         k = 0
      end if
      if (size(class_id) /= n .or. k < 2) then
         call initialize_invalid(result, n, size(x, 2))
         return
      end if
      if (any(class_id < 1) .or. any(class_id > k)) then
         call initialize_invalid(result, n, size(x, 2))
         return
      end if
      allocate(y(n, k))
      y = 0.0_dp
      do i = 1, n
         y(i, class_id(i)) = 1.0_dp
      end do
      call fit_multinomial_matrix_path(x, y, result, control, weights_in, offset_in, &
         lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      allocate(result%class_levels(k))
      result%class_levels = [(real(i, dp), i=1,k)]
   end subroutine fit_multinomial_path

   subroutine fit_multinomial_matrix_path(x, y_in, result, control, weights_in, &
      offset_in, lambda_in, penalty_factor_in, lower_in, upper_in, excluded_in)
      real(dp), intent(in) :: x(:,:), y_in(:,:)
      type(glmnet_path_result), intent(out) :: result
      type(glmnet_control_type), intent(in), optional :: control
      real(dp), intent(in), optional :: weights_in(:), offset_in(:,:), lambda_in(:)
      real(dp), intent(in), optional :: penalty_factor_in(:), lower_in(:), upper_in(:)
      logical, intent(in), optional :: excluded_in(:)
      type(glmnet_control_type) :: ctl
      real(dp), allocatable :: weights(:), row_total(:), y(:,:), offset(:,:)
      real(dp), allocatable :: xw(:,:), xm(:), xs(:), penalty(:), lower(:), upper(:)
      real(dp), allocatable :: lambda(:), beta(:,:), a0(:), probs(:,:), eta(:,:)
      real(dp), allocatable :: grad_beta(:,:), grad_a0(:), beta_new(:,:), a0_new(:)
      logical, allocatable :: usable(:), excluded(:)
      real(dp) :: lambda_max, nulldev, deviance, objective, objective_new
      real(dp) :: step, max_change, alpha_sequence, norm_gradient
      integer :: n, p, k, i, j, c, l, iter, status, backtrack
      logical :: converged

      ctl = glmnet_control_type()
      if (present(control)) ctl = control
      n = size(x, 1)
      p = size(x, 2)
      k = size(y_in, 2)
      call initialize_result(result, n, p, k)
      if (n < 2 .or. p < 1 .or. k < 2 .or. size(y_in, 1) /= n) then
         result%status = glmnet_invalid_argument
         return
      end if
      if (.not. all_finite_matrix(y_in) .or. any(y_in < 0.0_dp)) then
         result%status = glmnet_invalid_response
         return
      end if
      allocate(row_total(n), y(n, k))
      row_total = sum(y_in, dim=2)
      if (any(row_total <= 0.0_dp)) then
         result%status = glmnet_invalid_response
         return
      end if
      do i = 1, n
         y(i, :) = y_in(i, :) / row_total(i)
      end do
      call normalize_weights(n, weights_in, weights, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      weights = weights * row_total
      weights = weights / sum(weights)
      allocate(offset(n, k))
      if (present(offset_in)) then
         if (size(offset_in, 1) /= n .or. size(offset_in, 2) /= k) then
            result%status = glmnet_invalid_argument
            return
         end if
         if (.not. all_finite_matrix(offset_in)) then
            result%status = glmnet_nonfinite_input
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
      allocate(beta(p, k), a0(k), probs(n, k), eta(n, k), grad_beta(p, k), &
         grad_a0(k), beta_new(p, k), a0_new(k))
      beta = 0.0_dp
      do c = 1, k
         a0(c) = log(max(sum(weights * y(:, c)), ctl%probability_min))
      end do
      a0 = a0 - sum(a0) / real(k, dp)
      call softmax_matrix(spread(a0, 1, n) + offset, probs)
      nulldev = multinomial_deviance(y, probs, weights, ctl%probability_min)
      if (nulldev <= 100.0_dp * glmnet_eps) then
         result%status = glmnet_invalid_response
         return
      end if
      call multinomial_gradient(xw, y, probs, weights, grad_beta, grad_a0)
      alpha_sequence = max(ctl%alpha, 1.0e-3_dp)
      lambda_max = 0.0_dp
      do j = 1, p
         if (excluded(j)) cycle
         if (ctl%grouped) then
            norm_gradient = sqrt(sum(grad_beta(j, :) ** 2))
         else
            norm_gradient = maxval(abs(grad_beta(j, :)))
         end if
         lambda_max = max(lambda_max, norm_gradient / &
            max(alpha_sequence * penalty(j), glmnet_eps))
      end do
      call make_lambda_sequence(lambda_max, n, p, ctl, lambda_in, lambda, status)
      if (status /= glmnet_success) then
         result%status = status
         return
      end if
      call allocate_result(result, size(lambda), p, k)
      result%lambda = lambda
      result%nulldev = nulldev
      result%x_mean = xm
      result%x_scale = xs
      result%standardize = ctl%standardize
      result%intercept_fitted = ctl%intercept
      deviance = nulldev
      objective_new = 0.5_dp * nulldev
      do l = 1, size(lambda)
         converged = .false.
         eta = matmul(xw, beta) + spread(a0, 1, n) + offset
         call softmax_matrix(eta, probs)
         objective = 0.5_dp * multinomial_deviance(y, probs, weights, &
            ctl%probability_min) + multinomial_penalty(beta, penalty, lambda(l), &
            ctl%alpha, ctl%grouped)
         step = initial_step_size(xw, weights, lambda(l), ctl%alpha, penalty)
         do iter = 1, ctl%max_iterations
            call multinomial_gradient(xw, y, probs, weights, grad_beta, grad_a0)
            a0_new = a0
            if (ctl%intercept) a0_new = a0 - step * grad_a0
            a0_new = a0_new - sum(a0_new) / real(k, dp)
            beta_new = beta - step * grad_beta
            call apply_multinomial_prox(beta_new, step, lambda(l), ctl%alpha, &
               penalty, lower, upper, excluded, ctl%grouped)
            do j = 1, p
               beta_new(j, :) = beta_new(j, :) - sum(beta_new(j, :)) / real(k, dp)
            end do
            do backtrack = 1, 40
               eta = matmul(xw, beta_new) + spread(a0_new, 1, n) + offset
               call softmax_matrix(eta, probs)
               deviance = multinomial_deviance(y, probs, weights, ctl%probability_min)
               objective_new = 0.5_dp * deviance + multinomial_penalty(beta_new, &
                  penalty, lambda(l), ctl%alpha, ctl%grouped)
               if (objective_new <= objective + 1.0e-12_dp) exit
               step = 0.5_dp * step
               a0_new = a0
               if (ctl%intercept) a0_new = a0 - step * grad_a0
               a0_new = a0_new - sum(a0_new) / real(k, dp)
               beta_new = beta - step * grad_beta
               call apply_multinomial_prox(beta_new, step, lambda(l), ctl%alpha, &
                  penalty, lower, upper, excluded, ctl%grouped)
               do j = 1, p
                  beta_new(j, :) = beta_new(j, :) - sum(beta_new(j, :)) / real(k, dp)
               end do
            end do
            max_change = max(maxval(abs(beta_new - beta)), maxval(abs(a0_new - a0)))
            beta = beta_new
            a0 = a0_new
            objective = objective_new
            result%npasses = result%npasses + 1
            if (max_change <= ctl%threshold * (1.0_dp + maxval(abs(beta)))) then
               converged = .true.
               exit
            end if
            step = min(step * 1.2_dp, 1.0_dp)
         end do
         result%objective(l) = objective
         result%dev_ratio(l) = max(0.0_dp, min(1.0_dp, 1.0_dp - deviance / nulldev))
         result%iterations(l) = min(iter, ctl%max_iterations)
         result%converged(l) = converged
         do c = 1, k
            result%beta(:, c, l) = beta(:, c) / xs
            result%intercept(c, l) = a0(c) - dot_product(xm, result%beta(:, c, l))
         end do
         result%df(l) = count_rows(result%beta(:, :, l))
         if (.not. converged .and. result%status == glmnet_success) &
            result%status = glmnet_max_iterations
      end do
   end subroutine fit_multinomial_matrix_path

   subroutine softmax_matrix(eta, probabilities)
      real(dp), intent(in) :: eta(:,:)
      real(dp), intent(out) :: probabilities(size(eta, 1), size(eta, 2))
      real(dp) :: maximum, total
      integer :: i
      do i = 1, size(eta, 1)
         maximum = maxval(eta(i, :))
         probabilities(i, :) = exp(max(eta(i, :) - maximum, -700.0_dp))
         total = sum(probabilities(i, :))
         probabilities(i, :) = probabilities(i, :) / max(total, glmnet_eps)
      end do
   end subroutine softmax_matrix

   subroutine multinomial_gradient(x, y, probabilities, weights, grad_beta, grad_a0)
      real(dp), intent(in) :: x(:,:), y(:,:), probabilities(:,:), weights(:)
      real(dp), intent(out) :: grad_beta(size(x, 2), size(y, 2)), grad_a0(size(y, 2))
      real(dp), allocatable :: error(:,:)
      allocate(error(size(y, 1), size(y, 2)))
      error = spread(weights, 2, size(y, 2)) * (probabilities - y)
      grad_beta = matmul(transpose(x), error)
      grad_a0 = sum(error, dim=1)
   end subroutine multinomial_gradient

   subroutine apply_multinomial_prox(beta, step, lambda, alpha, penalty, lower, &
      upper, excluded, grouped)
      real(dp), intent(inout) :: beta(:,:)
      real(dp), intent(in) :: step, lambda, alpha, penalty(:), lower(:), upper(:)
      logical, intent(in) :: excluded(:), grouped
      real(dp), allocatable :: row(:)
      real(dp) :: denominator
      integer :: j, c
      allocate(row(size(beta, 2)))
      do j = 1, size(beta, 1)
         if (excluded(j)) then
            beta(j, :) = 0.0_dp
            cycle
         end if
         denominator = 1.0_dp + step * lambda * (1.0_dp - alpha) * penalty(j)
         if (grouped) then
            call group_soft_threshold(beta(j, :), step * lambda * alpha * penalty(j), row)
            beta(j, :) = row / denominator
         else
            do c = 1, size(beta, 2)
               beta(j, c) = soft_threshold(beta(j, c), step * lambda * alpha * penalty(j)) / denominator
            end do
         end if
         beta(j, :) = min(max(beta(j, :), lower(j)), upper(j))
      end do
   end subroutine apply_multinomial_prox

   pure function multinomial_deviance(y, probabilities, weights, pmin) result(value)
      real(dp), intent(in) :: y(:,:), probabilities(:,:), weights(:), pmin
      real(dp) :: value
      integer :: i, c
      value = 0.0_dp
      do i = 1, size(y, 1)
         do c = 1, size(y, 2)
            if (y(i, c) > 0.0_dp) value = value - 2.0_dp * weights(i) * &
               y(i, c) * log(max(probabilities(i, c), pmin))
         end do
      end do
   end function multinomial_deviance

   pure function multinomial_penalty(beta, penalty, lambda, alpha, grouped) result(value)
      real(dp), intent(in) :: beta(:,:), penalty(:), lambda, alpha
      logical, intent(in) :: grouped
      real(dp) :: value
      integer :: j
      value = 0.0_dp
      do j = 1, size(beta, 1)
         if (grouped) then
            value = value + lambda * penalty(j) * (alpha * sqrt(sum(beta(j, :) ** 2)) + &
               0.5_dp * (1.0_dp - alpha) * sum(beta(j, :) ** 2))
         else
            value = value + lambda * penalty(j) * (alpha * sum(abs(beta(j, :))) + &
               0.5_dp * (1.0_dp - alpha) * sum(beta(j, :) ** 2))
         end if
      end do
   end function multinomial_penalty

   pure function initial_step_size(x, weights, lambda, alpha, penalty) result(step)
      real(dp), intent(in) :: x(:,:), weights(:), lambda, alpha, penalty(:)
      real(dp) :: step, lipschitz
      integer :: j
      lipschitz = 0.0_dp
      do j = 1, size(x, 2)
         lipschitz = max(lipschitz, sum(weights * x(:, j) ** 2))
      end do
      lipschitz = 0.5_dp * max(lipschitz, glmnet_eps) + &
         lambda * (1.0_dp - alpha) * maxval(penalty)
      step = 1.0_dp / max(lipschitz, glmnet_eps)
   end function initial_step_size

   pure function count_rows(beta) result(value)
      real(dp), intent(in) :: beta(:,:)
      integer :: value, j
      value = 0
      do j = 1, size(beta, 1)
         if (maxval(abs(beta(j, :))) > 1.0e-12_dp) value = value + 1
      end do
   end function count_rows

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

   subroutine initialize_result(result, n, p, k)
      type(glmnet_path_result), intent(out) :: result
      integer, intent(in) :: n, p, k
      result%family_code = glmnet_family_multinomial
      result%family = family_name(glmnet_family_multinomial)
      result%status = glmnet_success
      result%nobs = n
      result%nvars = p
      result%nout = k
      result%nlambda = 0
   end subroutine initialize_result

   subroutine initialize_invalid(result, n, p)
      type(glmnet_path_result), intent(out) :: result
      integer, intent(in) :: n, p
      call initialize_result(result, n, p, 0)
      result%status = glmnet_invalid_argument
   end subroutine initialize_invalid

   subroutine allocate_result(result, nlambda, p, k)
      type(glmnet_path_result), intent(inout) :: result
      integer, intent(in) :: nlambda, p, k
      result%nlambda = nlambda
      allocate(result%lambda(nlambda), result%intercept(k, nlambda), &
         result%beta(p, k, nlambda), result%dev_ratio(nlambda), &
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
end module glmnet_multinomial
