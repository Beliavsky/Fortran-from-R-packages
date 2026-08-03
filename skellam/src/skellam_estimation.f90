! SPDX-License-Identifier: GPL-2.0-or-later
module skellam_estimation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use skellam_kinds, only : dp, i8, sqrt_two
   use skellam_distribution, only : skellam_log_pmf
   use skellam_optimization, only : optimizer_result, minimize_bfgs, minimize_nelder_mead, &
      numerical_hessian, invert_matrix, solve_linear_system
   implicit none
   private

   type, public :: skellam_mle_result
      real(dp) :: lambda1 = 0.0_dp
      real(dp) :: lambda2 = 0.0_dp
      real(dp) :: standard_error1 = 0.0_dp
      real(dp) :: standard_error2 = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: covariance(2,2) = 0.0_dp
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 0
      logical :: converged = .false.
   end type skellam_mle_result

   type, public :: skellam_regression_result
      real(dp), allocatable :: beta1(:)
      real(dp), allocatable :: beta2(:)
      real(dp), allocatable :: standard_error1(:)
      real(dp), allocatable :: standard_error2(:)
      real(dp), allocatable :: wald1(:)
      real(dp), allocatable :: wald2(:)
      real(dp), allocatable :: p_value1(:)
      real(dp), allocatable :: p_value2(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: fitted_lambda1(:)
      real(dp), allocatable :: fitted_lambda2(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = 0
      logical :: converged = .false.
      logical :: includes_intercept = .true.
      logical :: response_is_lambda1_minus_lambda2 = .false.
   end type skellam_regression_result

   public :: fit_skellam_mle, fit_skellam_regression
   public :: skellam_log_likelihood

contains

   real(dp) function skellam_log_likelihood(observations, lambda1, lambda2, status) result(log_likelihood)
      integer(i8), intent(in) :: observations(:)
      real(dp), intent(in) :: lambda1, lambda2
      integer, intent(out), optional :: status
      integer :: i, local_status
      real(dp) :: term

      if (present(status)) status = 0
      log_likelihood = 0.0_dp
      do i = 1, size(observations)
         term = skellam_log_pmf(observations(i), lambda1, lambda2, local_status)
         if (local_status /= 0 .or. .not. ieee_is_finite(term)) then
            log_likelihood = -huge(1.0_dp)
            if (present(status)) status = 1
            return
         end if
         log_likelihood = log_likelihood + term
      end do
   end function skellam_log_likelihood

   subroutine fit_skellam_mle(observations, result, tolerance, max_iterations)
      integer(i8), intent(in) :: observations(:)
      type(skellam_mle_result), intent(out) :: result
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      type(optimizer_result) :: bfgs_result, simplex_result
      real(dp) :: sample_mean, sample_variance, lower_bound, lambda2_start
      real(dp) :: initial(1), optimum(1), tol
      real(dp) :: full_parameters(2), hessian(2,2), inverse_hessian(2,2)
      real(dp) :: jacobian(2,2)
      integer :: n, max_iter, hessian_status, inverse_status

      result = skellam_mle_result()
      n = size(observations)
      if (n < 2) then
         result%status = 1
         return
      end if
      sample_mean = sum(real(observations, dp))/real(n, dp)
      sample_variance = sum((real(observations, dp) - sample_mean)**2)/real(n - 1, dp)
      lower_bound = max(0.0_dp, -sample_mean) + 1.0e-10_dp
      lambda2_start = max(lower_bound + 0.05_dp, 0.5_dp*(sample_variance - sample_mean))
      initial(1) = log(max(lambda2_start - lower_bound, 1.0e-6_dp))
      tol = 1.0e-8_dp
      max_iter = 500
      if (present(tolerance)) tol = tolerance
      if (present(max_iterations)) max_iter = max_iterations

      call minimize_bfgs(reduced_objective, initial, bfgs_result, tol, max_iter)
      if (.not. bfgs_result%converged) then
         call minimize_nelder_mead(reduced_objective, initial, simplex_result, tol, 4*max_iter)
         if (simplex_result%objective < bfgs_result%objective) then
            call minimize_bfgs(reduced_objective, simplex_result%parameters, bfgs_result, tol, max_iter)
            if (.not. bfgs_result%converged .and. simplex_result%converged) bfgs_result = simplex_result
         end if
      end if

      optimum = bfgs_result%parameters
      result%lambda2 = lower_bound + exp(optimum(1))
      result%lambda1 = result%lambda2 + sample_mean
      result%log_likelihood = -bfgs_result%objective
      result%iterations = bfgs_result%iterations
      result%evaluations = bfgs_result%evaluations
      result%converged = bfgs_result%converged
      result%status = bfgs_result%status

      full_parameters = [log(max(result%lambda1, tiny(1.0_dp))), &
                         log(max(result%lambda2, tiny(1.0_dp)))]
      call numerical_hessian(full_objective, full_parameters, hessian, hessian_status)
      call invert_matrix(hessian, inverse_hessian, inverse_status)
      if (hessian_status == 0 .and. inverse_status == 0) then
         jacobian = 0.0_dp
         jacobian(1,1) = result%lambda1
         jacobian(2,2) = result%lambda2
         result%covariance = matmul(jacobian, matmul(inverse_hessian, transpose(jacobian)))
         result%standard_error1 = sqrt(max(0.0_dp, result%covariance(1,1)))
         result%standard_error2 = sqrt(max(0.0_dp, result%covariance(2,2)))
      else
         result%covariance = quiet_nan()
         result%standard_error1 = quiet_nan()
         result%standard_error2 = quiet_nan()
      end if

   contains

      real(dp) function reduced_objective(parameters) result(value)
         real(dp), intent(in) :: parameters(:)
         real(dp) :: lambda1, lambda2
         lambda2 = lower_bound + exp(min(700.0_dp, parameters(1)))
         lambda1 = lambda2 + sample_mean
         if (lambda1 <= 0.0_dp .or. .not. ieee_is_finite(lambda1 + lambda2)) then
            value = huge(1.0_dp)
         else
            value = -skellam_log_likelihood(observations, lambda1, lambda2)
         end if
      end function reduced_objective

      real(dp) function full_objective(parameters) result(value)
         real(dp), intent(in) :: parameters(:)
         real(dp) :: lambda1, lambda2
         if (maxval(parameters) > 700.0_dp) then
            value = huge(1.0_dp)
            return
         end if
         lambda1 = exp(parameters(1))
         lambda2 = exp(parameters(2))
         value = -skellam_log_likelihood(observations, lambda1, lambda2)
      end function full_objective

   end subroutine fit_skellam_mle

   subroutine fit_skellam_regression(response, predictors, result, include_intercept, &
         response_is_lambda1_minus_lambda2, tolerance, max_iterations)
      integer(i8), intent(in) :: response(:)
      real(dp), intent(in) :: predictors(:,:)
      type(skellam_regression_result), intent(out) :: result
      logical, intent(in), optional :: include_intercept
      logical, intent(in), optional :: response_is_lambda1_minus_lambda2
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations
      real(dp), allocatable :: design(:,:), initial(:), hessian(:,:), inverse_hessian(:,:)
      real(dp), allocatable :: start1(:), start2(:), target1(:), target2(:)
      type(optimizer_result) :: bfgs_result, simplex_result
      real(dp) :: sample_mean, sample_variance, base1, base2, tol
      integer :: n, p_raw, p, i, max_iter, hessian_status, inverse_status, solve_status
      logical :: add_intercept, direct_order

      result = skellam_regression_result()
      n = size(response)
      p_raw = size(predictors,2)
      if (size(predictors,1) /= n .or. n < 2) then
         result%status = 1
         return
      end if
      add_intercept = .true.
      direct_order = .false.
      if (present(include_intercept)) add_intercept = include_intercept
      if (present(response_is_lambda1_minus_lambda2)) direct_order = response_is_lambda1_minus_lambda2
      p = p_raw + merge(1, 0, add_intercept)
      if (p < 1) then
         result%status = 1
         return
      end if
      allocate(design(n,p))
      if (add_intercept) then
         design(:,1) = 1.0_dp
         if (p_raw > 0) design(:,2:p) = predictors
      else
         design = predictors
      end if

      sample_mean = sum(real(response, dp))/real(n, dp)
      sample_variance = sum((real(response, dp) - sample_mean)**2)/real(n - 1, dp)
      if (direct_order) then
         base1 = max(0.05_dp, 0.5_dp*(sample_variance + sample_mean))
         base2 = max(0.05_dp, 0.5_dp*(sample_variance - sample_mean))
      else
         base1 = max(0.05_dp, 0.5_dp*(sample_variance - sample_mean))
         base2 = max(0.05_dp, 0.5_dp*(sample_variance + sample_mean))
      end if

      allocate(initial(2*p), start1(p), start2(p), target1(n), target2(n))
      if (direct_order) then
         target1 = log(max(real(response, dp), 0.0_dp) + 0.5_dp)
         target2 = log(max(-real(response, dp), 0.0_dp) + 0.5_dp)
      else
         target1 = log(max(-real(response, dp), 0.0_dp) + 0.5_dp)
         target2 = log(max(real(response, dp), 0.0_dp) + 0.5_dp)
      end if
      call ridge_least_squares(design, target1, start1, solve_status)
      if (solve_status /= 0) start1 = 0.0_dp
      call ridge_least_squares(design, target2, start2, solve_status)
      if (solve_status /= 0) start2 = 0.0_dp
      if (add_intercept) then
         start1(1) = log(base1)
         start2(1) = log(base2)
      else
         start1 = 0.25_dp*start1
         start2 = 0.25_dp*start2
      end if
      if (p > merge(1, 0, add_intercept)) then
         if (add_intercept) then
            start1(2:p) = 0.25_dp*start1(2:p)
            start2(2:p) = 0.25_dp*start2(2:p)
         end if
      end if
      initial(1:p) = start1
      initial(p + 1:2*p) = start2

      tol = 1.0e-7_dp
      max_iter = 1000
      if (present(tolerance)) tol = tolerance
      if (present(max_iterations)) max_iter = max_iterations
      call minimize_bfgs(regression_objective, initial, bfgs_result, tol, max_iter)
      if (.not. bfgs_result%converged) then
         call minimize_nelder_mead(regression_objective, bfgs_result%parameters, simplex_result, &
            10.0_dp*tol, 3*max_iter)
         if (simplex_result%objective < bfgs_result%objective) then
            call minimize_bfgs(regression_objective, simplex_result%parameters, bfgs_result, tol, max_iter)
            if (.not. bfgs_result%converged .and. simplex_result%converged) bfgs_result = simplex_result
         end if
      end if

      allocate(result%beta1(p), result%beta2(p), result%standard_error1(p), result%standard_error2(p))
      allocate(result%wald1(p), result%wald2(p), result%p_value1(p), result%p_value2(p))
      allocate(result%covariance(2*p,2*p), result%fitted_lambda1(n), result%fitted_lambda2(n))
      result%beta1 = bfgs_result%parameters(1:p)
      result%beta2 = bfgs_result%parameters(p + 1:2*p)
      result%fitted_lambda1 = exp(matmul(design, result%beta1))
      result%fitted_lambda2 = exp(matmul(design, result%beta2))
      result%log_likelihood = -bfgs_result%objective
      result%iterations = bfgs_result%iterations
      result%evaluations = bfgs_result%evaluations
      result%converged = bfgs_result%converged
      result%status = bfgs_result%status
      result%includes_intercept = add_intercept
      result%response_is_lambda1_minus_lambda2 = direct_order

      allocate(hessian(2*p,2*p), inverse_hessian(2*p,2*p))
      call numerical_hessian(regression_objective, bfgs_result%parameters, hessian, hessian_status)
      call invert_matrix(hessian, inverse_hessian, inverse_status)
      if (hessian_status == 0 .and. inverse_status == 0) then
         result%covariance = inverse_hessian
         do i = 1, p
            result%standard_error1(i) = sqrt(max(0.0_dp, inverse_hessian(i,i)))
            result%standard_error2(i) = sqrt(max(0.0_dp, inverse_hessian(p + i,p + i)))
         end do
         result%wald1 = result%beta1/result%standard_error1
         result%wald2 = result%beta2/result%standard_error2
         result%p_value1 = erfc(abs(result%wald1)/sqrt_two)
         result%p_value2 = erfc(abs(result%wald2)/sqrt_two)
      else
         result%covariance = quiet_nan()
         result%standard_error1 = quiet_nan()
         result%standard_error2 = quiet_nan()
         result%wald1 = quiet_nan()
         result%wald2 = quiet_nan()
         result%p_value1 = quiet_nan()
         result%p_value2 = quiet_nan()
      end if

   contains

      real(dp) function regression_objective(parameters) result(value)
         real(dp), intent(in) :: parameters(:)
         real(dp), allocatable :: eta1(:), eta2(:), lambda1(:), lambda2(:)
         real(dp) :: log_probability
         integer :: j

         allocate(eta1(n), eta2(n), lambda1(n), lambda2(n))
         eta1 = matmul(design, parameters(1:p))
         eta2 = matmul(design, parameters(p + 1:2*p))
         if (maxval(abs(eta1)) > 60.0_dp .or. maxval(abs(eta2)) > 60.0_dp) then
            value = huge(1.0_dp)
            return
         end if
         lambda1 = exp(eta1)
         lambda2 = exp(eta2)
         value = 0.0_dp
         do j = 1, n
            if (direct_order) then
               log_probability = skellam_log_pmf(response(j), lambda1(j), lambda2(j))
            else
               log_probability = skellam_log_pmf(response(j), lambda2(j), lambda1(j))
            end if
            if (.not. ieee_is_finite(log_probability)) then
               value = huge(1.0_dp)
               return
            end if
            value = value - log_probability
         end do
      end function regression_objective

   end subroutine fit_skellam_regression

   subroutine ridge_least_squares(design, target, coefficients, status)
      real(dp), intent(in) :: design(:,:), target(:)
      real(dp), intent(out) :: coefficients(:)
      integer, intent(out) :: status
      real(dp), allocatable :: crossproduct(:,:), rhs(:)
      integer :: i, p

      p = size(design,2)
      allocate(crossproduct(p,p), rhs(p))
      crossproduct = matmul(transpose(design), design)
      do i = 1, p
         crossproduct(i,i) = crossproduct(i,i) + 1.0e-6_dp*max(1.0_dp, crossproduct(i,i))
      end do
      rhs = matmul(transpose(design), target)
      call solve_linear_system(crossproduct, rhs, coefficients, status)
   end subroutine ridge_least_squares

   pure real(dp) function quiet_nan() result(value)
      value = ieee_value(value, ieee_quiet_nan)
   end function quiet_nan

end module skellam_estimation
