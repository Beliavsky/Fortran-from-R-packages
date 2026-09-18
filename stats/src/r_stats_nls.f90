! SPDX-License-Identifier: MIT
! SPDX-FileComment: Pure formula-free nonlinear least squares corresponding to R stats nls.
module r_stats_nls
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_distributions, only: r_qt
   use r_kinds, only: dp
   use r_linalg, only: inverse_matrix, solve_system
   use r_optional, only: optval
   use r_status, only: r_invalid_input, r_ok, r_singular
   use r_stats_types, only: nls_fit_t
   implicit none
   private

   integer, parameter, public :: nls_max_iterations = 2

   public :: nls_basis_interface, nls_confint_profile, nls_fit, nls_fit_plinear
   public :: nls_jacobian_interface, nls_model_interface

   abstract interface
      pure subroutine nls_model_interface(x, parameters, values)
         !! Evaluates a nonlinear regression model at every observation.
         import dp
         real(dp), intent(in) :: x(:, :) !! Predictor matrix with shape `(n, q)`.
         real(dp), intent(in) :: parameters(:) !! Model parameter vector with size `p`.
         real(dp), intent(out) :: values(size(x, 1)) !! Model values with size `n`.
      end subroutine nls_model_interface

      pure subroutine nls_jacobian_interface(x, parameters, jacobian)
         !! Evaluates derivatives of model values with respect to parameters.
         import dp
         real(dp), intent(in) :: x(:, :) !! Predictor matrix with shape `(n, q)`.
         real(dp), intent(in) :: parameters(:) !! Model parameter vector with size `p`.
         real(dp), intent(out) :: jacobian(size(x, 1), size(parameters)) !! Model Jacobian `(n,p)`.
      end subroutine nls_jacobian_interface

      pure subroutine nls_basis_interface(x, nonlinear_parameters, basis_values)
         !! Evaluates the linear basis of a separable nonlinear model.
         import dp
         real(dp), intent(in) :: x(:, :) !! Predictor matrix with shape `(n, q)`.
         real(dp), intent(in) :: nonlinear_parameters(:) !! Nonlinear parameter vector.
         real(dp), intent(out) :: basis_values(:, :) !! Linear basis matrix with shape `(n,k)`.
      end subroutine nls_basis_interface
   end interface

contains

   pure function nls_fit_plinear(x, y, start, n_linear, basis, tolerance, &
                                 max_iterations, lower, upper) result(fit)
      !! Fits a separable nonlinear model by profiling its linear coefficients.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one observation per row.
      real(dp), intent(in) :: y(:) !! Finite response vector conformable with `x`.
      real(dp), intent(in) :: start(:) !! Starting values for nonlinear parameters.
      integer, intent(in) :: n_linear !! Positive number of linear coefficients.
      procedure(nls_basis_interface) :: basis !! Pure callback returning the linear basis.
      real(dp), intent(in), optional :: tolerance !! Positive relative convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Positive iteration limit.
      real(dp), intent(in), optional :: lower(:) !! Lower bounds for nonlinear parameters.
      real(dp), intent(in), optional :: upper(:) !! Upper bounds for nonlinear parameters.
      type(nls_fit_t) :: fit
      type(nls_fit_t) :: nonlinear_fit
      real(dp), allocatable :: basis_values(:, :), full_jacobian(:, :), linear_parameters(:)
      real(dp), allocatable :: lower_all(:), normal_matrix(:, :), upper_all(:)
      real(dp) :: lower_basis(size(y), n_linear), shifted(size(start))
      real(dp) :: step, upper_basis(size(y), n_linear)
      integer :: info, j, n, n_nonlinear, n_total

      n = size(y)
      n_nonlinear = size(start)
      n_total = n_nonlinear + n_linear
      if (n_linear <= 0 .or. n_nonlinear <= 0 .or. n <= n_total .or. size(x, 1) /= n) then
         fit%status = r_invalid_input
         return
      end if
      nonlinear_fit = nls_fit(x, y, start, profiled_model, tolerance=tolerance, &
                              max_iterations=max_iterations, lower=lower, upper=upper)
      if (.not. allocated(nonlinear_fit%coefficients)) then
         fit = nonlinear_fit
         return
      end if
      call evaluate_profiled_model(x, y, nonlinear_fit%coefficients, basis, n_linear, &
                                   basis_values, linear_parameters, fit%fitted_values, info)
      if (info /= 0) then
         fit%status = r_singular
         return
      end if
      fit = nonlinear_fit
      fit%coefficients = [nonlinear_fit%coefficients, linear_parameters]
      fit%fitted_values = matmul(basis_values, linear_parameters)
      fit%residuals = y - fit%fitted_values
      fit%rss = sum(fit%residuals**2)
      fit%degrees_freedom = n - n_total
      fit%residual_standard_error = sqrt(fit%rss/real(fit%degrees_freedom, dp))
      allocate (full_jacobian(n, n_total))
      do j = 1, n_nonlinear
         step = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(fit%coefficients(j)))
         shifted = fit%coefficients(:n_nonlinear)
         shifted(j) = shifted(j) + step
         call basis(x, shifted, upper_basis)
         shifted(j) = fit%coefficients(j) - step
         call basis(x, shifted, lower_basis)
         full_jacobian(:, j) = matmul(upper_basis - lower_basis, linear_parameters)/(2.0_dp*step)
      end do
      full_jacobian(:, n_nonlinear + 1:) = basis_values
      fit%rank = numerical_rank(full_jacobian)
      if (fit%rank < n_total .or. .not. all(ieee_is_finite(full_jacobian))) then
         fit%status = r_singular
         fit%converged = .false.
         return
      end if
      normal_matrix = matmul(transpose(full_jacobian), full_jacobian)
      call inverse_matrix(normal_matrix, fit%covariance, info)
      if (info /= 0) then
         fit%status = r_singular
         fit%converged = .false.
         return
      end if
      fit%covariance = fit%covariance*fit%residual_standard_error**2
      if (allocated(fit%standard_errors)) deallocate (fit%standard_errors)
      allocate (fit%standard_errors(n_total))
      do j = 1, n_total
         fit%standard_errors(j) = sqrt(max(0.0_dp, fit%covariance(j, j)))
      end do
      allocate (lower_all(n_total), upper_all(n_total))
      lower_all = -huge(1.0_dp)
      upper_all = huge(1.0_dp)
      lower_all(:n_nonlinear) = nonlinear_fit%lower_bounds
      upper_all(:n_nonlinear) = nonlinear_fit%upper_bounds
      fit%lower_bounds = lower_all
      fit%upper_bounds = upper_all

   contains

      pure subroutine profiled_model(model_x, nonlinear_parameters, values)
         !! Profiles the linear coefficients and evaluates the resulting fitted values.
         real(dp), intent(in) :: model_x(:, :) !! Predictor matrix with shape `(n,q)`.
         real(dp), intent(in) :: nonlinear_parameters(:) !! Nonlinear parameter vector.
         real(dp), intent(out) :: values(size(model_x, 1)) !! Profiled fitted values.
         real(dp), allocatable :: local_basis(:, :), local_linear(:)
         integer :: local_info

         call evaluate_profiled_model(model_x, y, nonlinear_parameters, basis, n_linear, &
                                      local_basis, local_linear, values, local_info)
         if (local_info /= 0) values = huge(1.0_dp)
      end subroutine profiled_model

   end function nls_fit_plinear

   pure subroutine evaluate_profiled_model(x, y, nonlinear_parameters, basis, n_linear, &
                                           basis_values, linear_parameters, values, info)
      !! Solves the linear least-squares problem at fixed nonlinear parameters.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with shape `(n,q)`.
      real(dp), intent(in) :: y(:) !! Response vector with size `n`.
      real(dp), intent(in) :: nonlinear_parameters(:) !! Fixed nonlinear parameters.
      procedure(nls_basis_interface) :: basis !! Pure callback returning the linear basis.
      integer, intent(in) :: n_linear !! Number of columns in the linear basis.
      real(dp), allocatable, intent(out) :: basis_values(:, :) !! Evaluated basis `(n,k)`.
      real(dp), allocatable, intent(out) :: linear_parameters(:) !! Profiled coefficients.
      real(dp), intent(out) :: values(size(x, 1)) !! Profiled fitted values.
      integer, intent(out) :: info !! Zero on success or a linear-system error code.
      real(dp), allocatable :: normal_matrix(:, :), right_hand_side(:)

      allocate (basis_values(size(x, 1), n_linear), linear_parameters(n_linear))
      call basis(x, nonlinear_parameters, basis_values)
      if (.not. all(ieee_is_finite(basis_values))) then
         info = r_invalid_input
         return
      end if
      normal_matrix = matmul(transpose(basis_values), basis_values)
      right_hand_side = matmul(transpose(basis_values), y)
      call solve_system(normal_matrix, right_hand_side, linear_parameters, info)
      if (info == 0) values = matmul(basis_values, linear_parameters)
   end subroutine evaluate_profiled_model

   pure function nls_confint_profile(x, y, fit, model, jacobian, level, tolerance, &
                                     max_iterations) result(intervals)
      !! Computes profile-RSS confidence intervals corresponding to R `confint.nls`.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one observation per row.
      real(dp), intent(in) :: y(:) !! Finite response vector conformable with `x`.
      type(nls_fit_t), intent(in) :: fit !! Converged full-model NLS fit to profile.
      procedure(nls_model_interface) :: model !! Pure model-evaluation callback.
      procedure(nls_jacobian_interface), optional :: jacobian !! Optional analytic Jacobian.
      real(dp), intent(in), optional :: level !! Confidence level; defaults to 0.95.
      real(dp), intent(in), optional :: tolerance !! Relative endpoint tolerance.
      integer, intent(in), optional :: max_iterations !! Iteration limit for each refit.
      real(dp), allocatable :: intervals(:, :)
      real(dp) :: confidence_level, endpoint_tolerance, target_rss, t_quantile
      integer :: j, iteration_limit, p

      p = size(fit%coefficients)
      allocate (intervals(p, 2), source=0.0_dp)
      confidence_level = optval(level, 0.95_dp)
      endpoint_tolerance = optval(tolerance, 1.0e-7_dp)
      iteration_limit = optval(max_iterations, 200)
      if (.not. fit%converged .or. fit%degrees_freedom <= 0 .or. p == 0 .or. &
          size(x, 1) /= size(y) .or. confidence_level <= 0.0_dp .or. &
          confidence_level >= 1.0_dp .or. endpoint_tolerance <= 0.0_dp .or. &
          iteration_limit <= 0 .or. .not. allocated(fit%standard_errors) .or. &
          .not. allocated(fit%lower_bounds) .or. .not. allocated(fit%upper_bounds)) then
         intervals = huge(1.0_dp)
         return
      end if
      t_quantile = r_qt(0.5_dp*(1.0_dp + confidence_level), &
                        real(fit%degrees_freedom, dp))
      target_rss = fit%rss + (fit%residual_standard_error*t_quantile)**2
      do j = 1, p
         call profile_endpoint(x, y, fit, model, jacobian, j, -1, target_rss, &
                               endpoint_tolerance, iteration_limit, intervals(j, 1))
         call profile_endpoint(x, y, fit, model, jacobian, j, 1, target_rss, &
                               endpoint_tolerance, iteration_limit, intervals(j, 2))
      end do
   end function nls_confint_profile

   pure subroutine profile_endpoint(x, y, fit, model, jacobian, parameter_index, direction, &
                                    target_rss, tolerance, max_iterations, endpoint)
      !! Locates one profile-RSS confidence bound while refitting all other parameters.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one observation per row.
      real(dp), intent(in) :: y(:) !! Response vector conformable with `x`.
      type(nls_fit_t), intent(in) :: fit !! Converged full-model NLS fit.
      procedure(nls_model_interface) :: model !! Pure model-evaluation callback.
      procedure(nls_jacobian_interface), optional :: jacobian !! Optional analytic Jacobian.
      integer, intent(in) :: parameter_index !! Parameter whose endpoint is sought.
      integer, intent(in) :: direction !! Search direction, either minus one or plus one.
      real(dp), intent(in) :: target_rss !! Profile RSS defining the confidence boundary.
      real(dp), intent(in) :: tolerance !! Relative endpoint tolerance.
      integer, intent(in) :: max_iterations !! Iteration limit for each profile refit.
      real(dp), intent(out) :: endpoint !! Located parameter endpoint.
      real(dp) :: candidate, candidate_rss, center, inside, midpoint, outside, scale
      real(dp) :: inside_rss, outside_rss
      integer :: iteration
      logical :: valid

      center = fit%coefficients(parameter_index)
      scale = max(fit%standard_errors(parameter_index), &
                  sqrt(epsilon(1.0_dp))*(1.0_dp + abs(center)))
      inside = center
      inside_rss = fit%rss
      outside = center
      outside_rss = fit%rss
      do iteration = 1, 60
         candidate = center + real(direction, dp)*scale
         candidate = max(fit%lower_bounds(parameter_index), &
                         min(fit%upper_bounds(parameter_index), candidate))
         call profile_rss(x, y, fit, model, jacobian, parameter_index, candidate, &
                          max_iterations, candidate_rss, valid)
         if (.not. valid) then
            scale = 0.5_dp*scale
            cycle
         end if
         if (candidate_rss >= target_rss) then
            outside = candidate
            outside_rss = candidate_rss
            exit
         end if
         inside = candidate
         inside_rss = candidate_rss
         if (candidate == fit%lower_bounds(parameter_index) .or. &
             candidate == fit%upper_bounds(parameter_index)) then
            endpoint = candidate
            return
         end if
         scale = 1.8_dp*scale
      end do
      if (outside_rss < target_rss) then
         endpoint = outside
         return
      end if
      do iteration = 1, 80
         midpoint = 0.5_dp*(inside + outside)
         call profile_rss(x, y, fit, model, jacobian, parameter_index, midpoint, &
                          max_iterations, candidate_rss, valid)
         if (.not. valid .or. candidate_rss >= target_rss) then
            outside = midpoint
            if (valid) outside_rss = candidate_rss
         else
            inside = midpoint
            inside_rss = candidate_rss
         end if
         if (abs(outside - inside) <= tolerance*(1.0_dp + abs(midpoint))) exit
      end do
      if (outside_rss > inside_rss) then
         endpoint = inside + (outside - inside)*(target_rss - inside_rss)/ &
                    (outside_rss - inside_rss)
      else
         endpoint = 0.5_dp*(inside + outside)
      end if
   end subroutine profile_endpoint

   pure subroutine profile_rss(x, y, fit, model, jacobian, parameter_index, fixed_value, &
                               max_iterations, rss, valid)
      !! Refits an NLS model with one parameter fixed and returns its residual sum of squares.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one observation per row.
      real(dp), intent(in) :: y(:) !! Response vector conformable with `x`.
      type(nls_fit_t), intent(in) :: fit !! Full-model fit supplying starts and bounds.
      procedure(nls_model_interface) :: model !! Pure model-evaluation callback.
      procedure(nls_jacobian_interface), optional :: jacobian !! Optional analytic Jacobian.
      integer, intent(in) :: parameter_index !! Parameter held fixed during the refit.
      real(dp), intent(in) :: fixed_value !! Fixed trial value for that parameter.
      integer, intent(in) :: max_iterations !! Refit iteration limit.
      real(dp), intent(out) :: rss !! Profile residual sum of squares.
      logical, intent(out) :: valid !! Whether the refit produced a finite RSS.
      real(dp) :: lower(size(fit%coefficients)), start(size(fit%coefficients))
      real(dp) :: upper(size(fit%coefficients))
      type(nls_fit_t) :: trial

      start = fit%coefficients
      lower = fit%lower_bounds
      upper = fit%upper_bounds
      start(parameter_index) = fixed_value
      lower(parameter_index) = fixed_value
      upper(parameter_index) = fixed_value
      trial = nls_fit(x, y, start, model, jacobian, tolerance=1.0e-10_dp, &
                      max_iterations=max_iterations, lower=lower, upper=upper)
      rss = trial%rss
      valid = allocated(trial%coefficients) .and. ieee_is_finite(rss)
   end subroutine profile_rss

   pure function nls_fit(x, y, start, model, jacobian, tolerance, max_iterations, &
                         minimum_factor, lower, upper) result(fit)
      !! Fits a bounded or unbounded nonlinear model using projected Gauss-Newton steps.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one observation per row.
      real(dp), intent(in) :: y(:) !! Finite response vector conformable with `x`.
      real(dp), intent(in) :: start(:) !! Finite starting parameter values.
      procedure(nls_model_interface) :: model !! Pure model-evaluation callback.
      procedure(nls_jacobian_interface), optional :: jacobian !! Optional analytic Jacobian.
      real(dp), intent(in), optional :: tolerance !! Positive relative convergence tolerance.
      integer, intent(in), optional :: max_iterations !! Positive iteration limit.
      real(dp), intent(in), optional :: minimum_factor !! Smallest accepted step factor.
      real(dp), intent(in), optional :: lower(:) !! Lower parameter bounds with size `p`.
      real(dp), intent(in), optional :: upper(:) !! Upper parameter bounds with size `p`.
      type(nls_fit_t) :: fit
      real(dp), allocatable :: candidate_values(:), delta(:), gradient(:), jacobian_values(:, :)
      real(dp), allocatable :: candidate_parameters(:), normal_matrix(:, :), parameters(:)
      real(dp), allocatable :: reduced_gradient(:), reduced_matrix(:, :), reduced_step(:), values(:)
      real(dp) :: bound_tolerance, candidate_rss, factor, minimum_step, old_rss
      real(dp) :: step_tolerance, tol
      integer :: free_indices(size(start)), info, iteration, iteration_limit, j, n, n_free, p
      logical :: active(size(start))

      n = size(y)
      p = size(start)
      tol = optval(tolerance, 1.0e-8_dp)
      iteration_limit = optval(max_iterations, 50)
      minimum_step = optval(minimum_factor, 1.0_dp/1024.0_dp)
      if (size(x, 1) /= n .or. n <= p .or. p == 0 .or. tol <= 0.0_dp .or. &
          iteration_limit <= 0 .or. minimum_step <= 0.0_dp .or. minimum_step > 1.0_dp .or. &
          .not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y)) .or. &
          .not. all(ieee_is_finite(start))) then
         fit%status = r_invalid_input
         return
      end if
      if (present(lower)) then
         if (size(lower) /= p .or. any(ieee_is_nan(lower))) then
            fit%status = r_invalid_input
            return
         end if
         fit%lower_bounds = lower
         fit%bounded = .true.
      else
         allocate (fit%lower_bounds(p), source=-huge(1.0_dp))
      end if
      if (present(upper)) then
         if (size(upper) /= p .or. any(ieee_is_nan(upper))) then
            fit%status = r_invalid_input
            return
         end if
         fit%upper_bounds = upper
         fit%bounded = .true.
      else
         allocate (fit%upper_bounds(p), source=huge(1.0_dp))
      end if
      if (any(fit%lower_bounds > fit%upper_bounds) .or. &
          any(start < fit%lower_bounds) .or. any(start > fit%upper_bounds)) then
         fit%status = r_invalid_input
         return
      end if
      allocate (parameters(p), candidate_parameters(p), values(n), candidate_values(n))
      allocate (delta(p), gradient(p))
      allocate (jacobian_values(n, p), normal_matrix(p, p))
      parameters = start
      call model(x, parameters, values)
      if (.not. all(ieee_is_finite(values))) then
         fit%status = r_invalid_input
         return
      end if
      old_rss = sum((y - values)**2)
      fit%analytic_jacobian = present(jacobian)
      do iteration = 1, iteration_limit
         call evaluate_jacobian(x, parameters, model, jacobian, jacobian_values, info)
         if (info /= r_ok) then
            fit%status = info
            return
         end if
         fit%rank = numerical_rank(jacobian_values)
         if (fit%rank < p) then
            fit%status = r_singular
            fit%coefficients = parameters
            fit%fitted_values = values
            fit%residuals = y - values
            fit%rss = sum(fit%residuals**2)
            return
         end if
         normal_matrix = matmul(transpose(jacobian_values), jacobian_values)
         gradient = matmul(transpose(jacobian_values), y - values)
         active = .false.
         do j = 1, p
            bound_tolerance = sqrt(epsilon(1.0_dp))*(1.0_dp + abs(parameters(j)))
            if (parameters(j) <= fit%lower_bounds(j) + bound_tolerance .and. &
                gradient(j) < 0.0_dp) active(j) = .true.
            if (parameters(j) >= fit%upper_bounds(j) - bound_tolerance .and. &
                gradient(j) > 0.0_dp) active(j) = .true.
         end do
         n_free = 0
         do j = 1, p
            if (.not. active(j)) then
               n_free = n_free + 1
               free_indices(n_free) = j
            end if
         end do
         if (n_free == 0) then
            fit%converged = .true.
            fit%iterations = iteration - 1
            exit
         end if
         reduced_matrix = normal_matrix(free_indices(:n_free), free_indices(:n_free))
         reduced_gradient = gradient(free_indices(:n_free))
         if (allocated(reduced_step)) deallocate (reduced_step)
         allocate (reduced_step(n_free))
         call solve_system(reduced_matrix, reduced_gradient, reduced_step, info)
         if (info /= 0) then
            fit%status = r_singular
            return
         end if
         delta = 0.0_dp
         do j = 1, n_free
            delta(free_indices(j)) = reduced_step(j)
         end do
         factor = 1.0_dp
         do
            candidate_parameters = max(fit%lower_bounds, &
                                       min(fit%upper_bounds, parameters + factor*delta))
            call model(x, candidate_parameters, candidate_values)
            if (all(ieee_is_finite(candidate_values))) then
               candidate_rss = sum((y - candidate_values)**2)
            else
               candidate_rss = huge(1.0_dp)
            end if
            if (candidate_rss <= old_rss .or. factor < minimum_step) exit
            factor = 0.5_dp*factor
         end do
         if (factor < minimum_step .and. candidate_rss >= old_rss) then
            fit%iterations = iteration - 1
            exit
         end if
         delta = candidate_parameters - parameters
         parameters = candidate_parameters
         values = candidate_values
         fit%iterations = iteration
         step_tolerance = tol*(sqrt(sum(parameters**2)) + tol)
         if (sqrt(sum(delta**2)) <= step_tolerance .or. &
             abs(old_rss - candidate_rss) <= tol*(1.0_dp + candidate_rss)) then
            fit%converged = .true.
            exit
         end if
         old_rss = candidate_rss
      end do
      if (.not. fit%converged .and. fit%iterations == iteration_limit) then
         fit%status = nls_max_iterations
      else if (.not. fit%converged) then
         fit%status = nls_max_iterations
      else
         fit%status = r_ok
      end if
      fit%coefficients = parameters
      fit%fitted_values = values
      fit%residuals = y - values
      fit%rss = sum(fit%residuals**2)
      fit%degrees_freedom = n - p
      fit%residual_standard_error = sqrt(fit%rss/real(fit%degrees_freedom, dp))
      call evaluate_jacobian(x, parameters, model, jacobian, jacobian_values, info)
      if (info /= r_ok) then
         fit%status = info
         return
      end if
      fit%rank = numerical_rank(jacobian_values)
      if (fit%rank < p) then
         fit%status = r_singular
         return
      end if
      normal_matrix = matmul(transpose(jacobian_values), jacobian_values)
      call inverse_matrix(normal_matrix, fit%covariance, info)
      if (info /= 0) then
         fit%status = r_singular
         return
      end if
      fit%covariance = fit%covariance*fit%residual_standard_error**2
      allocate (fit%standard_errors(p))
      do info = 1, p
         fit%standard_errors(info) = sqrt(max(0.0_dp, fit%covariance(info, info)))
      end do
   end function nls_fit

   pure integer function numerical_rank(matrix) result(rank)
      !! Estimates matrix rank by tolerance-scaled Gaussian elimination with pivoting.
      real(dp), intent(in) :: matrix(:, :) !! Matrix whose numerical rank is required.
      real(dp), allocatable :: work(:, :)
      real(dp) :: pivot_size, scale, temporary, threshold
      integer :: column, i, pivot_row, row

      work = matrix
      scale = maxval(abs(work))
      if (scale == 0.0_dp) then
         rank = 0
         return
      end if
      threshold = real(max(size(work, 1), size(work, 2)), dp)*epsilon(1.0_dp)*scale
      rank = 0
      row = 1
      do column = 1, size(work, 2)
         if (row > size(work, 1)) exit
         pivot_row = row - 1 + maxloc(abs(work(row:, column)), dim=1)
         pivot_size = abs(work(pivot_row, column))
         if (pivot_size <= threshold) cycle
         if (pivot_row /= row) then
            do i = column, size(work, 2)
               temporary = work(row, i)
               work(row, i) = work(pivot_row, i)
               work(pivot_row, i) = temporary
            end do
         end if
         do i = row + 1, size(work, 1)
            work(i, column:) = work(i, column:) - &
                               work(i, column)/work(row, column)*work(row, column:)
         end do
         rank = rank + 1
         row = row + 1
      end do
   end function numerical_rank

   pure subroutine evaluate_jacobian(x, parameters, model, jacobian, values, status)
      !! Evaluates an analytic Jacobian or a central finite-difference approximation.
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with shape `(n, q)`.
      real(dp), intent(in) :: parameters(:) !! Current parameter vector with size `p`.
      procedure(nls_model_interface) :: model !! Pure model-evaluation callback.
      procedure(nls_jacobian_interface), optional :: jacobian !! Optional analytic Jacobian.
      real(dp), intent(out) :: values(size(x, 1), size(parameters)) !! Model Jacobian `(n,p)`.
      integer, intent(out) :: status !! Zero on success or invalid input for nonfinite values.
      real(dp) :: lower_values(size(x, 1)), step, upper_values(size(x, 1))
      real(dp) :: shifted(size(parameters))
      integer :: j

      if (present(jacobian)) then
         call jacobian(x, parameters, values)
      else
         do j = 1, size(parameters)
            step = sqrt(epsilon(1.0_dp))*max(1.0_dp, abs(parameters(j)))
            shifted = parameters
            shifted(j) = parameters(j) + step
            call model(x, shifted, upper_values)
            shifted(j) = parameters(j) - step
            call model(x, shifted, lower_values)
            values(:, j) = (upper_values - lower_values)/(2.0_dp*step)
         end do
      end if
      status = merge(r_ok, r_invalid_input, all(ieee_is_finite(values)))
   end subroutine evaluate_jacobian

end module r_stats_nls
