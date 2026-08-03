! SPDX-License-Identifier: GPL-2.0-or-later
module fints_arima
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use fints_kinds, only : dp, pi_dp
   use fints_status, only : fints_ok, fints_invalid_input, fints_iteration_limit, &
      fints_numerical_failure
   use fints_types, only : arima_result
   use fints_linalg, only : least_squares
   use fints_summary_mod, only : sample_mean, sample_variance
   use fints_time_series, only : autocor_test
   implicit none
   private
   public :: arima_fit

contains

   subroutine arima_fit(x, order, result, seasonal_order, seasonal_period, xreg, &
      include_mean, transform_pars, method, n_cond, box_test_lag, box_test_df, &
      box_test_type, tolerance, max_iterations, initial)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: order(3)
      type(arima_result), intent(out) :: result
      integer, intent(in), optional :: seasonal_order(3), seasonal_period
      real(dp), intent(in), optional :: xreg(:,:)
      logical, intent(in), optional :: include_mean, transform_pars
      character(len=*), intent(in), optional :: method, box_test_type
      integer, intent(in), optional :: n_cond, box_test_lag, max_iterations
      real(dp), intent(in), optional :: box_test_df, tolerance, initial(:)

      integer :: p, d, q, sp, sd, sq, period, regressors, npar, n, offset
      integer :: maxit, iter, stat, likelihood_start, effective_n, lag_value, k_parameters
      real(dp) :: tol, objective_value, sse, df_value, nan_value
      real(dp), allocatable :: parameters(:), steps(:), regression_design(:,:), beta(:)
      real(dp), allocatable :: residual_ols(:), residual_work(:), fitted_work(:)
      real(dp), allocatable :: effective_ar(:), effective_ma(:), usable_residuals(:)
      real(dp) :: initial_sse, variance_x, variance_explained
      logical :: mean_flag, transform_flag, converged
      character(len=16) :: method_name, test_name

      result = arima_result()
      n = size(x)
      p = order(1)
      d = order(2)
      q = order(3)
      sp = 0
      sd = 0
      sq = 0
      if (present(seasonal_order)) then
         sp = seasonal_order(1)
         sd = seasonal_order(2)
         sq = seasonal_order(3)
      end if
      period = 1
      if (present(seasonal_period)) period = seasonal_period
      regressors = 0
      if (present(xreg)) regressors = size(xreg, 2)
      mean_flag = .true.
      if (present(include_mean)) mean_flag = include_mean
      if (d + sd > 0) mean_flag = .false.
      transform_flag = .true.
      if (present(transform_pars)) transform_flag = transform_pars
      method_name = 'CSS-ML'
      if (present(method)) method_name = uppercase(adjustl(method))
      tol = 1.0e-8_dp
      if (present(tolerance)) tol = max(tolerance, 10.0_dp * epsilon(1.0_dp))
      maxit = 2000
      if (present(max_iterations)) maxit = max_iterations

      if (n < 3 .or. any(.not. ieee_is_finite(x)) .or. any(order < 0) .or. &
         sp < 0 .or. sd < 0 .or. sq < 0 .or. maxit < 1) then
         result%status = fints_invalid_input
         return
      end if
      if ((sp + sd + sq > 0) .and. period < 2) then
         result%status = fints_invalid_input
         return
      end if
      if (present(xreg)) then
         if (size(xreg, 1) /= n .or. any(.not. ieee_is_finite(xreg))) then
            result%status = fints_invalid_input
            return
         end if
      end if
      if (trim(method_name) /= 'CSS' .and. trim(method_name) /= 'ML' .and. &
         trim(method_name) /= 'CSS-ML') then
         result%status = fints_invalid_input
         return
      end if
      offset = d + sd * period
      if (n - offset <= max(p + sp * period, q + sq * period) + 2) then
         result%status = fints_invalid_input
         return
      end if

      npar = p + q + sp + sq + merge(1, 0, mean_flag) + regressors
      allocate(parameters(npar), steps(npar))
      parameters = 0.0_dp
      steps = 0.08_dp
      call initialize_regression(x, xreg, mean_flag, parameters, p + q + sp + sq, &
         regression_design, beta, residual_ols, initial_sse)
      if (present(initial)) then
         if (size(initial) /= npar) then
            result%status = fints_invalid_input
            return
         end if
         parameters = initial
      end if
      call set_parameter_steps(parameters, p, q, sp, sq, mean_flag, regressors, &
         sqrt(max(sample_variance(x), epsilon(1.0_dp))), steps)

      if (npar > 0) then
         call nelder_mead_arima(parameters, steps, x, p, d, q, sp, sd, sq, period, &
            xreg, mean_flag, transform_flag, n_cond, tol, maxit, objective_value, &
            iter, converged)
      else
         call arima_objective(parameters, x, p, d, q, sp, sd, sq, period, xreg, &
            mean_flag, transform_flag, n_cond, objective_value)
         iter = 0
         converged = .true.
      end if
      if (.not. ieee_is_finite(objective_value) .or. objective_value >= huge(1.0_dp) / 100.0_dp) then
         result%status = fints_numerical_failure
         return
      end if

      call decode_parameters(parameters, p, q, sp, sq, mean_flag, regressors, &
         transform_flag, result%ar, result%ma, result%seasonal_ar, result%seasonal_ma, &
         result%intercept, result%regression)
      call effective_arma(result%ar, result%ma, result%seasonal_ar, result%seasonal_ma, &
         period, effective_ar, effective_ma)
      call compute_arima_residuals(x, xreg, result%intercept, result%regression, d, sd, &
         period, effective_ar, effective_ma, n_cond, residual_work, fitted_work, sse, &
         likelihood_start, stat)
      if (stat /= fints_ok) then
         result%status = stat
         return
      end if
      effective_n = size(residual_work) - likelihood_start + 1
      result%n_used = effective_n
      result%sigma2 = sse / real(effective_n, dp)
      result%log_likelihood = -0.5_dp * real(effective_n, dp) * &
         (log(2.0_dp * pi_dp * result%sigma2) + 1.0_dp)
      result%aic = -2.0_dp * result%log_likelihood + 2.0_dp * real(npar + 1, dp)
      result%method = trim(method_name)
      result%iterations = iter
      result%converged = converged
      result%status = merge(fints_ok, fints_iteration_limit, converged)
      allocate(result%coefficients(npar))
      call pack_coefficients(result%ar, result%ma, result%seasonal_ar, result%seasonal_ma, &
         mean_flag, result%intercept, result%regression, result%coefficients)

      nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate(result%residuals(n), result%fitted(n))
      result%residuals = nan_value
      result%fitted = nan_value
      result%residuals(offset + 1:n) = residual_work
      result%fitted(offset + 1:n) = x(offset + 1:n) - residual_work

      if (regressors > 0) then
         variance_x = sample_variance(x)
         variance_explained = sample_variance(matmul(xreg, result%regression))
         if (variance_x > 0.0_dp) result%r_squared = min(1.0_dp, variance_explained / variance_x)
      end if

      k_parameters = p + q + sp + sq + regressors
      lag_value = max(k_parameters + 1, nint(log(real(n, dp))))
      if (present(box_test_lag)) lag_value = max(k_parameters + 1, box_test_lag)
      lag_value = min(lag_value, effective_n - 1)
      df_value = real(max(1, lag_value - k_parameters), dp)
      if (present(box_test_df)) df_value = max(1.0_dp, box_test_df)
      test_name = 'Ljung-Box'
      if (present(box_test_type)) test_name = box_test_type
      allocate(usable_residuals(effective_n))
      usable_residuals = residual_work(likelihood_start:)
      if (lag_value >= 1) then
         call autocor_test(usable_residuals, result%box_test, lag=lag_value, &
            test_type=test_name, degrees_freedom=df_value)
      end if
   end subroutine arima_fit


   subroutine pack_coefficients(ar, ma, seasonal_ar, seasonal_ma, include_mean, intercept, &
      regression, coefficients)
      real(dp), intent(in) :: ar(:), ma(:), seasonal_ar(:), seasonal_ma(:), intercept
      logical, intent(in) :: include_mean
      real(dp), intent(in) :: regression(:)
      real(dp), intent(out) :: coefficients(:)
      integer :: position

      position = 0
      if (size(ar) > 0) then
         coefficients(position + 1:position + size(ar)) = ar
         position = position + size(ar)
      end if
      if (size(ma) > 0) then
         coefficients(position + 1:position + size(ma)) = ma
         position = position + size(ma)
      end if
      if (size(seasonal_ar) > 0) then
         coefficients(position + 1:position + size(seasonal_ar)) = seasonal_ar
         position = position + size(seasonal_ar)
      end if
      if (size(seasonal_ma) > 0) then
         coefficients(position + 1:position + size(seasonal_ma)) = seasonal_ma
         position = position + size(seasonal_ma)
      end if
      if (include_mean) then
         position = position + 1
         coefficients(position) = intercept
      end if
      if (size(regression) > 0) coefficients(position + 1:) = regression
   end subroutine pack_coefficients

   subroutine initialize_regression(x, xreg, include_mean, parameters, dynamic_count, &
      design, beta, residuals, sse)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: xreg(:,:)
      logical, intent(in) :: include_mean
      real(dp), intent(inout) :: parameters(:)
      integer, intent(in) :: dynamic_count
      real(dp), allocatable, intent(out) :: design(:,:), beta(:), residuals(:)
      real(dp), intent(out) :: sse
      integer :: columns, n, stat, position

      n = size(x)
      columns = merge(1, 0, include_mean)
      if (present(xreg)) columns = columns + size(xreg, 2)
      allocate(design(n, max(0, columns)))
      allocate(beta(max(0, columns)), residuals(n))
      residuals = x
      sse = dot_product(x, x)
      if (columns == 0) return
      position = 0
      if (include_mean) then
         position = 1
         design(:, 1) = 1.0_dp
      end if
      if (present(xreg)) design(:, position + 1:columns) = xreg
      call least_squares(design, x, beta, residuals, sse, stat)
      if (stat /= fints_ok) then
         beta = 0.0_dp
         if (include_mean) beta(1) = sample_mean(x)
      end if
      parameters(dynamic_count + 1:dynamic_count + columns) = beta
   end subroutine initialize_regression

   subroutine set_parameter_steps(parameters, p, q, sp, sq, include_mean, regressors, &
      data_scale, steps)
      real(dp), intent(in) :: parameters(:), data_scale
      integer, intent(in) :: p, q, sp, sq, regressors
      logical, intent(in) :: include_mean
      real(dp), intent(out) :: steps(:)
      integer :: dynamic_count, position, i

      dynamic_count = p + q + sp + sq
      if (dynamic_count > 0) steps(1:dynamic_count) = 0.08_dp
      position = dynamic_count
      if (include_mean) then
         position = position + 1
         steps(position) = max(0.05_dp * data_scale, 0.02_dp * abs(parameters(position)), 1.0e-4_dp)
      end if
      do i = 1, regressors
         position = position + 1
         steps(position) = max(0.05_dp * max(1.0_dp, abs(parameters(position))), 1.0e-4_dp)
      end do
   end subroutine set_parameter_steps

   subroutine nelder_mead_arima(parameters, steps, x, p, d, q, sp, sd, sq, period, &
      xreg, include_mean, transform_pars, n_cond, tolerance, max_iterations, best_value, &
      iterations, converged)
      real(dp), intent(inout) :: parameters(:)
      real(dp), intent(in) :: steps(:), x(:)
      integer, intent(in) :: p, d, q, sp, sd, sq, period
      real(dp), intent(in), optional :: xreg(:,:)
      logical, intent(in) :: include_mean, transform_pars
      integer, intent(in), optional :: n_cond
      real(dp), intent(in) :: tolerance
      integer, intent(in) :: max_iterations
      real(dp), intent(out) :: best_value
      integer, intent(out) :: iterations
      logical, intent(out) :: converged
      real(dp), allocatable :: simplex(:,:), values(:), centroid(:), reflected(:)
      real(dp), allocatable :: expanded(:), contracted(:)
      real(dp) :: reflected_value, expanded_value, contracted_value
      real(dp) :: coordinate_spread, value_spread
      integer :: npar, i, iter

      npar = size(parameters)
      allocate(simplex(npar, npar + 1), values(npar + 1), centroid(npar))
      allocate(reflected(npar), expanded(npar), contracted(npar))
      simplex(:, 1) = parameters
      do i = 1, npar
         simplex(:, i + 1) = parameters
         simplex(i, i + 1) = simplex(i, i + 1) + steps(i)
      end do
      do i = 1, npar + 1
         call arima_objective(simplex(:, i), x, p, d, q, sp, sd, sq, period, xreg, &
            include_mean, transform_pars, n_cond, values(i))
      end do

      converged = .false.
      do iter = 1, max_iterations
         call sort_simplex(simplex, values)
         coordinate_spread = maxval(abs(simplex(:, 2:) - spread(simplex(:, 1), 2, npar)))
         value_spread = maxval(abs(values(2:) - values(1)))
         if (coordinate_spread <= tolerance * (1.0_dp + maxval(abs(simplex(:, 1)))) .and. &
            value_spread <= tolerance * (1.0_dp + abs(values(1)))) then
            converged = .true.
            exit
         end if

         centroid = sum(simplex(:, 1:npar), dim=2) / real(npar, dp)
         reflected = centroid + (centroid - simplex(:, npar + 1))
         call arima_objective(reflected, x, p, d, q, sp, sd, sq, period, xreg, &
            include_mean, transform_pars, n_cond, reflected_value)

         if (reflected_value < values(1)) then
            expanded = centroid + 2.0_dp * (reflected - centroid)
            call arima_objective(expanded, x, p, d, q, sp, sd, sq, period, xreg, &
               include_mean, transform_pars, n_cond, expanded_value)
            if (expanded_value < reflected_value) then
               simplex(:, npar + 1) = expanded
               values(npar + 1) = expanded_value
            else
               simplex(:, npar + 1) = reflected
               values(npar + 1) = reflected_value
            end if
         else if (reflected_value < values(npar)) then
            simplex(:, npar + 1) = reflected
            values(npar + 1) = reflected_value
         else
            if (reflected_value < values(npar + 1)) then
               contracted = centroid + 0.5_dp * (reflected - centroid)
            else
               contracted = centroid + 0.5_dp * (simplex(:, npar + 1) - centroid)
            end if
            call arima_objective(contracted, x, p, d, q, sp, sd, sq, period, xreg, &
               include_mean, transform_pars, n_cond, contracted_value)
            if (contracted_value < min(reflected_value, values(npar + 1))) then
               simplex(:, npar + 1) = contracted
               values(npar + 1) = contracted_value
            else
               do i = 2, npar + 1
                  simplex(:, i) = simplex(:, 1) + 0.5_dp * (simplex(:, i) - simplex(:, 1))
                  call arima_objective(simplex(:, i), x, p, d, q, sp, sd, sq, period, &
                     xreg, include_mean, transform_pars, n_cond, values(i))
               end do
            end if
         end if
      end do
      call sort_simplex(simplex, values)
      parameters = simplex(:, 1)
      best_value = values(1)
      iterations = min(iter, max_iterations)
   end subroutine nelder_mead_arima

   subroutine arima_objective(parameters, x, p, d, q, sp, sd, sq, period, xreg, &
      include_mean, transform_pars, n_cond, value)
      real(dp), intent(in) :: parameters(:), x(:)
      integer, intent(in) :: p, d, q, sp, sd, sq, period
      real(dp), intent(in), optional :: xreg(:,:)
      logical, intent(in) :: include_mean, transform_pars
      integer, intent(in), optional :: n_cond
      real(dp), intent(out) :: value
      real(dp), allocatable :: ar(:), ma(:), seasonal_ar(:), seasonal_ma(:), regression(:)
      real(dp), allocatable :: effective_ar(:), effective_ma(:), residuals(:), fitted(:)
      real(dp) :: intercept, sse
      integer :: regressors, likelihood_start, stat, effective_n

      value = huge(1.0_dp)
      if (any(.not. ieee_is_finite(parameters)) .or. any(abs(parameters) > 1.0e8_dp)) return
      regressors = 0
      if (present(xreg)) regressors = size(xreg, 2)
      call decode_parameters(parameters, p, q, sp, sq, include_mean, regressors, &
         transform_pars, ar, ma, seasonal_ar, seasonal_ma, intercept, regression)
      call effective_arma(ar, ma, seasonal_ar, seasonal_ma, period, effective_ar, effective_ma)
      call compute_arima_residuals(x, xreg, intercept, regression, d, sd, period, &
         effective_ar, effective_ma, n_cond, residuals, fitted, sse, likelihood_start, stat)
      if (stat /= fints_ok .or. .not. ieee_is_finite(sse) .or. sse <= tiny(1.0_dp)) return
      effective_n = size(residuals) - likelihood_start + 1
      if (effective_n < 2) return
      value = real(effective_n, dp) * log(sse / real(effective_n, dp))
   end subroutine arima_objective

   subroutine decode_parameters(parameters, p, q, sp, sq, include_mean, regressors, &
      transform_pars, ar, ma, seasonal_ar, seasonal_ma, intercept, regression)
      real(dp), intent(in) :: parameters(:)
      integer, intent(in) :: p, q, sp, sq, regressors
      logical, intent(in) :: include_mean, transform_pars
      real(dp), allocatable, intent(out) :: ar(:), ma(:), seasonal_ar(:), seasonal_ma(:)
      real(dp), allocatable, intent(out) :: regression(:)
      real(dp), intent(out) :: intercept
      real(dp), allocatable :: reflection(:)
      integer :: position

      allocate(ar(p), ma(q), seasonal_ar(sp), seasonal_ma(sq), regression(regressors))
      position = 0
      if (p > 0) then
         if (transform_pars) then
            reflection = tanh(parameters(position + 1:position + p))
            call reflection_to_ar(reflection, ar)
         else
            ar = parameters(position + 1:position + p)
         end if
         position = position + p
      end if
      if (q > 0) then
         if (transform_pars) then
            reflection = tanh(parameters(position + 1:position + q))
            call reflection_to_ar(reflection, ma)
            ma = -ma
         else
            ma = parameters(position + 1:position + q)
         end if
         position = position + q
      end if
      if (sp > 0) then
         if (transform_pars) then
            reflection = tanh(parameters(position + 1:position + sp))
            call reflection_to_ar(reflection, seasonal_ar)
         else
            seasonal_ar = parameters(position + 1:position + sp)
         end if
         position = position + sp
      end if
      if (sq > 0) then
         if (transform_pars) then
            reflection = tanh(parameters(position + 1:position + sq))
            call reflection_to_ar(reflection, seasonal_ma)
            seasonal_ma = -seasonal_ma
         else
            seasonal_ma = parameters(position + 1:position + sq)
         end if
         position = position + sq
      end if
      intercept = 0.0_dp
      if (include_mean) then
         position = position + 1
         intercept = parameters(position)
      end if
      if (regressors > 0) regression = parameters(position + 1:position + regressors)
   end subroutine decode_parameters

   subroutine reflection_to_ar(reflection, coefficients)
      real(dp), intent(in) :: reflection(:)
      real(dp), intent(out) :: coefficients(:)
      real(dp), allocatable :: previous(:)
      integer :: order, j

      coefficients = 0.0_dp
      do order = 1, size(reflection)
         if (order > 1) then
            allocate(previous(order - 1))
            previous = coefficients(1:order - 1)
            do j = 1, order - 1
               coefficients(j) = previous(j) - reflection(order) * previous(order - j)
            end do
            deallocate(previous)
         end if
         coefficients(order) = reflection(order)
      end do
   end subroutine reflection_to_ar

   subroutine effective_arma(ar, ma, seasonal_ar, seasonal_ma, period, &
      effective_ar, effective_ma)
      real(dp), intent(in) :: ar(:), ma(:), seasonal_ar(:), seasonal_ma(:)
      integer, intent(in) :: period
      real(dp), allocatable, intent(out) :: effective_ar(:), effective_ma(:)
      real(dp), allocatable :: ar_factor(:), seasonal_ar_factor(:), ar_polynomial(:)
      real(dp), allocatable :: ma_factor(:), seasonal_ma_factor(:), ma_polynomial(:)
      integer :: i

      allocate(ar_factor(0:size(ar)), seasonal_ar_factor(0:size(seasonal_ar) * period))
      ar_factor = 0.0_dp
      seasonal_ar_factor = 0.0_dp
      ar_factor(0) = 1.0_dp
      if (size(ar) > 0) ar_factor(1:) = -ar
      seasonal_ar_factor(0) = 1.0_dp
      do i = 1, size(seasonal_ar)
         seasonal_ar_factor(i * period) = -seasonal_ar(i)
      end do
      call convolve(ar_factor, seasonal_ar_factor, ar_polynomial)
      allocate(effective_ar(ubound(ar_polynomial, 1)))
      if (size(effective_ar) > 0) effective_ar = -ar_polynomial(1:)

      allocate(ma_factor(0:size(ma)), seasonal_ma_factor(0:size(seasonal_ma) * period))
      ma_factor = 0.0_dp
      seasonal_ma_factor = 0.0_dp
      ma_factor(0) = 1.0_dp
      if (size(ma) > 0) ma_factor(1:) = ma
      seasonal_ma_factor(0) = 1.0_dp
      do i = 1, size(seasonal_ma)
         seasonal_ma_factor(i * period) = seasonal_ma(i)
      end do
      call convolve(ma_factor, seasonal_ma_factor, ma_polynomial)
      allocate(effective_ma(ubound(ma_polynomial, 1)))
      if (size(effective_ma) > 0) effective_ma = ma_polynomial(1:)
   end subroutine effective_arma

   subroutine convolve(a, b, c)
      real(dp), intent(in) :: a(0:), b(0:)
      real(dp), allocatable, intent(out) :: c(:)
      integer :: i, j

      allocate(c(0:ubound(a, 1) + ubound(b, 1)))
      c = 0.0_dp
      do i = 0, ubound(a, 1)
         do j = 0, ubound(b, 1)
            c(i + j) = c(i + j) + a(i) * b(j)
         end do
      end do
   end subroutine convolve

   subroutine compute_arima_residuals(x, xreg, intercept, regression, d, sd, period, &
      ar, ma, n_cond, residuals, fitted, sse, likelihood_start, status)
      real(dp), intent(in) :: x(:), intercept, regression(:), ar(:), ma(:)
      real(dp), intent(in), optional :: xreg(:,:)
      integer, intent(in) :: d, sd, period
      integer, intent(in), optional :: n_cond
      real(dp), allocatable, intent(out) :: residuals(:), fitted(:)
      real(dp), intent(out) :: sse
      integer, intent(out) :: likelihood_start, status
      real(dp), allocatable :: work(:), next(:)
      real(dp) :: prediction
      integer :: t, j, iteration, nwork, start_value

      allocate(work(size(x)))
      work = x - intercept
      if (present(xreg)) then
         if (size(regression) > 0) work = work - matmul(xreg, regression)
      end if
      do iteration = 1, sd
         nwork = size(work)
         if (nwork <= period) then
            status = fints_invalid_input
            allocate(residuals(0), fitted(0))
            sse = huge(1.0_dp)
            likelihood_start = 1
            return
         end if
         allocate(next(nwork - period))
         next = work(period + 1:nwork) - work(1:nwork - period)
         call move_alloc(next, work)
      end do
      do iteration = 1, d
         nwork = size(work)
         if (nwork <= 1) then
            status = fints_invalid_input
            allocate(residuals(0), fitted(0))
            sse = huge(1.0_dp)
            likelihood_start = 1
            return
         end if
         allocate(next(nwork - 1))
         next = work(2:nwork) - work(1:nwork - 1)
         call move_alloc(next, work)
      end do

      allocate(residuals(size(work)), fitted(size(work)))
      residuals = 0.0_dp
      fitted = 0.0_dp
      do t = 1, size(work)
         prediction = 0.0_dp
         do j = 1, min(size(ar), t - 1)
            prediction = prediction + ar(j) * work(t - j)
         end do
         do j = 1, min(size(ma), t - 1)
            prediction = prediction + ma(j) * residuals(t - j)
         end do
         fitted(t) = prediction
         residuals(t) = work(t) - prediction
      end do
      start_value = max(size(ar), size(ma)) + 1
      if (present(n_cond)) start_value = max(start_value, n_cond + 1)
      likelihood_start = max(1, min(start_value, size(work)))
      sse = sum(residuals(likelihood_start:) ** 2)
      if (.not. ieee_is_finite(sse)) then
         status = fints_numerical_failure
      else
         status = fints_ok
      end if
   end subroutine compute_arima_residuals

   subroutine sort_simplex(simplex, values)
      real(dp), intent(inout) :: simplex(:,:), values(:)
      real(dp) :: value_temp
      real(dp) :: point_temp(size(simplex, 1))
      integer :: i, j, k

      do i = 1, size(values) - 1
         k = i
         do j = i + 1, size(values)
            if (values(j) < values(k)) k = j
         end do
         if (k /= i) then
            value_temp = values(i)
            values(i) = values(k)
            values(k) = value_temp
            point_temp = simplex(:, i)
            simplex(:, i) = simplex(:, k)
            simplex(:, k) = point_temp
         end if
      end do
   end subroutine sort_simplex

   pure function uppercase(text) result(upper)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: upper
      integer :: i, code

      upper = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) upper(i:i) = achar(code - 32)
      end do
   end function uppercase

end module fints_arima
