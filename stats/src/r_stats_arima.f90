! SPDX-License-Identifier: MIT
! SPDX-FileComment: Conditional-sum-of-squares seasonal ARIMA fitting and forecasting.
module r_stats_arima
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_distributions, only: r_qnorm
   use r_linalg, only: inverse_matrix, signed_log_determinant, solve_system
   use r_status, only: r_ok
   use r_stats_types, only: arima_fit_t, arima_forecast_t
   use r_time_series, only: r_arma_autocovariance, r_durbin_levinson
   implicit none
   private

   real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

   public :: arima_css, arima_ml, predict_arima, predict_arima_interval

contains

   pure function arima_css(x, order, include_mean, initial, max_iterations, tolerance, &
                           seasonal_order, period, fixed_values, fixed_mask, xreg) result(fit)
      !! Fits a multiplicative seasonal ARIMA model by R-compatible conditional sum of squares.
      real(dp), intent(in) :: x(:) !! Finite time-series observations in chronological order.
      integer, intent(in) :: order(3) !! Nonnegative `(p, d, q)` model orders.
      logical, intent(in), optional :: include_mean !! Estimate a mean when `d=0`; defaults to true.
      real(dp), intent(in), optional :: initial(:) !! Initial `[ar, ma, sar, sma, mean?, xreg]` values.
      integer, intent(in), optional :: max_iterations !! Simplex iteration limit; defaults to 4000.
      real(dp), intent(in), optional :: tolerance !! Convergence tolerance; defaults to `1e-10`.
      integer, intent(in), optional :: seasonal_order(3) !! Nonnegative seasonal `(P, D, Q)` orders.
      integer, intent(in), optional :: period !! Period greater than one for seasonal terms.
      real(dp), intent(in), optional :: fixed_values(:) !! Values in `[ar, ma, sar, sma, mean?, xreg]` order.
      logical, intent(in), optional :: fixed_mask(:) !! True for each parameter held fixed.
      real(dp), intent(in), optional :: xreg(:, :) !! Regressors with one row per observation.
      type(arima_fit_t) :: fit
      real(dp), allocatable :: adjusted(:), candidate(:), full_ar(:), full_ma(:), parameters(:)
      real(dp), allocatable :: regressors(:, :), residuals(:), working(:)
      real(dp) :: best_objective, candidate_objective, objective, requested_tolerance
      integer :: arma_parameter_count, best_iterations, best_status, condition_order
      integer :: iteration_limit, parameter_count, regression_count
      integer :: p, d, q
      integer :: candidate_iterations, candidate_status
      integer :: seasonal_p, seasonal_d, seasonal_q, seasonal_period
      integer :: specified_seasonal_order(3)
      logical :: estimate_mean
      logical, allocatable :: parameter_is_fixed(:)

      p = order(1)
      d = order(2)
      q = order(3)
      specified_seasonal_order = 0
      if (present(seasonal_order)) specified_seasonal_order = seasonal_order
      seasonal_p = specified_seasonal_order(1)
      seasonal_d = specified_seasonal_order(2)
      seasonal_q = specified_seasonal_order(3)
      seasonal_period = 1
      if (present(period)) seasonal_period = period
      estimate_mean = .true.
      if (present(include_mean)) estimate_mean = include_mean
      estimate_mean = estimate_mean .and. d == 0 .and. seasonal_d == 0
      iteration_limit = 4000
      if (present(max_iterations)) iteration_limit = max_iterations
      requested_tolerance = 1.0e-10_dp
      if (present(tolerance)) requested_tolerance = tolerance
      condition_order = p + seasonal_p*seasonal_period
      if (any(order < 0) .or. any(specified_seasonal_order < 0) .or. &
          seasonal_period < 1 .or. &
          (any(specified_seasonal_order /= 0) .and. seasonal_period < 2) .or. &
          size(x) <= d + seasonal_d*seasonal_period + condition_order .or. &
          .not. all(ieee_is_finite(x)) .or. &
          iteration_limit < 1 .or. &
          requested_tolerance <= 0.0_dp) then
         fit%status = 1
         return
      end if
      regression_count = 0
      if (present(xreg)) then
         if (size(xreg, 1) /= size(x) .or. .not. all(ieee_is_finite(xreg))) then
            fit%status = 1
            return
         end if
         regression_count = size(xreg, 2)
      end if

      call difference_series(x, d, seasonal_d, seasonal_period, working)
      call difference_regressors(xreg, size(x), d, seasonal_d, seasonal_period, regressors)
      arma_parameter_count = p + q + seasonal_p + seasonal_q
      parameter_count = arma_parameter_count + merge(1, 0, estimate_mean) + regression_count
      allocate (parameters(parameter_count))
      allocate (parameter_is_fixed(parameter_count), source=.false.)
      parameters = 0.0_dp
      call initialize_regression_parameters(working, regressors, estimate_mean, &
                                            arma_parameter_count, parameters)
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            fit%status = 1
            return
         end if
         parameters = initial
      end if
      if (present(fixed_values) .neqv. present(fixed_mask)) then
         fit%status = 1
         return
      end if
      if (present(fixed_values)) then
         if (size(fixed_values) /= parameter_count .or. size(fixed_mask) /= parameter_count) then
            fit%status = 1
            return
         end if
         parameter_is_fixed = fixed_mask
         if (.not. all(ieee_is_finite(fixed_values) .or. .not. parameter_is_fixed)) then
            fit%status = 1
            return
         end if
         where (parameter_is_fixed) parameters = fixed_values
      end if

      if (count(.not. parameter_is_fixed) > 0) then
         call optimize_parameters(working, p, q, seasonal_p, seasonal_q, seasonal_period, &
                                  estimate_mean, regressors, parameters, iteration_limit, &
                                  requested_tolerance, 1, 0, 0, parameter_is_fixed, &
                                  fit%iterations, fit%status)
         if (.not. present(initial) .and. seasonal_p + seasonal_q > 0) then
            best_objective = css_objective(working, p, q, seasonal_p, seasonal_q, &
                                            seasonal_period, estimate_mean, regressors, parameters)
            best_iterations = fit%iterations
            best_status = fit%status
            allocate (candidate(parameter_count))
            candidate = 0.5_dp
            candidate(arma_parameter_count + 1:) = parameters(arma_parameter_count + 1:)
            where (parameter_is_fixed) candidate = parameters
            call optimize_parameters(working, p, q, seasonal_p, seasonal_q, seasonal_period, &
                                     estimate_mean, regressors, candidate, iteration_limit, &
                                     requested_tolerance, 1, 0, 0, parameter_is_fixed, &
                                     candidate_iterations, candidate_status)
            candidate_objective = css_objective(working, p, q, seasonal_p, seasonal_q, &
                                                 seasonal_period, estimate_mean, regressors, candidate)
            if (candidate_objective < best_objective) then
               parameters = candidate
               best_objective = candidate_objective
               best_iterations = candidate_iterations
               best_status = candidate_status
            end if
            candidate = -0.5_dp
            candidate(arma_parameter_count + 1:) = parameters(arma_parameter_count + 1:)
            where (parameter_is_fixed) candidate = parameters
            call optimize_parameters(working, p, q, seasonal_p, seasonal_q, seasonal_period, &
                                     estimate_mean, regressors, candidate, iteration_limit, &
                                     requested_tolerance, 1, 0, 0, parameter_is_fixed, &
                                     candidate_iterations, candidate_status)
            candidate_objective = css_objective(working, p, q, seasonal_p, seasonal_q, &
                                                 seasonal_period, estimate_mean, regressors, candidate)
            if (candidate_objective < best_objective) then
               parameters = candidate
               best_iterations = candidate_iterations
               best_status = candidate_status
            end if
            fit%iterations = best_iterations
            fit%status = best_status
         end if
      end if
      call expanded_coefficients(p, q, seasonal_p, seasonal_q, seasonal_period, &
                                 parameters, full_ar, full_ma)
      call regression_adjusted_series(working, regressors, estimate_mean, &
                                      arma_parameter_count, parameters, adjusted)
      call css_residuals(adjusted, condition_order, estimate_mean, &
                         arma_parameter_count, parameters, &
                         full_ar, full_ma, residuals)
      call arima_parameter_covariance(working, p, q, seasonal_p, seasonal_q, &
                                      seasonal_period, estimate_mean, regressors, parameters, &
                                      parameter_is_fixed, 1, 0, 0, fit%coef_covariance, &
                                      fit%covariance_status)
      objective = sum(residuals(condition_order + 1:)**2)
      fit%n_used = size(working)
      fit%variance = objective/real(fit%n_used - condition_order, dp)
      if (fit%variance > 0.0_dp) then
         fit%log_likelihood = -0.5_dp*real(size(working), dp)* &
                              (log(two_pi) + 1.0_dp + log(fit%variance))
      end if
      allocate (fit%ar(p), fit%ma(q), fit%seasonal_ar(seasonal_p), &
                fit%seasonal_ma(seasonal_q), fit%regression(regression_count), &
                fit%residuals(size(x)), &
                fit%innovation_variance(size(x)), fit%series(size(x)))
      allocate (fit%adjusted_series(size(x)), fit%working_series(size(working)))
      allocate (fit%fixed(parameter_count))
      if (p > 0) fit%ar = parameters(:p)
      if (q > 0) fit%ma = parameters(p + 1:p + q)
      if (seasonal_p > 0) then
         fit%seasonal_ar = parameters(p + q + 1:p + q + seasonal_p)
      end if
      if (seasonal_q > 0) then
         fit%seasonal_ma = parameters(p + q + seasonal_p + 1:p + q + seasonal_p + seasonal_q)
      end if
      if (estimate_mean) fit%mean = parameters(arma_parameter_count + 1)
      if (regression_count > 0) then
         fit%regression = parameters(arma_parameter_count + merge(1, 0, estimate_mean) + 1:)
      end if
      fit%residuals = 0.0_dp
      fit%innovation_variance = 1.0_dp
      fit%residuals(d + seasonal_d*seasonal_period + 1:) = residuals
      fit%series = x
      call regression_adjusted_series(x, xreg, estimate_mean, arma_parameter_count, &
                                      parameters, fit%adjusted_series)
      fit%working_series = adjusted
      fit%fixed = parameter_is_fixed
      fit%difference_order = d
      fit%seasonal_difference_order = seasonal_d
      fit%period = seasonal_period
      fit%n_conditioned = d + seasonal_d*seasonal_period + condition_order
      fit%include_mean = estimate_mean
   end function arima_css

   pure function arima_ml(x, order, include_mean, initial, max_iterations, tolerance, &
                          seasonal_order, period, fixed_values, fixed_mask, xreg) result(fit)
      !! Fits a multiplicative ARIMA model by exact Gaussian likelihood.
      real(dp), intent(in) :: x(:) !! Observations; exact ML permits quiet NaNs for gaps.
      integer, intent(in) :: order(3) !! Nonnegative ordinary `(p, d, q)` orders.
      logical, intent(in), optional :: include_mean !! Estimate a mean when undifferenced; defaults true.
      real(dp), intent(in), optional :: initial(:) !! Initial `[ar, ma, sar, sma, mean?, xreg]` values.
      integer, intent(in), optional :: max_iterations !! Simplex iteration limit; defaults to 4000.
      real(dp), intent(in), optional :: tolerance !! Convergence tolerance; defaults to `1e-10`.
      integer, intent(in), optional :: seasonal_order(3) !! Nonnegative seasonal `(P, D, Q)` orders.
      integer, intent(in), optional :: period !! Period greater than one for seasonal terms.
      real(dp), intent(in), optional :: fixed_values(:) !! Values in `[ar, ma, sar, sma, mean?, xreg]` order.
      logical, intent(in), optional :: fixed_mask(:) !! True for each parameter held fixed.
      real(dp), intent(in), optional :: xreg(:, :) !! Regressors with one row per observation.
      type(arima_fit_t) :: fit
      real(dp), allocatable :: adjusted(:), adjusted_original(:), exact_residuals(:)
      real(dp), allocatable :: full_ar(:), full_ma(:), objective_regressors(:, :)
      real(dp), allocatable :: objective_series(:), parameters(:), regressors(:, :)
      real(dp), allocatable :: variance_factor(:), working(:)
      real(dp) :: mean_value, requested_tolerance
      integer :: arma_parameter_count, d, iteration_limit, likelihood_status, loss
      integer :: parameter_count, p, q, regression_count
      integer :: seasonal_d, seasonal_p, seasonal_period, seasonal_q
      integer :: specified_seasonal_order(3)
      logical :: estimate_mean, use_diffuse_missing
      logical, allocatable :: parameter_is_fixed(:)
      type(arima_fit_t) :: css_fit

      p = order(1)
      d = order(2)
      q = order(3)
      specified_seasonal_order = 0
      if (present(seasonal_order)) specified_seasonal_order = seasonal_order
      seasonal_p = specified_seasonal_order(1)
      seasonal_d = specified_seasonal_order(2)
      seasonal_q = specified_seasonal_order(3)
      seasonal_period = 1
      if (present(period)) seasonal_period = period
      loss = d + seasonal_d*seasonal_period
      estimate_mean = .true.
      if (present(include_mean)) estimate_mean = include_mean
      estimate_mean = estimate_mean .and. loss == 0
      iteration_limit = 4000
      if (present(max_iterations)) iteration_limit = max_iterations
      requested_tolerance = 1.0e-10_dp
      if (present(tolerance)) requested_tolerance = tolerance
      if (any(order < 0) .or. any(specified_seasonal_order < 0) .or. &
          seasonal_period < 1 .or. &
          (any(specified_seasonal_order /= 0) .and. seasonal_period < 2) .or. &
          size(x) - loss < 2 .or. &
          any(.not. ieee_is_finite(x) .and. .not. ieee_is_nan(x)) .or. &
          iteration_limit < 1 .or. &
          requested_tolerance <= 0.0_dp) then
         fit%status = 1
         return
      end if
      regression_count = 0
      if (present(xreg)) then
         if (size(xreg, 1) /= size(x) .or. .not. all(ieee_is_finite(xreg))) then
            fit%status = 1
            return
         end if
         regression_count = size(xreg, 2)
      end if

      call difference_series(x, d, seasonal_d, seasonal_period, working)
      use_diffuse_missing = loss > 0 .and. any(ieee_is_nan(x))
      if (count(ieee_is_finite(x)) - loss < 1) then
         fit%status = 1
         return
      end if
      call difference_regressors(xreg, size(x), d, seasonal_d, seasonal_period, regressors)
      if (use_diffuse_missing) then
         objective_series = x
         call difference_regressors(xreg, size(x), 0, 0, 1, objective_regressors)
      else
         objective_series = working
         objective_regressors = regressors
      end if
      arma_parameter_count = p + q + seasonal_p + seasonal_q
      parameter_count = arma_parameter_count + merge(1, 0, estimate_mean) + regression_count
      allocate (parameters(parameter_count), source=0.0_dp)
      allocate (parameter_is_fixed(parameter_count), source=.false.)
      call initialize_regression_parameters(working, regressors, estimate_mean, &
                                            arma_parameter_count, parameters)
      if (present(initial)) then
         if (size(initial) /= parameter_count) then
            fit%status = 1
            return
         end if
         parameters = initial
      end if
      if (present(fixed_values) .neqv. present(fixed_mask)) then
         fit%status = 1
         return
      end if
      if (present(fixed_values)) then
         if (size(fixed_values) /= parameter_count .or. size(fixed_mask) /= parameter_count) then
            fit%status = 1
            return
         end if
         parameter_is_fixed = fixed_mask
         if (.not. all(ieee_is_finite(fixed_values) .or. .not. parameter_is_fixed)) then
            fit%status = 1
            return
         end if
         where (parameter_is_fixed) parameters = fixed_values
      end if
      if (.not. present(initial) .and. .not. any(ieee_is_nan(x))) then
         css_fit = arima_css(x, order, include_mean=estimate_mean, &
                             max_iterations=iteration_limit, tolerance=requested_tolerance, &
                             seasonal_order=specified_seasonal_order, period=seasonal_period, &
                             fixed_values=parameters, fixed_mask=parameter_is_fixed, xreg=xreg)
         if (css_fit%status /= 0) then
            fit%status = css_fit%status
            return
         end if
         if (p > 0) parameters(:p) = css_fit%ar
         if (q > 0) parameters(p + 1:p + q) = css_fit%ma
         if (seasonal_p > 0) then
            parameters(p + q + 1:p + q + seasonal_p) = css_fit%seasonal_ar
         end if
         if (seasonal_q > 0) then
            parameters(p + q + seasonal_p + 1:p + q + seasonal_p + seasonal_q) = &
               css_fit%seasonal_ma
         end if
         if (estimate_mean) parameters(arma_parameter_count + 1) = css_fit%mean
         if (regression_count > 0) then
            parameters(arma_parameter_count + merge(1, 0, estimate_mean) + 1:) = &
               css_fit%regression
         end if
      end if

      if (count(.not. parameter_is_fixed) > 0) then
         call optimize_parameters(objective_series, p, q, seasonal_p, seasonal_q, &
                                  seasonal_period, estimate_mean, objective_regressors, &
                                  parameters, iteration_limit, requested_tolerance, &
                                  merge(3, 2, use_diffuse_missing), &
                                  merge(d, 0, use_diffuse_missing), &
                                  merge(seasonal_d, 0, use_diffuse_missing), &
                                  parameter_is_fixed, &
                                  fit%iterations, fit%status)
      end if
      call expanded_coefficients(p, q, seasonal_p, seasonal_q, seasonal_period, &
                                 parameters, full_ar, full_ma)
      call regression_adjusted_series(working, regressors, estimate_mean, &
                                      arma_parameter_count, parameters, adjusted)
      call regression_adjusted_series(x, xreg, estimate_mean, arma_parameter_count, &
                                      parameters, adjusted_original)
      mean_value = 0.0_dp
      if (estimate_mean) mean_value = parameters(arma_parameter_count + 1)
      if (use_diffuse_missing) then
         call exact_diffuse_arma_likelihood(adjusted_original, full_ar, full_ma, d, &
                                            seasonal_d, seasonal_period, fit%log_likelihood, &
                                            fit%variance, exact_residuals, variance_factor, &
                                            likelihood_status)
      else
         call exact_arma_likelihood(adjusted, full_ar, full_ma, mean_value, &
                                    fit%log_likelihood, fit%variance, exact_residuals, &
                                    variance_factor, likelihood_status)
      end if
      if (likelihood_status /= r_ok) then
         fit%status = 1
         return
      end if
      call arima_parameter_covariance(objective_series, p, q, seasonal_p, seasonal_q, &
                                      seasonal_period, estimate_mean, objective_regressors, &
                                      parameters, parameter_is_fixed, &
                                      merge(3, 2, use_diffuse_missing), &
                                      merge(d, 0, use_diffuse_missing), &
                                      merge(seasonal_d, 0, use_diffuse_missing), &
                                      fit%coef_covariance, &
                                      fit%covariance_status)
      allocate (fit%ar(p), fit%ma(q), fit%seasonal_ar(seasonal_p), &
                fit%seasonal_ma(seasonal_q), fit%regression(regression_count))
      allocate (fit%residuals(size(x)), fit%innovation_variance(size(x)))
      allocate (fit%series(size(x)), fit%adjusted_series(size(x)))
      allocate (fit%working_series(size(working)))
      allocate (fit%fixed(parameter_count))
      if (p > 0) fit%ar = parameters(:p)
      if (q > 0) fit%ma = parameters(p + 1:p + q)
      if (seasonal_p > 0) then
         fit%seasonal_ar = parameters(p + q + 1:p + q + seasonal_p)
      end if
      if (seasonal_q > 0) then
         fit%seasonal_ma = parameters(p + q + seasonal_p + 1:p + q + seasonal_p + seasonal_q)
      end if
      if (estimate_mean) fit%mean = parameters(arma_parameter_count + 1)
      if (regression_count > 0) then
         fit%regression = parameters(arma_parameter_count + merge(1, 0, estimate_mean) + 1:)
      end if
      if (use_diffuse_missing) then
         fit%residuals = exact_residuals/sqrt(variance_factor)
         fit%innovation_variance = variance_factor
      else
         fit%residuals = 0.0_dp
         fit%innovation_variance = 1.0_dp
         fit%residuals(loss + 1:) = exact_residuals/sqrt(variance_factor)
         fit%innovation_variance(loss + 1:) = variance_factor
      end if
      fit%series = x
      fit%adjusted_series = adjusted_original
      fit%working_series = adjusted
      fit%fixed = parameter_is_fixed
      fit%difference_order = d
      fit%seasonal_difference_order = seasonal_d
      fit%period = seasonal_period
      fit%n_conditioned = loss
      if (use_diffuse_missing) then
         fit%n_used = count(ieee_is_finite(x)) - loss
      else
         fit%n_used = count(ieee_is_finite(working))
      end if
      fit%include_mean = estimate_mean
      fit%method = "ML"
   end function arima_ml

   pure function predict_arima(fit, n_ahead, new_xreg) result(prediction)
      !! Computes deterministic point forecasts from a fitted CSS or exact ARIMA model.
      type(arima_fit_t), intent(in) :: fit !! Successful fit containing histories and coefficients.
      integer, intent(in) :: n_ahead !! Positive number of future observations requested.
      real(dp), intent(in), optional :: new_xreg(:, :) !! Future regressors with shape `(n_ahead, k)`.
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: difference_weights(:), errors(:), extended_series(:)
      real(dp), allocatable :: full_ar(:), full_ma(:), ml_forecast(:), parameters(:), working(:)
      real(dp) :: next_value
      integer :: d, h, i, loss, n, p, q, seasonal_d, seasonal_p, seasonal_q, status

      if (fit%status /= 0 .or. n_ahead < 1 .or. .not. allocated(fit%adjusted_series) .or. &
          .not. allocated(fit%working_series)) then
         allocate (prediction(0))
         return
      end if
      if (size(fit%regression) > 0) then
         if (.not. present(new_xreg)) then
            allocate (prediction(0))
            return
         end if
         if (size(new_xreg, 1) /= n_ahead .or. &
             size(new_xreg, 2) /= size(fit%regression) .or. &
             .not. all(ieee_is_finite(new_xreg))) then
            allocate (prediction(0))
            return
         end if
      end if
      p = size(fit%ar)
      q = size(fit%ma)
      seasonal_p = size(fit%seasonal_ar)
      seasonal_q = size(fit%seasonal_ma)
      d = fit%difference_order
      seasonal_d = fit%seasonal_difference_order
      loss = d + seasonal_d*fit%period
      n = size(fit%working_series)
      allocate (prediction(n_ahead), working(n + n_ahead), errors(n + n_ahead))
      allocate (extended_series(size(fit%adjusted_series) + n_ahead))
      allocate (parameters(p + q + seasonal_p + seasonal_q))
      if (p > 0) parameters(:p) = fit%ar
      if (q > 0) parameters(p + 1:p + q) = fit%ma
      if (seasonal_p > 0) parameters(p + q + 1:p + q + seasonal_p) = fit%seasonal_ar
      if (seasonal_q > 0) then
         parameters(p + q + seasonal_p + 1:) = fit%seasonal_ma
      end if
      call expanded_coefficients(p, q, seasonal_p, seasonal_q, fit%period, &
                                 parameters, full_ar, full_ma)
      call difference_weights_for(d, seasonal_d, fit%period, difference_weights)
      working(:n) = fit%working_series
      errors = 0.0_dp
      errors(:n) = fit%residuals(loss + 1:)
      if (trim(fit%method) == "ML") then
         errors(:n) = errors(:n)*sqrt(fit%innovation_variance(loss + 1:))
         call exact_arma_forecast(fit%working_series, full_ar, full_ma, fit%mean, &
                                  n_ahead, ml_forecast, status)
         if (status /= r_ok) then
            deallocate (prediction)
            allocate (prediction(0))
            return
         end if
      end if
      extended_series(:size(fit%adjusted_series)) = fit%adjusted_series

      do h = 1, n_ahead
         if (trim(fit%method) == "ML") then
            next_value = ml_forecast(h)
         else
            next_value = fit%mean
            do i = 1, min(size(full_ar), n + h - 1)
               next_value = next_value + full_ar(i)*(working(n + h - i) - fit%mean)
            end do
            do i = 1, min(size(full_ma), n + h - 1)
               next_value = next_value + full_ma(i)*errors(n + h - i)
            end do
         end if
         working(n + h) = next_value
         if (loss == 0) then
            prediction(h) = next_value
         else
            do i = 1, min(size(difference_weights), size(fit%adjusted_series) + h - 1)
               next_value = next_value + difference_weights(i)* &
                            extended_series(size(fit%adjusted_series) + h - i)
            end do
            prediction(h) = next_value
         end if
         extended_series(size(fit%adjusted_series) + h) = prediction(h)
         if (size(fit%regression) > 0) then
            prediction(h) = prediction(h) + dot_product(new_xreg(h, :), fit%regression)
         end if
      end do
   end function predict_arima

   pure function predict_arima_interval(fit, n_ahead, confidence_level, new_xreg) result(forecast)
      !! Computes ARIMA point forecasts and innovation-based normal prediction intervals.
      type(arima_fit_t), intent(in) :: fit !! Successful CSS or ML fit.
      integer, intent(in) :: n_ahead !! Positive number of future observations requested.
      real(dp), intent(in), optional :: confidence_level !! Central probability; defaults to 0.95.
      real(dp), intent(in), optional :: new_xreg(:, :) !! Future regressors with shape `(n_ahead, k)`.
      type(arima_forecast_t) :: forecast
      real(dp), allocatable :: impulse(:)
      real(dp) :: level, multiplier, variance_sum
      integer :: h

      level = 0.95_dp
      if (present(confidence_level)) level = confidence_level
      if (fit%status /= 0 .or. n_ahead < 1 .or. level <= 0.0_dp .or. level >= 1.0_dp) then
         forecast%status = 1
         return
      end if
      if (present(new_xreg)) then
         forecast%fit = predict_arima(fit, n_ahead, new_xreg)
      else
         forecast%fit = predict_arima(fit, n_ahead)
      end if
      if (size(forecast%fit) /= n_ahead) then
         forecast%status = 1
         return
      end if
      call arima_impulse_weights(fit, n_ahead, impulse)
      allocate (forecast%standard_error(n_ahead), forecast%lower(n_ahead), &
                forecast%upper(n_ahead))
      variance_sum = 0.0_dp
      do h = 1, n_ahead
         variance_sum = variance_sum + impulse(h - 1)**2
         forecast%standard_error(h) = sqrt(max(0.0_dp, fit%variance*variance_sum))
      end do
      forecast%confidence_level = level
      multiplier = r_qnorm(0.5_dp*(1.0_dp + level))
      forecast%lower = forecast%fit - multiplier*forecast%standard_error
      forecast%upper = forecast%fit + multiplier*forecast%standard_error
   end function predict_arima_interval

   pure subroutine arima_impulse_weights(fit, n_ahead, weights)
      !! Generates the impulse response of a fitted multiplicative ARIMA model.
      type(arima_fit_t), intent(in) :: fit !! Fitted ARIMA coefficients and differencing orders.
      integer, intent(in) :: n_ahead !! Number of impulse coefficients requested.
      real(dp), allocatable, intent(out) :: weights(:) !! Coefficients indexed from lag zero.
      real(dp), allocatable :: ar_polynomial(:), denominator(:), difference_polynomial(:)
      real(dp), allocatable :: difference_weights(:), full_ar(:), full_ma(:), parameters(:)
      real(dp), allocatable :: recurrence(:)
      integer :: h, j, p, q, seasonal_p, seasonal_q

      p = size(fit%ar)
      q = size(fit%ma)
      seasonal_p = size(fit%seasonal_ar)
      seasonal_q = size(fit%seasonal_ma)
      allocate (parameters(p + q + seasonal_p + seasonal_q))
      if (p > 0) parameters(:p) = fit%ar
      if (q > 0) parameters(p + 1:p + q) = fit%ma
      if (seasonal_p > 0) parameters(p + q + 1:p + q + seasonal_p) = fit%seasonal_ar
      if (seasonal_q > 0) parameters(p + q + seasonal_p + 1:) = fit%seasonal_ma
      call expanded_coefficients(p, q, seasonal_p, seasonal_q, fit%period, &
                                 parameters, full_ar, full_ma)
      call difference_weights_for(fit%difference_order, fit%seasonal_difference_order, &
                                  fit%period, difference_weights)

      allocate (ar_polynomial(0:size(full_ar)), source=0.0_dp)
      allocate (difference_polynomial(0:size(difference_weights)), source=0.0_dp)
      ar_polynomial(0) = 1.0_dp
      difference_polynomial(0) = 1.0_dp
      if (size(full_ar) > 0) ar_polynomial(1:) = -full_ar
      if (size(difference_weights) > 0) difference_polynomial(1:) = -difference_weights
      call multiply_polynomials(ar_polynomial, difference_polynomial, denominator)
      allocate (recurrence(size(denominator) - 1))
      if (size(recurrence) > 0) recurrence = -denominator(1:)

      allocate (weights(0:n_ahead - 1), source=0.0_dp)
      weights(0) = 1.0_dp
      do h = 1, n_ahead - 1
         if (h <= size(full_ma)) weights(h) = full_ma(h)
         do j = 1, min(h, size(recurrence))
            weights(h) = weights(h) + recurrence(j)*weights(h - j)
         end do
      end do
   end subroutine arima_impulse_weights

   pure subroutine exact_arma_forecast(x, ar, ma, mean_value, n_ahead, forecast, status)
      !! Computes exact finite-history Gaussian ARMA point forecasts.
      real(dp), intent(in) :: x(:) !! Stationary observation history.
      real(dp), intent(in) :: ar(:) !! Effective autoregressive coefficients by lag.
      real(dp), intent(in) :: ma(:) !! Effective moving-average coefficients by lag.
      real(dp), intent(in) :: mean_value !! Stationary process mean.
      integer, intent(in) :: n_ahead !! Number of future observations requested.
      real(dp), allocatable, intent(out) :: forecast(:) !! Exact conditional point forecasts.
      integer, intent(out) :: status !! Zero on success or a shared numerical status.
      real(dp), allocatable :: ar_table(:, :), autocovariance(:), extended(:), partial(:)
      real(dp), allocatable :: prediction_variance(:)
      integer :: h, j, k, n

      n = size(x)
      if (any(ieee_is_nan(x))) then
         call exact_arma_forecast_missing(x, ar, ma, mean_value, n_ahead, forecast, status)
         return
      end if
      call r_arma_autocovariance(ar, ma, n + n_ahead - 1, autocovariance, status)
      if (status /= r_ok) return
      call r_durbin_levinson(autocovariance, partial, status, ar_table, prediction_variance)
      if (status /= r_ok) return
      allocate (forecast(n_ahead), extended(n + n_ahead))
      extended(:n) = x
      do h = 1, n_ahead
         k = n + h - 1
         forecast(h) = mean_value
         do j = 1, k
            forecast(h) = forecast(h) + ar_table(k, j)*(extended(n + h - j) - mean_value)
         end do
         extended(n + h) = forecast(h)
      end do
   end subroutine exact_arma_forecast

   pure subroutine exact_arma_forecast_missing(x, ar, ma, mean_value, n_ahead, forecast, status)
      !! Conditions stationary Gaussian forecasts on all available historical observations.
      real(dp), intent(in) :: x(:) !! Stationary history containing quiet NaNs for gaps.
      real(dp), intent(in) :: ar(:) !! Effective autoregressive coefficients by lag.
      real(dp), intent(in) :: ma(:) !! Effective moving-average coefficients by lag.
      real(dp), intent(in) :: mean_value !! Stationary process mean.
      integer, intent(in) :: n_ahead !! Number of future observations requested.
      real(dp), allocatable, intent(out) :: forecast(:) !! Conditional point forecasts.
      integer, intent(out) :: status !! Zero on success or one for invalid covariance.
      real(dp), allocatable :: autocovariance(:), centered(:), covariance(:, :), solution(:)
      integer, allocatable :: observed_indices(:)
      integer :: h, i, info, j, k, n, observed_count

      n = size(x)
      observed_count = count(ieee_is_finite(x))
      if (observed_count < 1) then
         allocate (forecast(0))
         status = 1
         return
      end if
      allocate (observed_indices(observed_count), centered(observed_count))
      k = 0
      do i = 1, n
         if (.not. ieee_is_finite(x(i))) cycle
         k = k + 1
         observed_indices(k) = i
         centered(k) = x(i) - mean_value
      end do
      call r_arma_autocovariance(ar, ma, n + n_ahead - 1, autocovariance, status)
      if (status /= r_ok) return
      allocate (covariance(observed_count, observed_count), solution(observed_count))
      do j = 1, observed_count
         do i = 1, observed_count
            covariance(i, j) = autocovariance(abs(observed_indices(i) - observed_indices(j)))
         end do
      end do
      call solve_system(covariance, centered, solution, info)
      if (info /= 0) then
         status = 1
         return
      end if
      allocate (forecast(n_ahead))
      do h = 1, n_ahead
         forecast(h) = mean_value
         do j = 1, observed_count
            forecast(h) = forecast(h) + &
               autocovariance(n + h - observed_indices(j))*solution(j)
         end do
      end do
      status = r_ok
   end subroutine exact_arma_forecast_missing

   pure subroutine difference_series(x, difference_order, seasonal_difference_order, period, values)
      !! Applies repeated ordinary and seasonal differences to a time series.
      real(dp), intent(in) :: x(:) !! Original observations.
      integer, intent(in) :: difference_order !! Number of first differences.
      integer, intent(in) :: seasonal_difference_order !! Number of seasonal differences.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), allocatable, intent(out) :: values(:) !! Differenced observations.
      real(dp), allocatable :: temporary(:)
      integer :: k

      values = x
      do k = 1, difference_order
         temporary = values(2:) - values(:size(values) - 1)
         call move_alloc(temporary, values)
      end do
      do k = 1, seasonal_difference_order
         temporary = values(period + 1:) - values(:size(values) - period)
         call move_alloc(temporary, values)
      end do
   end subroutine difference_series

   pure subroutine difference_regressors(xreg, observation_count, difference_order, &
                                         seasonal_difference_order, period, values)
      !! Applies the response differencing operators to every regressor column.
      real(dp), intent(in), optional :: xreg(:, :) !! Original regressor matrix.
      integer, intent(in) :: observation_count !! Number of original response observations.
      integer, intent(in) :: difference_order !! Number of ordinary differences.
      integer, intent(in) :: seasonal_difference_order !! Number of seasonal differences.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), allocatable, intent(out) :: values(:, :) !! Differenced regressor matrix.
      real(dp), allocatable :: column(:)
      integer :: j, retained_count

      retained_count = observation_count - difference_order - seasonal_difference_order*period
      if (.not. present(xreg)) then
         allocate (values(retained_count, 0))
         return
      end if
      allocate (values(retained_count, size(xreg, 2)))
      do j = 1, size(xreg, 2)
         call difference_series(xreg(:, j), difference_order, seasonal_difference_order, &
                                period, column)
         values(:, j) = column
      end do
   end subroutine difference_regressors

   pure subroutine initialize_regression_parameters(x, regressors, estimate_mean, &
                                                    arma_parameter_count, parameters)
      !! Initializes the mean and regression coefficients by ordinary least squares.
      real(dp), intent(in) :: x(:) !! Differenced response observations.
      real(dp), intent(in) :: regressors(:, :) !! Differenced user regressors.
      logical, intent(in) :: estimate_mean !! Whether an intercept-like mean is included.
      integer, intent(in) :: arma_parameter_count !! Number of preceding AR and MA parameters.
      real(dp), intent(inout) :: parameters(:) !! Parameter vector receiving initial estimates.
      real(dp), allocatable :: design(:, :), full_design(:, :), normal_matrix(:, :)
      real(dp), allocatable :: rhs(:), solution(:)
      logical, allocatable :: observed(:)
      integer :: deterministic_count, info, j, mean_count, observed_count

      mean_count = merge(1, 0, estimate_mean)
      deterministic_count = mean_count + size(regressors, 2)
      if (deterministic_count == 0) return
      observed = ieee_is_finite(x)
      observed_count = count(observed)
      if (observed_count == 0) return
      allocate (full_design(size(x), deterministic_count))
      allocate (design(observed_count, deterministic_count), rhs(deterministic_count))
      allocate (solution(deterministic_count))
      if (estimate_mean) full_design(:, 1) = 1.0_dp
      if (size(regressors, 2) > 0) full_design(:, mean_count + 1:) = regressors
      do j = 1, deterministic_count
         design(:, j) = pack(full_design(:, j), observed)
      end do
      normal_matrix = matmul(transpose(design), design)
      rhs = matmul(transpose(design), pack(x, observed))
      call solve_system(normal_matrix, rhs, solution, info)
      if (info == 0) then
         parameters(arma_parameter_count + 1:) = solution
      else if (estimate_mean) then
         parameters(arma_parameter_count + 1) = &
            sum(pack(x, observed))/real(observed_count, dp)
      end if
   end subroutine initialize_regression_parameters

   pure subroutine regression_adjusted_series(x, regressors, estimate_mean, &
                                              arma_parameter_count, parameters, adjusted)
      !! Removes user-regressor effects while retaining the separately modeled mean.
      real(dp), intent(in) :: x(:) !! Original or differenced response observations.
      real(dp), intent(in), optional :: regressors(:, :) !! Conformable user regressors.
      logical, intent(in) :: estimate_mean !! Whether a mean precedes user coefficients.
      integer, intent(in) :: arma_parameter_count !! Number of preceding AR and MA parameters.
      real(dp), intent(in) :: parameters(:) !! Full candidate parameter vector.
      real(dp), allocatable, intent(out) :: adjusted(:) !! Response minus regression effects.
      integer :: coefficient_start

      adjusted = x
      if (.not. present(regressors)) return
      if (size(regressors, 2) == 0) return
      coefficient_start = arma_parameter_count + merge(1, 0, estimate_mean) + 1
      adjusted = x - matmul(regressors, parameters(coefficient_start:))
   end subroutine regression_adjusted_series

   pure subroutine css_residuals(x, condition_order, estimate_mean, arma_parameter_count, parameters, &
                                 full_ar, full_ma, residuals)
      !! Evaluates conditional ARMA innovations for one parameter vector.
      real(dp), intent(in) :: x(:) !! Stationary or differenced series.
      integer, intent(in) :: condition_order !! Number of leading values excluded from CSS.
      logical, intent(in) :: estimate_mean !! Whether the final parameter is the series mean.
      integer, intent(in) :: arma_parameter_count !! Number of AR and MA parameters.
      real(dp), intent(in) :: parameters(:) !! Values ordered as `[ar, ma, sar, sma, mean?, xreg]`.
      real(dp), intent(in) :: full_ar(:) !! Expanded multiplicative AR coefficients by lag.
      real(dp), intent(in) :: full_ma(:) !! Expanded multiplicative MA coefficients by lag.
      real(dp), allocatable, intent(out) :: residuals(:) !! Conditional innovations.
      real(dp) :: center
      integer :: i, j

      allocate (residuals(size(x)), source=0.0_dp)
      center = 0.0_dp
      if (estimate_mean) center = parameters(arma_parameter_count + 1)
      do i = condition_order + 1, size(x)
         residuals(i) = x(i) - center
         do j = 1, min(size(full_ar), i - 1)
            residuals(i) = residuals(i) - full_ar(j)*(x(i - j) - center)
         end do
         do j = 1, min(size(full_ma), i - 1)
            residuals(i) = residuals(i) - full_ma(j)*residuals(i - j)
         end do
      end do
   end subroutine css_residuals

   pure function css_objective(x, p, q, seasonal_p, seasonal_q, period, &
                               estimate_mean, regressors, parameters) result(value)
      !! Returns the logarithmic conditional residual-variance objective.
      real(dp), intent(in) :: x(:) !! Stationary or differenced series.
      integer, intent(in) :: p !! Autoregressive order.
      integer, intent(in) :: q !! Moving-average order.
      integer, intent(in) :: seasonal_p !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_q !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      logical, intent(in) :: estimate_mean !! Whether the final parameter is the series mean.
      real(dp), intent(in) :: regressors(:, :) !! Differenced user regressors.
      real(dp), intent(in) :: parameters(:) !! Candidate `[ar, ma, sar, sma, mean?, xreg]` values.
      real(dp) :: value
      real(dp), allocatable :: adjusted(:), full_ar(:), full_ma(:), residuals(:)
      integer :: arma_parameter_count, condition_order

      arma_parameter_count = p + q + seasonal_p + seasonal_q
      condition_order = p + seasonal_p*period
      call expanded_coefficients(p, q, seasonal_p, seasonal_q, period, &
                                 parameters, full_ar, full_ma)
      call regression_adjusted_series(x, regressors, estimate_mean, &
                                      arma_parameter_count, parameters, adjusted)
      call css_residuals(adjusted, condition_order, estimate_mean, &
                         arma_parameter_count, parameters, &
                         full_ar, full_ma, residuals)
      value = sum(residuals(condition_order + 1:)**2)/real(size(x) - condition_order, dp)
      value = 0.5_dp*log(max(value, tiny(1.0_dp)))
   end function css_objective

   pure function arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                 estimate_mean, regressors, parameters, objective_method, &
                                 difference_order, seasonal_difference_order) result(value)
      !! Evaluates either the CSS or profiled exact Gaussian optimization objective.
      real(dp), intent(in) :: x(:) !! Stationary observations supplied to the objective.
      integer, intent(in) :: p !! Ordinary autoregressive order.
      integer, intent(in) :: q !! Ordinary moving-average order.
      integer, intent(in) :: seasonal_p !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_q !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      logical, intent(in) :: estimate_mean !! Whether the final parameter is a mean.
      real(dp), intent(in) :: regressors(:, :) !! Differenced user regressors.
      real(dp), intent(in) :: parameters(:) !! Candidate model parameters.
      integer, intent(in) :: objective_method !! One for CSS, two for stationary ML, or three for diffuse ML.
      integer, intent(in) :: difference_order !! Diffuse ordinary differencing order.
      integer, intent(in) :: seasonal_difference_order !! Diffuse seasonal differencing order.
      real(dp) :: value
      real(dp), allocatable :: adjusted(:), full_ar(:), full_ma(:), residuals(:)
      real(dp), allocatable :: variance_factor(:)
      real(dp) :: log_likelihood, mean_value, variance
      integer :: arma_parameter_count, status

      arma_parameter_count = p + q + seasonal_p + seasonal_q
      if (objective_method == 1) then
         value = css_objective(x, p, q, seasonal_p, seasonal_q, period, &
                               estimate_mean, regressors, parameters)
         return
      end if
      call expanded_coefficients(p, q, seasonal_p, seasonal_q, period, &
                                 parameters, full_ar, full_ma)
      mean_value = 0.0_dp
      if (estimate_mean) mean_value = parameters(arma_parameter_count + 1)
      call regression_adjusted_series(x, regressors, estimate_mean, &
                                      arma_parameter_count, parameters, adjusted)
      if (objective_method == 3) then
         call exact_diffuse_arma_likelihood(adjusted, full_ar, full_ma, difference_order, &
                                            seasonal_difference_order, period, log_likelihood, &
                                            variance, residuals, variance_factor, status)
      else
         call exact_arma_likelihood(adjusted, full_ar, full_ma, mean_value, log_likelihood, &
                                    variance, residuals, variance_factor, status)
      end if
      if (status == r_ok .and. ieee_is_finite(log_likelihood)) then
         value = -log_likelihood/real(count(ieee_is_finite(adjusted)) - &
                 difference_order - seasonal_difference_order*period, dp)
      else
         value = huge(1.0_dp)
      end if
   end function arima_objective

   pure subroutine exact_arma_likelihood(x, ar, ma, mean_value, log_likelihood, variance, &
                                         residuals, variance_factor, status)
      !! Evaluates the profiled exact Gaussian likelihood through innovations recursions.
      real(dp), intent(in) :: x(:) !! Stationary observations, optionally containing quiet NaNs.
      real(dp), intent(in) :: ar(:) !! Effective autoregressive coefficients by lag.
      real(dp), intent(in) :: ma(:) !! Effective moving-average coefficients by lag.
      real(dp), intent(in) :: mean_value !! Stationary process mean.
      real(dp), intent(out) :: log_likelihood !! Profiled Gaussian log likelihood.
      real(dp), intent(out) :: variance !! Profiled innovation variance.
      real(dp), allocatable, intent(out) :: residuals(:) !! Exact one-step innovations.
      real(dp), allocatable, intent(out) :: variance_factor(:) !! Unit-scale innovation variances.
      integer, intent(out) :: status !! Zero on success or a shared numerical status.
      real(dp), allocatable :: ar_table(:, :), autocovariance(:), partial(:), prediction_variance(:)
      real(dp) :: sum_squares
      integer :: i, j, n

      n = size(x)
      log_likelihood = -huge(1.0_dp)
      variance = huge(1.0_dp)
      if (any(ieee_is_nan(x))) then
         call exact_arma_likelihood_missing(x, ar, ma, mean_value, log_likelihood, &
                                            variance, residuals, variance_factor, status)
         return
      end if
      call r_arma_autocovariance(ar, ma, n - 1, autocovariance, status)
      if (status /= r_ok) return
      call r_durbin_levinson(autocovariance, partial, status, ar_table, prediction_variance)
      if (status /= r_ok) return
      allocate (residuals(n), variance_factor(n))
      residuals(1) = x(1) - mean_value
      variance_factor(1) = prediction_variance(0)
      do i = 2, n
         residuals(i) = x(i) - mean_value
         do j = 1, i - 1
            residuals(i) = residuals(i) - ar_table(i - 1, j)*(x(i - j) - mean_value)
         end do
         variance_factor(i) = prediction_variance(i - 1)
      end do
      if (any(variance_factor <= 0.0_dp)) then
         status = 1
         return
      end if
      sum_squares = sum(residuals**2/variance_factor)
      variance = max(sum_squares/real(n, dp), tiny(1.0_dp))
      log_likelihood = -0.5_dp*(real(n, dp)*(log(two_pi*variance) + 1.0_dp) + &
                                  sum(log(variance_factor)))
   end subroutine exact_arma_likelihood

   pure subroutine exact_arma_likelihood_missing(x, ar, ma, mean_value, log_likelihood, &
                                                 variance, residuals, variance_factor, status)
      !! Evaluates exact stationary Gaussian likelihood using the observed covariance submatrix.
      real(dp), intent(in) :: x(:) !! Stationary observations containing quiet NaNs for gaps.
      real(dp), intent(in) :: ar(:) !! Effective autoregressive coefficients by lag.
      real(dp), intent(in) :: ma(:) !! Effective moving-average coefficients by lag.
      real(dp), intent(in) :: mean_value !! Stationary process mean.
      real(dp), intent(out) :: log_likelihood !! Profiled observed-data log likelihood.
      real(dp), intent(out) :: variance !! Profiled innovation variance.
      real(dp), allocatable, intent(out) :: residuals(:) !! Standardization-ready innovations.
      real(dp), allocatable, intent(out) :: variance_factor(:) !! Unit-scale innovation variances.
      integer, intent(out) :: status !! Zero on success or one for invalid covariance.
      real(dp), allocatable :: autocovariance(:), centered(:), covariance(:, :)
      real(dp), allocatable :: covariance_vector(:), solution(:), weights(:)
      integer, allocatable :: observed_indices(:)
      real(dp) :: determinant_sign, log_determinant, quiet_nan, quadratic
      integer :: i, info, j, k, n, observed_count

      n = size(x)
      observed_count = count(ieee_is_finite(x))
      quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate (residuals(n), variance_factor(n))
      residuals = quiet_nan
      variance_factor = quiet_nan
      log_likelihood = -huge(1.0_dp)
      variance = huge(1.0_dp)
      if (observed_count < 2) then
         status = 1
         return
      end if
      allocate (observed_indices(observed_count), centered(observed_count))
      k = 0
      do i = 1, n
         if (.not. ieee_is_finite(x(i))) cycle
         k = k + 1
         observed_indices(k) = i
         centered(k) = x(i) - mean_value
      end do
      call r_arma_autocovariance(ar, ma, n - 1, autocovariance, status)
      if (status /= r_ok) return
      allocate (covariance(observed_count, observed_count))
      do j = 1, observed_count
         do i = 1, observed_count
            covariance(i, j) = autocovariance(abs(observed_indices(i) - observed_indices(j)))
         end do
      end do
      allocate (solution(observed_count))
      call solve_system(covariance, centered, solution, info)
      if (info /= 0) then
         status = 1
         return
      end if
      call signed_log_determinant(covariance, determinant_sign, log_determinant, info)
      if (info /= 0 .or. determinant_sign <= 0.0_dp) then
         status = 1
         return
      end if
      quadratic = dot_product(centered, solution)
      variance = max(quadratic/real(observed_count, dp), tiny(1.0_dp))
      log_likelihood = -0.5_dp*(real(observed_count, dp)* &
                       (log(two_pi*variance) + 1.0_dp) + log_determinant)

      residuals(observed_indices(1)) = centered(1)
      variance_factor(observed_indices(1)) = autocovariance(0)
      allocate (covariance_vector(observed_count - 1), weights(observed_count - 1))
      do k = 2, observed_count
         do j = 1, k - 1
            covariance_vector(j) = &
               autocovariance(observed_indices(k) - observed_indices(j))
         end do
         call solve_system(covariance(:k - 1, :k - 1), covariance_vector(:k - 1), &
                           weights(:k - 1), info)
         if (info /= 0) then
            status = 1
            return
         end if
         residuals(observed_indices(k)) = centered(k) - &
            dot_product(weights(:k - 1), centered(:k - 1))
         variance_factor(observed_indices(k)) = autocovariance(0) - &
            dot_product(covariance_vector(:k - 1), weights(:k - 1))
      end do
      if (any(pack(variance_factor, ieee_is_finite(x)) <= 0.0_dp)) then
         status = 1
         return
      end if
      status = r_ok
   end subroutine exact_arma_likelihood_missing

   pure subroutine exact_diffuse_arma_likelihood(x, ar, ma, difference_order, &
                                                 seasonal_difference_order, period, &
                                                 log_likelihood, variance, residuals, &
                                                 variance_factor, status)
      !! Evaluates an integrated ARIMA likelihood with exact diffuse initialization.
      real(dp), intent(in) :: x(:) !! Regression-adjusted levels containing optional quiet NaNs.
      real(dp), intent(in) :: ar(:) !! Effective stationary autoregressive coefficients.
      real(dp), intent(in) :: ma(:) !! Effective stationary moving-average coefficients.
      integer, intent(in) :: difference_order !! Number of ordinary differences.
      integer, intent(in) :: seasonal_difference_order !! Number of seasonal differences.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(out) :: log_likelihood !! Profiled exact diffuse log likelihood.
      real(dp), intent(out) :: variance !! Profiled innovation variance.
      real(dp), allocatable, intent(out) :: residuals(:) !! One-step innovations on the level scale.
      real(dp), allocatable, intent(out) :: variance_factor(:) !! Unit-scale innovation variances.
      integer, intent(out) :: status !! Zero on success or one for an invalid state covariance.
      real(dp), allocatable :: difference_weights(:), finite_covariance(:, :), finite_new(:, :)
      real(dp), allocatable :: finite_stationary(:, :), finite_stationary_new(:, :)
      real(dp), allocatable :: diffuse_covariance(:, :), gain_diffuse(:), gain_finite(:)
      real(dp), allocatable :: state(:), state_noise(:), stationary_noise(:)
      real(dp), allocatable :: stationary_transition(:, :), transition(:, :)
      real(dp) :: diffuse_log_determinant, finite_log_determinant, finite_variance
      real(dp) :: diffuse_variance, innovation, quiet_nan, scale, sum_squares
      integer :: ar_state_count, diffuse_count, i, iteration, nstate, nuse
      integer :: p, q, stationary_state_count, t

      p = size(ar)
      q = size(ma)
      call difference_weights_for(difference_order, seasonal_difference_order, period, &
                                  difference_weights)
      diffuse_count = size(difference_weights)
      ar_state_count = max(1, p)
      stationary_state_count = ar_state_count + q
      allocate (stationary_transition(stationary_state_count, stationary_state_count), &
                stationary_noise(stationary_state_count), source=0.0_dp)
      if (p > 0) stationary_transition(1, :p) = ar
      do i = 2, ar_state_count
         stationary_transition(i, i - 1) = 1.0_dp
      end do
      if (q > 0) then
         stationary_transition(1, ar_state_count + 1:) = ma
         do i = 2, q
            stationary_transition(ar_state_count + i, ar_state_count + i - 1) = 1.0_dp
         end do
         stationary_noise(ar_state_count + 1) = 1.0_dp
      end if
      stationary_noise(1) = 1.0_dp
      allocate (finite_stationary(stationary_state_count, stationary_state_count), source=0.0_dp)
      allocate (finite_stationary_new(stationary_state_count, stationary_state_count))
      do iteration = 1, 20000
         finite_stationary_new = &
            matmul(stationary_transition, &
                   matmul(finite_stationary, transpose(stationary_transition))) + &
            spread(stationary_noise, 2, stationary_state_count)* &
            spread(stationary_noise, 1, stationary_state_count)
         scale = max(1.0_dp, maxval(abs(finite_stationary_new)))
         if (maxval(abs(finite_stationary_new - finite_stationary)) < &
             1.0e-13_dp*scale) exit
         finite_stationary = finite_stationary_new
      end do
      if (iteration > 20000) then
         status = 1
         return
      end if
      finite_stationary = finite_stationary_new

      nstate = diffuse_count + stationary_state_count
      allocate (transition(nstate, nstate), state_noise(nstate), state(nstate), &
                finite_covariance(nstate, nstate), diffuse_covariance(nstate, nstate), &
                source=0.0_dp)
      transition(1, :diffuse_count) = difference_weights
      transition(1, diffuse_count + 1:) = stationary_transition(1, :)
      state_noise(1) = stationary_noise(1)
      do i = 2, diffuse_count
         transition(i, i - 1) = 1.0_dp
      end do
      transition(diffuse_count + 1:, diffuse_count + 1:) = stationary_transition
      state_noise(diffuse_count + 1:) = stationary_noise
      do i = 1, diffuse_count
         diffuse_covariance(i, i) = 1.0_dp
      end do
      finite_covariance(diffuse_count + 1:, diffuse_count + 1:) = finite_stationary

      quiet_nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate (residuals(size(x)), variance_factor(size(x)))
      residuals = quiet_nan
      variance_factor = quiet_nan
      diffuse_log_determinant = 0.0_dp
      finite_log_determinant = 0.0_dp
      sum_squares = 0.0_dp
      nuse = 0
      do t = 1, size(x)
         diffuse_variance = max(0.0_dp, diffuse_covariance(1, 1))
         finite_variance = max(0.0_dp, finite_covariance(1, 1))
         if (ieee_is_finite(x(t))) then
            innovation = x(t) - state(1)
            residuals(t) = innovation
            if (diffuse_variance > 1.0e-10_dp) then
               diffuse_log_determinant = diffuse_log_determinant + log(diffuse_variance)
               gain_diffuse = diffuse_covariance(:, 1)/diffuse_variance
               gain_finite = (finite_covariance(:, 1) - &
                              gain_diffuse*finite_variance)/diffuse_variance
               state = state + gain_diffuse*innovation
               finite_covariance = finite_covariance - &
                  spread(gain_diffuse, 2, nstate)*spread(finite_covariance(1, :), 1, nstate) - &
                  spread(gain_finite, 2, nstate)*spread(diffuse_covariance(1, :), 1, nstate)
               diffuse_covariance = diffuse_covariance - &
                  spread(gain_diffuse, 2, nstate)*spread(diffuse_covariance(1, :), 1, nstate)
               finite_covariance = 0.5_dp*(finite_covariance + transpose(finite_covariance))
               diffuse_covariance = 0.5_dp*(diffuse_covariance + transpose(diffuse_covariance))
            else
               if (finite_variance <= tiny(1.0_dp)) then
                  status = 1
                  return
               end if
               variance_factor(t) = finite_variance
               gain_finite = finite_covariance(:, 1)/finite_variance
               state = state + gain_finite*innovation
               finite_covariance = finite_covariance - &
                  spread(gain_finite, 2, nstate)*spread(finite_covariance(1, :), 1, nstate)
               finite_covariance = 0.5_dp*(finite_covariance + transpose(finite_covariance))
               sum_squares = sum_squares + innovation**2/finite_variance
               finite_log_determinant = finite_log_determinant + log(finite_variance)
               nuse = nuse + 1
            end if
         end if
         if (t < size(x)) then
            state = matmul(transition, state)
            finite_new = matmul(transition, &
                                matmul(finite_covariance, transpose(transition))) + &
                         spread(state_noise, 2, nstate)*spread(state_noise, 1, nstate)
            finite_covariance = 0.5_dp*(finite_new + transpose(finite_new))
            diffuse_covariance = matmul(transition, &
                                        matmul(diffuse_covariance, transpose(transition)))
            diffuse_covariance = 0.5_dp*(diffuse_covariance + transpose(diffuse_covariance))
         end if
      end do
      if (maxval(abs(diffuse_covariance)) > 1.0e-7_dp .or. nuse < 1 .or. &
          sum_squares <= 0.0_dp) then
         status = 1
         return
      end if
      variance = max(sum_squares/real(nuse, dp), tiny(1.0_dp))
      log_likelihood = -0.5_dp*(real(nuse, dp)*(log(two_pi*variance) + 1.0_dp) + &
                       finite_log_determinant + diffuse_log_determinant)
      status = r_ok
   end subroutine exact_diffuse_arma_likelihood

   pure subroutine optimize_parameters(x, p, q, seasonal_p, seasonal_q, period, estimate_mean, &
                                       regressors, parameters, max_iterations, tolerance, objective_method, &
                                       difference_order, seasonal_difference_order, fixed_mask, &
                                       iterations, status)
      !! Minimizes a selected ARIMA objective using a deterministic Nelder-Mead simplex.
      real(dp), intent(in) :: x(:) !! Stationary or differenced series.
      integer, intent(in) :: p !! Autoregressive order.
      integer, intent(in) :: q !! Moving-average order.
      integer, intent(in) :: seasonal_p !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_q !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      logical, intent(in) :: estimate_mean !! Whether the final parameter is the series mean.
      real(dp), intent(in) :: regressors(:, :) !! Differenced user regressors.
      real(dp), intent(inout) :: parameters(:) !! Initial and optimized parameter vector.
      integer, intent(in) :: max_iterations !! Maximum simplex iterations.
      real(dp), intent(in) :: tolerance !! Parameter and objective convergence tolerance.
      integer, intent(in) :: objective_method !! One for CSS, two for stationary ML, or three for diffuse ML.
      integer, intent(in) :: difference_order !! Diffuse ordinary differencing order.
      integer, intent(in) :: seasonal_difference_order !! Diffuse seasonal differencing order.
      logical, intent(in) :: fixed_mask(:) !! True for each parameter excluded from optimization.
      integer, intent(out) :: iterations !! Completed simplex iterations.
      integer, intent(out) :: status !! Zero on convergence or two at the iteration limit.
      real(dp), allocatable :: centroid(:), contracted(:), expanded(:), reflected(:)
      real(dp), allocatable :: simplex(:, :), values(:)
      real(dp) :: contracted_value, expanded_value, reflected_value, step
      integer, allocatable :: free_indices(:)
      integer :: i, j, n, n_free

      n = size(parameters)
      n_free = count(.not. fixed_mask)
      allocate (free_indices(n_free))
      j = 0
      do i = 1, n
         if (fixed_mask(i)) cycle
         j = j + 1
         free_indices(j) = i
      end do
      allocate (simplex(n, n_free + 1), values(n_free + 1))
      allocate (centroid(n), contracted(n), expanded(n), reflected(n))
      simplex(:, 1) = parameters
      do i = 1, n_free
         simplex(:, i + 1) = parameters
         step = 0.1_dp
         j = free_indices(i)
         if (estimate_mean .and. j == n) step = max(0.1_dp, 0.1_dp*maxval(abs(x)))
         simplex(j, i + 1) = simplex(j, i + 1) + step
      end do
      do i = 1, n_free + 1
         values(i) = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                     estimate_mean, regressors, simplex(:, i), objective_method, &
                                     difference_order, seasonal_difference_order)
      end do

      status = 2
      do iterations = 1, max_iterations
         call sort_simplex(simplex, values)
         if (maxval(abs(simplex(:, 2:) - spread(simplex(:, 1), 2, n_free))) <= tolerance .and. &
             maxval(abs(values(2:) - values(1))) <= tolerance) then
            status = 0
            exit
         end if
         centroid = sum(simplex(:, :n_free), dim=2)/real(n_free, dp)
         reflected = 2.0_dp*centroid - simplex(:, n_free + 1)
         where (fixed_mask) reflected = parameters
         reflected_value = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                           estimate_mean, regressors, reflected, objective_method, &
                                           difference_order, seasonal_difference_order)
         if (reflected_value < values(1)) then
            expanded = centroid + 2.0_dp*(reflected - centroid)
            where (fixed_mask) expanded = parameters
            expanded_value = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                             estimate_mean, regressors, expanded, objective_method, &
                                             difference_order, seasonal_difference_order)
            if (expanded_value < reflected_value) then
               simplex(:, n_free + 1) = expanded
               values(n_free + 1) = expanded_value
            else
               simplex(:, n_free + 1) = reflected
               values(n_free + 1) = reflected_value
            end if
         else if (reflected_value < values(n_free)) then
            simplex(:, n_free + 1) = reflected
            values(n_free + 1) = reflected_value
         else
            if (reflected_value < values(n_free + 1)) then
               contracted = centroid + 0.5_dp*(reflected - centroid)
            else
               contracted = centroid + 0.5_dp*(simplex(:, n_free + 1) - centroid)
            end if
            where (fixed_mask) contracted = parameters
            contracted_value = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                               estimate_mean, regressors, contracted, objective_method, &
                                               difference_order, seasonal_difference_order)
            if (contracted_value < min(reflected_value, values(n_free + 1))) then
               simplex(:, n_free + 1) = contracted
               values(n_free + 1) = contracted_value
            else
               do i = 2, n_free + 1
                  simplex(:, i) = simplex(:, 1) + 0.5_dp*(simplex(:, i) - simplex(:, 1))
                  where (fixed_mask) simplex(:, i) = parameters
                  values(i) = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                              estimate_mean, regressors, simplex(:, i), objective_method, &
                                              difference_order, seasonal_difference_order)
               end do
            end if
         end if
      end do
      if (status /= 0) iterations = max_iterations
      call sort_simplex(simplex, values)
      parameters = simplex(:, 1)
   end subroutine optimize_parameters

   pure subroutine arima_parameter_covariance(x, p, q, seasonal_p, seasonal_q, period, &
                                              estimate_mean, regressors, parameters, fixed_mask, &
                                              objective_method, difference_order, &
                                              seasonal_difference_order, covariance, status)
      !! Approximates the observed free-parameter covariance using R's Hessian convention.
      real(dp), intent(in) :: x(:) !! Stationary or differenced observations.
      integer, intent(in) :: p !! Ordinary autoregressive order.
      integer, intent(in) :: q !! Ordinary moving-average order.
      integer, intent(in) :: seasonal_p !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_q !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      logical, intent(in) :: estimate_mean !! Whether the final parameter is a mean.
      real(dp), intent(in) :: regressors(:, :) !! Differenced user regressors.
      real(dp), intent(in) :: parameters(:) !! Fitted full parameter vector.
      logical, intent(in) :: fixed_mask(:) !! True for parameters excluded from the Hessian.
      integer, intent(in) :: objective_method !! One for CSS, two for stationary ML, or three for diffuse ML.
      integer, intent(in) :: difference_order !! Diffuse ordinary differencing order.
      integer, intent(in) :: seasonal_difference_order !! Diffuse seasonal differencing order.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Full covariance with fixed rows zero.
      integer, intent(out) :: status !! Zero on success or one when the Hessian is singular.
      real(dp), parameter :: finite_difference_step = 1.0e-3_dp
      real(dp), allocatable :: free_covariance(:, :), hessian(:, :), step(:), trial(:)
      integer, allocatable :: free_indices(:)
      real(dp) :: f_mm, f_mp, f_pm, f_pp
      integer :: a, b, i, info, j, n_free

      allocate (covariance(size(parameters), size(parameters)), source=0.0_dp)
      n_free = count(.not. fixed_mask)
      if (n_free == 0) then
         status = 0
         return
      end if
      allocate (free_indices(n_free), hessian(n_free, n_free), step(n_free))
      allocate (trial(size(parameters)))
      j = 0
      do i = 1, size(parameters)
         if (fixed_mask(i)) cycle
         j = j + 1
         free_indices(j) = i
         step(j) = finite_difference_step
      end do

      do b = 1, n_free
         j = free_indices(b)
         do a = 1, b
            i = free_indices(a)
            trial = parameters
            trial(i) = trial(i) + step(a)
            trial(j) = trial(j) + step(b)
            f_pp = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                   estimate_mean, regressors, trial, objective_method, &
                                   difference_order, seasonal_difference_order)
            trial = parameters
            trial(i) = trial(i) + step(a)
            trial(j) = trial(j) - step(b)
            f_pm = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                   estimate_mean, regressors, trial, objective_method, &
                                   difference_order, seasonal_difference_order)
            trial = parameters
            trial(i) = trial(i) - step(a)
            trial(j) = trial(j) + step(b)
            f_mp = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                   estimate_mean, regressors, trial, objective_method, &
                                   difference_order, seasonal_difference_order)
            trial = parameters
            trial(i) = trial(i) - step(a)
            trial(j) = trial(j) - step(b)
            f_mm = arima_objective(x, p, q, seasonal_p, seasonal_q, period, &
                                   estimate_mean, regressors, trial, objective_method, &
                                   difference_order, seasonal_difference_order)
            hessian(a, b) = (f_pp - f_pm - f_mp + f_mm)/(4.0_dp*step(a)*step(b))
            hessian(b, a) = hessian(a, b)
         end do
      end do
      if (.not. all(ieee_is_finite(hessian))) then
         status = 1
         return
      end if
      call inverse_matrix(real(count(ieee_is_finite(x)) - difference_order - &
                          seasonal_difference_order*period, dp)*hessian, &
                          free_covariance, info)
      if (info /= 0 .or. .not. all(ieee_is_finite(free_covariance))) then
         status = 1
         return
      end if
      do b = 1, n_free
         do a = 1, n_free
            covariance(free_indices(a), free_indices(b)) = free_covariance(a, b)
         end do
      end do
      status = 0
   end subroutine arima_parameter_covariance

   pure subroutine sort_simplex(simplex, values)
      !! Sorts simplex columns in ascending objective order.
      real(dp), intent(inout) :: simplex(:, :) !! Parameter vertices arranged by column.
      real(dp), intent(inout) :: values(:) !! Objective value for each vertex.
      real(dp), allocatable :: temporary(:)
      real(dp) :: temporary_value
      integer :: i, j

      allocate (temporary(size(simplex, 1)))
      do i = 2, size(values)
         j = i
         do while (j > 1)
            if (values(j) >= values(j - 1)) exit
            temporary_value = values(j)
            values(j) = values(j - 1)
            values(j - 1) = temporary_value
            temporary = simplex(:, j)
            simplex(:, j) = simplex(:, j - 1)
            simplex(:, j - 1) = temporary
            j = j - 1
         end do
      end do
   end subroutine sort_simplex

   pure subroutine expanded_coefficients(p, q, seasonal_p, seasonal_q, period, &
                                         parameters, full_ar, full_ma)
      !! Expands multiplicative ordinary and seasonal AR and MA polynomials.
      integer, intent(in) :: p !! Ordinary autoregressive order.
      integer, intent(in) :: q !! Ordinary moving-average order.
      integer, intent(in) :: seasonal_p !! Seasonal autoregressive order.
      integer, intent(in) :: seasonal_q !! Seasonal moving-average order.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), intent(in) :: parameters(:) !! Vector beginning with `[ar, ma, sar, sma]`.
      real(dp), allocatable, intent(out) :: full_ar(:) !! Expanded AR coefficients by lag.
      real(dp), allocatable, intent(out) :: full_ma(:) !! Expanded MA coefficients by lag.
      real(dp), allocatable :: ordinary(:), polynomial(:), seasonal(:)
      integer :: i

      allocate (ordinary(0:p), source=0.0_dp)
      allocate (seasonal(0:seasonal_p*period), source=0.0_dp)
      ordinary(0) = 1.0_dp
      seasonal(0) = 1.0_dp
      do i = 1, p
         ordinary(i) = -parameters(i)
      end do
      do i = 1, seasonal_p
         seasonal(i*period) = -parameters(p + q + i)
      end do
      call multiply_polynomials(ordinary, seasonal, polynomial)
      allocate (full_ar(ubound(polynomial, 1)))
      if (size(full_ar) > 0) full_ar = -polynomial(1:)

      deallocate (ordinary, polynomial, seasonal)
      allocate (ordinary(0:q), source=0.0_dp)
      allocate (seasonal(0:seasonal_q*period), source=0.0_dp)
      ordinary(0) = 1.0_dp
      seasonal(0) = 1.0_dp
      do i = 1, q
         ordinary(i) = parameters(p + i)
      end do
      do i = 1, seasonal_q
         seasonal(i*period) = parameters(p + q + seasonal_p + i)
      end do
      call multiply_polynomials(ordinary, seasonal, polynomial)
      allocate (full_ma(ubound(polynomial, 1)))
      if (size(full_ma) > 0) full_ma = polynomial(1:)
   end subroutine expanded_coefficients

   pure subroutine multiply_polynomials(left, right, product)
      !! Multiplies two lag polynomials whose lower bounds represent lag zero.
      real(dp), intent(in) :: left(0:) !! First polynomial indexed by lag.
      real(dp), intent(in) :: right(0:) !! Second polynomial indexed by lag.
      real(dp), allocatable, intent(out) :: product(:) !! Product polynomial indexed by lag.
      integer :: i, j

      allocate (product(0:ubound(left, 1) + ubound(right, 1)), source=0.0_dp)
      do j = 0, ubound(right, 1)
         do i = 0, ubound(left, 1)
            product(i + j) = product(i + j) + left(i)*right(j)
         end do
      end do
   end subroutine multiply_polynomials

   pure subroutine difference_weights_for(difference_order, seasonal_difference_order, &
                                           period, weights)
      !! Builds coefficients that reconstruct observations from a differenced series.
      integer, intent(in) :: difference_order !! Number of ordinary differences.
      integer, intent(in) :: seasonal_difference_order !! Number of seasonal differences.
      integer, intent(in) :: period !! Seasonal period.
      real(dp), allocatable, intent(out) :: weights(:) !! Prior-observation reconstruction weights.
      real(dp), allocatable :: polynomial(:)
      integer :: i

      allocate (polynomial(0:0), source=1.0_dp)
      do i = 1, difference_order
         call multiply_polynomial_factor(polynomial, 1, -1.0_dp)
      end do
      do i = 1, seasonal_difference_order
         call multiply_polynomial_factor(polynomial, period, -1.0_dp)
      end do
      allocate (weights(ubound(polynomial, 1)))
      if (size(weights) > 0) weights = -polynomial(1:)
   end subroutine difference_weights_for

   pure subroutine multiply_polynomial_factor(polynomial, lag, coefficient)
      !! Multiplies a lag polynomial by `(1 + coefficient*B**lag)`.
      real(dp), allocatable, intent(inout) :: polynomial(:) !! Polynomial indexed from lag zero.
      integer, intent(in) :: lag !! Positive lag of the second factor term.
      real(dp), intent(in) :: coefficient !! Coefficient of the lagged factor term.
      real(dp), allocatable :: product(:)
      integer :: i, old_order

      old_order = ubound(polynomial, 1)
      allocate (product(0:old_order + lag), source=0.0_dp)
      do i = 0, old_order
         product(i) = product(i) + polynomial(i)
         product(i + lag) = product(i + lag) + coefficient*polynomial(i)
      end do
      call move_alloc(product, polynomial)
   end subroutine multiply_polynomial_factor

end module r_stats_arima
