! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Workflow algorithms translated from the R rugarch package.
module rugarch_extensions_mod
   !! Bootstrap, rolling, multi-model, parameter, and diagnostic workflows.
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_pseudoinverse
   use random_mod, only: set_random_seed, random_standard_normal
   use rugarch_mod, only: rugarch_spec_t, rugarch_parameters_t, rugarch_fit_t, &
      rugarch_filter_t, rugarch_forecast_t, rugarch_fit, rugarch_filter, &
      rugarch_forecast, rugarch_coefficients, &
      rugarch_model_sgarch, rugarch_model_igarch, rugarch_model_egarch, &
      rugarch_model_gjrgarch, rugarch_model_aparch, rugarch_model_figarch, &
      rugarch_model_csgarch, rugarch_model_realgarch, rugarch_model_fgarch, &
      random_standardized_innovation, rugarch_fgarch_exponent
   use stats_mod, only: normal_quantile
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   type, public :: rugarch_bootstrap_t
      !! Parametric predictive paths and pointwise forecast intervals.
      real(dp), allocatable :: paths(:, :)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: level = 0.95_dp
      integer :: info = 0
   end type rugarch_bootstrap_t

   type, public :: rugarch_roll_t
      !! One-step rolling forecasts with periodic model refitting.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: actual(:)
      real(dp), allocatable :: standardized_error(:)
      integer, allocatable :: origin(:)
      logical, allocatable :: refitted(:)
      integer :: info = 0
   end type rugarch_roll_t

   type, public :: rugarch_multifit_t
      !! Fits of multiple specifications and information-criterion selection.
      type(rugarch_fit_t), allocatable :: fits(:)
      integer :: best_aic = 0
      integer :: best_bic = 0
      integer :: info = 0
   end type rugarch_multifit_t

   type, public :: rugarch_multifilter_t
      !! Fixed-parameter filters for several volatility specifications.
      type(rugarch_filter_t), allocatable :: filters(:)
      integer :: info = 0
   end type rugarch_multifilter_t

   type, public :: rugarch_multiforecast_t
      !! Forecasts from every successful fit in a multi-model collection.
      type(rugarch_forecast_t), allocatable :: forecasts(:)
      integer :: info = 0
   end type rugarch_multiforecast_t

   type, public :: rugarch_parameter_distribution_t
      !! Gaussian approximation to the fitted coefficient distribution.
      real(dp), allocatable :: draws(:, :)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: standard_deviation(:)
      integer :: info = 0
   end type rugarch_parameter_distribution_t

   type, public :: rugarch_model_confidence_t
      !! Loss-based model-confidence summary using paired normal tests.
      real(dp), allocatable :: mean_loss(:)
      real(dp), allocatable :: p_value(:)
      logical, allocatable :: included(:)
      integer :: best_model = 0
      real(dp) :: level = 0.90_dp
      integer :: info = 0
   end type rugarch_model_confidence_t

   type, public :: rugarch_news_impact_t
      !! Conditional-variance response over a grid of standardized shocks.
      real(dp), allocatable :: shock(:)
      real(dp), allocatable :: variance(:)
      real(dp) :: previous_variance = 0.0_dp
      integer :: info = 0
   end type rugarch_news_impact_t

   type, public :: rugarch_sign_bias_test_t
      !! Engle-Ng sign and size-bias auxiliary-regression tests.
      real(dp) :: sign_statistic = 0.0_dp
      real(dp) :: negative_size_statistic = 0.0_dp
      real(dp) :: positive_size_statistic = 0.0_dp
      real(dp) :: joint_statistic = 0.0_dp
      real(dp) :: sign_p_value = 1.0_dp
      real(dp) :: negative_size_p_value = 1.0_dp
      real(dp) :: positive_size_p_value = 1.0_dp
      real(dp) :: joint_p_value = 1.0_dp
      integer :: info = 0
   end type rugarch_sign_bias_test_t

   interface display
      module procedure display_rugarch_bootstrap
      module procedure display_rugarch_roll
      module procedure display_rugarch_multifit
      module procedure display_rugarch_multifilter
      module procedure display_rugarch_multiforecast
      module procedure display_rugarch_parameter_distribution
      module procedure display_rugarch_model_confidence
      module procedure display_rugarch_news_impact
      module procedure display_rugarch_sign_bias_test
   end interface display

   public :: rugarch_bootstrap_forecast, rugarch_roll
   public :: rugarch_multifit, rugarch_parameter_distribution
   public :: rugarch_multifilter, rugarch_multiforecast
   public :: rugarch_model_confidence, rugarch_news_impact
   public :: rugarch_sign_bias_test, display

contains

   function rugarch_bootstrap_forecast(fit, horizon, draws, level, seed) &
      result(out)
      !! Draw conditional predictive paths from the fitted location and scale.
      type(rugarch_fit_t), intent(in) :: fit !! Fitted volatility model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      integer, intent(in) :: draws !! Number of predictive paths.
      real(dp), intent(in), optional :: level !! Central interval coverage.
      integer, intent(in), optional :: seed !! Optional random-number seed.
      type(rugarch_bootstrap_t) :: out
      type(rugarch_forecast_t) :: forecast
      type(rugarch_filter_t) :: filtered
      real(dp), allocatable :: ordered(:), history(:)
      real(dp) :: probability
      integer :: path, step, lower_index, upper_index, observations

      if (horizon < 1 .or. draws < 2 .or. fit%info /= 0) then
         out%info = 1
         return
      end if
      if (present(level)) out%level = level
      if (out%level <= 0.0_dp .or. out%level >= 1.0_dp) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      forecast = rugarch_forecast(fit, horizon)
      if (forecast%info /= 0) then
         out%info = 2
         return
      end if
      allocate(out%paths(horizon, draws), out%mean(horizon), &
         out%lower(horizon), out%upper(horizon), ordered(draws))
      observations = size(fit%filtered%residuals)
      allocate(history(observations + horizon))
      do path = 1, draws
         history(:observations) = fit%filtered%conditional_mean + &
            fit%filtered%residuals
         do step = 1, horizon
            if (fit%specification%variance_model /= rugarch_model_realgarch &
               .and. fit%specification%mean_regressor_count == 0 .and. &
               fit%specification%variance_regressor_count == 0) then
               history(observations + step) = 0.0_dp
               filtered = rugarch_filter(history(:observations + step), &
                  fit%specification, fit%parameters)
            else
               filtered%info = 1
            end if
            if (filtered%info == 0) then
               out%paths(step, path) = filtered%conditional_mean( &
                  observations + step) + filtered%conditional_sigma( &
                  observations + step)*random_standardized_innovation( &
                  fit%specification%distribution, fit%parameters%shape, &
                  fit%parameters%skew, fit%parameters%lambda)
            else
               out%paths(step, path) = forecast%mean(step) + &
                  forecast%sigma(step)*random_standardized_innovation( &
               fit%specification%distribution, fit%parameters%shape, &
               fit%parameters%skew, fit%parameters%lambda)
            end if
            history(observations + step) = out%paths(step, path)
         end do
      end do
      probability = 0.5_dp*(1.0_dp - out%level)
      lower_index = max(1, ceiling(probability*real(draws, dp)))
      upper_index = min(draws, ceiling((1.0_dp - probability)*real(draws, dp)))
      do step = 1, horizon
         ordered = out%paths(step, :)
         call insertion_sort(ordered)
         out%mean(step) = sum(ordered)/real(draws, dp)
         out%lower(step) = ordered(lower_index)
         out%upper(step) = ordered(upper_index)
      end do
   end function rugarch_bootstrap_forecast

   pure function rugarch_roll(series, specification, initial_window, &
      refit_every, moving_window, max_iterations) result(out)
      !! Produce one-step rolling forecasts with expanding or moving samples.
      real(dp), intent(in) :: series(:) !! Complete observed series.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      integer, intent(in) :: initial_window !! Observations in the first fit.
      integer, intent(in), optional :: refit_every !! Origins between refits.
      logical, intent(in), optional :: moving_window !! Use a fixed-length window.
      integer, intent(in), optional :: max_iterations !! Optimizer iteration limit.
      type(rugarch_roll_t) :: out
      type(rugarch_fit_t) :: fit
      type(rugarch_forecast_t) :: forecast
      integer :: count, every, origin, first, index, iterations
      logical :: moving

      count = size(series) - initial_window
      if (initial_window < 5 .or. count < 1) then
         out%info = 1
         return
      end if
      every = 1
      if (present(refit_every)) every = refit_every
      moving = .false.
      if (present(moving_window)) moving = moving_window
      iterations = 200
      if (present(max_iterations)) iterations = max_iterations
      if (every < 1 .or. iterations < 1) then
         out%info = 1
         return
      end if
      allocate(out%mean(count), out%sigma(count), out%actual(count), &
         out%standardized_error(count), out%origin(count), &
         out%refitted(count))
      do index = 1, count
         origin = initial_window + index - 1
         out%origin(index) = origin
         out%actual(index) = series(origin + 1)
         out%refitted(index) = index == 1 .or. mod(index - 1, every) == 0
         if (out%refitted(index)) then
            first = 1
            if (moving) first = origin - initial_window + 1
            fit = rugarch_fit(series(first:origin), specification, &
               max_iterations=iterations)
         end if
         if (fit%info /= 0) then
            out%info = 2
            return
         end if
         forecast = rugarch_forecast(fit, 1)
         if (forecast%info /= 0) then
            out%info = 2
            return
         end if
         out%mean(index) = forecast%mean(1)
         out%sigma(index) = forecast%sigma(1)
         out%standardized_error(index) = &
            (out%actual(index) - out%mean(index))/max(out%sigma(index), &
            tiny(1.0_dp))
      end do
   end function rugarch_roll

   pure function rugarch_multifit(series, specifications, max_iterations) &
      result(out)
      !! Fit several specifications and identify minimum-AIC and minimum-BIC fits.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(rugarch_spec_t), intent(in) :: specifications(:) !! Candidate models.
      integer, intent(in), optional :: max_iterations !! Optimizer iteration limit.
      type(rugarch_multifit_t) :: out
      integer :: model, iterations

      if (size(specifications) < 1) then
         out%info = 1
         return
      end if
      iterations = 200
      if (present(max_iterations)) iterations = max_iterations
      allocate(out%fits(size(specifications)))
      do model = 1, size(specifications)
         out%fits(model) = rugarch_fit(series, specifications(model), &
            max_iterations=iterations)
      end do
      out%best_aic = best_fit_index(out%fits, .true.)
      out%best_bic = best_fit_index(out%fits, .false.)
      if (out%best_aic == 0) out%info = 2
   end function rugarch_multifit

   pure function rugarch_multifilter(series, specifications, parameters) &
      result(out)
      !! Filter one series under several specifications and parameter vectors.
      real(dp), intent(in) :: series(:) !! Observed time series.
      type(rugarch_spec_t), intent(in) :: specifications(:) !! Model specifications.
      type(rugarch_parameters_t), intent(in) :: parameters(:) !! Physical parameters.
      type(rugarch_multifilter_t) :: out
      integer :: model

      if (size(specifications) < 1 .or. &
         size(parameters) /= size(specifications)) then
         out%info = 1
         return
      end if
      allocate(out%filters(size(specifications)))
      do model = 1, size(specifications)
         out%filters(model) = rugarch_filter(series, specifications(model), &
            parameters(model))
         if (out%filters(model)%info /= 0) out%info = 2
      end do
   end function rugarch_multifilter

   pure function rugarch_multiforecast(multifit, horizon) result(out)
      !! Forecast every successful member of a fitted model collection.
      type(rugarch_multifit_t), intent(in) :: multifit !! Multi-model fit result.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      type(rugarch_multiforecast_t) :: out
      integer :: model

      if (horizon < 1 .or. .not. allocated(multifit%fits) .or. &
         size(multifit%fits) < 1) then
         out%info = 1
         return
      end if
      allocate(out%forecasts(size(multifit%fits)))
      do model = 1, size(multifit%fits)
         out%forecasts(model) = rugarch_forecast(multifit%fits(model), horizon)
         if (out%forecasts(model)%info /= 0) out%info = 2
      end do
   end function rugarch_multiforecast

   function rugarch_parameter_distribution(fit, draws, seed) result(out)
      !! Simulate coefficients from the fitted asymptotic Gaussian distribution.
      type(rugarch_fit_t), intent(in) :: fit !! Fitted model and covariance.
      integer, intent(in) :: draws !! Number of coefficient draws.
      integer, intent(in), optional :: seed !! Optional random-number seed.
      type(rugarch_parameter_distribution_t) :: out
      real(dp), allocatable :: root(:, :), normal(:)
      integer :: draw, parameter, count

      count = size(fit%coefficients)
      if (fit%info /= 0 .or. draws < 1 .or. count < 1 .or. &
         size(fit%covariance, 1) /= count) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      root = covariance_root(fit%covariance)
      allocate(out%draws(count, draws), out%mean(count), &
         out%standard_deviation(count), normal(count))
      do draw = 1, draws
         do parameter = 1, count
            normal(parameter) = random_standard_normal()
         end do
         out%draws(:, draw) = fit%coefficients + matmul(root, normal)
      end do
      out%mean = sum(out%draws, dim=2)/real(draws, dp)
      if (draws > 1) then
         do parameter = 1, count
            out%standard_deviation(parameter) = sqrt(sum((out%draws(parameter, :) &
               - out%mean(parameter))**2)/real(draws - 1, dp))
         end do
      else
         out%standard_deviation = 0.0_dp
      end if
   end function rugarch_parameter_distribution

   pure function rugarch_model_confidence(loss, level) result(out)
      !! Compare model losses to the best model using paired normal statistics.
      real(dp), intent(in) :: loss(:, :) !! Losses by observation and model.
      real(dp), intent(in), optional :: level !! Confidence-set coverage.
      type(rugarch_model_confidence_t) :: out
      real(dp), allocatable :: difference(:)
      real(dp) :: average, standard_error, critical
      integer :: model, observations, models

      observations = size(loss, 1)
      models = size(loss, 2)
      if (observations < 2 .or. models < 1) then
         out%info = 1
         return
      end if
      if (present(level)) out%level = level
      if (out%level <= 0.0_dp .or. out%level >= 1.0_dp) then
         out%info = 1
         return
      end if
      allocate(out%mean_loss(models), out%p_value(models), &
         out%included(models), difference(observations))
      out%mean_loss = sum(loss, dim=1)/real(observations, dp)
      out%best_model = minloc(out%mean_loss, dim=1)
      critical = normal_quantile(0.5_dp*(1.0_dp + out%level))
      do model = 1, models
         difference = loss(:, model) - loss(:, out%best_model)
         average = sum(difference)/real(observations, dp)
         standard_error = sqrt(sum((difference - average)**2)/ &
            real(observations - 1, dp)/real(observations, dp))
         if (standard_error > tiny(1.0_dp)) then
            out%p_value(model) = 2.0_dp*normal_upper_tail(abs(average/standard_error))
            out%included(model) = average/standard_error <= critical
         else
            out%p_value(model) = merge(1.0_dp, 0.0_dp, abs(average) < tiny(1.0_dp))
            out%included(model) = average <= 0.0_dp
         end if
      end do
      out%included(out%best_model) = .true.
   end function rugarch_model_confidence

   pure function rugarch_news_impact(specification, parameters, shocks, &
      previous_variance) result(out)
      !! Evaluate the one-step variance response to standardized return shocks.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp), intent(in) :: shocks(:) !! Standardized shock grid.
      real(dp), intent(in), optional :: previous_variance !! Conditioning variance.
      type(rugarch_news_impact_t) :: out
      real(dp) :: h, residual, scale
      integer :: point

      h = max(parameters%omega/max(1.0_dp - sum(parameters%alpha) - &
         sum(parameters%beta), 0.05_dp), tiny(1.0_dp))
      if (present(previous_variance)) h = previous_variance
      if (h <= 0.0_dp .or. size(shocks) < 1 .or. &
         size(parameters%alpha) < 1) then
         out%info = 1
         return
      end if
      out%previous_variance = h
      allocate(out%shock(size(shocks)), out%variance(size(shocks)))
      out%shock = shocks
      do point = 1, size(shocks)
         residual = sqrt(h)*shocks(point)
         select case (specification%variance_model)
         case (rugarch_model_egarch)
            scale = parameters%omega + parameters%alpha(1)*abs(shocks(point)) &
               + parameters%asymmetry(1)*shocks(point) + &
               sum(parameters%beta)*log(h)
            out%variance(point) = exp(scale)
         case (rugarch_model_gjrgarch)
            out%variance(point) = parameters%omega + &
               parameters%alpha(1)*residual**2 + sum(parameters%beta)*h
            if (residual < 0.0_dp) out%variance(point) = out%variance(point) + &
               parameters%asymmetry(1)*residual**2
         case (rugarch_model_aparch)
            scale = parameters%omega + parameters%alpha(1)* &
               (abs(residual) - parameters%asymmetry(1)*residual)** &
               parameters%power + sum(parameters%beta)* &
               h**(0.5_dp*parameters%power)
            out%variance(point) = scale**(2.0_dp/parameters%power)
         case (rugarch_model_fgarch)
            scale = parameters%omega + parameters%alpha(1)* &
               (sqrt(1.0e-6_dp + (shocks(point) - &
               parameters%fgarch_shift(1))**2) - parameters%asymmetry(1)* &
               (shocks(point) - parameters%fgarch_shift(1)))** &
               rugarch_fgarch_exponent(specification, parameters)* &
               sqrt(h)**parameters%fgarch_lambda + &
               sum(parameters%beta)*sqrt(h)**parameters%fgarch_lambda
            out%variance(point) = max(scale, tiny(1.0_dp))** &
               (2.0_dp/parameters%fgarch_lambda)
         case default
            out%variance(point) = parameters%omega + &
               parameters%alpha(1)*residual**2 + sum(parameters%beta)*h
         end select
         out%variance(point) = max(out%variance(point), tiny(1.0_dp))
      end do
   end function rugarch_news_impact

   pure function rugarch_sign_bias_test(standardized_residuals) result(out)
      !! Apply the Engle-Ng sign, negative-size, positive-size, and joint tests.
      real(dp), intent(in) :: standardized_residuals(:) !! Standardized innovations.
      type(rugarch_sign_bias_test_t) :: out
      real(dp), allocatable :: x(:, :), y(:), coefficients(:), covariance(:, :)
      real(dp), allocatable :: gram(:, :), inverse(:, :), residual(:)
      real(dp) :: sigma2
      integer :: observations, time, inverse_info

      observations = size(standardized_residuals) - 1
      if (observations < 8) then
         out%info = 1
         return
      end if
      allocate(x(observations, 4), y(observations), gram(4, 4), inverse(4, 4))
      x(:, 1) = 1.0_dp
      do time = 1, observations
         y(time) = standardized_residuals(time + 1)**2
         x(time, 2) = merge(1.0_dp, 0.0_dp, standardized_residuals(time) < 0.0_dp)
         x(time, 3) = x(time, 2)*standardized_residuals(time)
         x(time, 4) = (1.0_dp - x(time, 2))*standardized_residuals(time)
      end do
      gram = matmul(transpose(x), x)
      call symmetric_pseudoinverse(gram, inverse, inverse_info)
      if (inverse_info /= 0) then
         out%info = 2
         return
      end if
      coefficients = matmul(inverse, matmul(transpose(x), y))
      residual = y - matmul(x, coefficients)
      sigma2 = sum(residual**2)/real(max(1, observations - 4), dp)
      covariance = sigma2*inverse
      out%sign_statistic = safe_t(coefficients(2), covariance(2, 2))
      out%negative_size_statistic = safe_t(coefficients(3), covariance(3, 3))
      out%positive_size_statistic = safe_t(coefficients(4), covariance(4, 4))
      out%sign_p_value = 2.0_dp*normal_upper_tail(abs(out%sign_statistic))
      out%negative_size_p_value = 2.0_dp* &
         normal_upper_tail(abs(out%negative_size_statistic))
      out%positive_size_p_value = 2.0_dp* &
         normal_upper_tail(abs(out%positive_size_statistic))
      out%joint_statistic = max(0.0_dp, sum(y**2) - sum(residual**2))/ &
         max(sum(residual**2)/real(observations, dp), tiny(1.0_dp))
      out%joint_p_value = exp(-0.5_dp*out%joint_statistic)* &
         (1.0_dp + sqrt(2.0_dp*out%joint_statistic/acos(-1.0_dp)))
      out%joint_p_value = min(1.0_dp, out%joint_p_value)
   end function rugarch_sign_bias_test

   pure integer function best_fit_index(fits, use_aic) result(index)
      !! Return the successful fit with the smallest requested criterion.
      type(rugarch_fit_t), intent(in) :: fits(:) !! Candidate fits.
      logical, intent(in) :: use_aic !! Select AIC instead of BIC.
      real(dp) :: best, criterion
      integer :: model

      index = 0
      best = huge(1.0_dp)
      do model = 1, size(fits)
         if (fits(model)%info /= 0) cycle
         criterion = merge(fits(model)%aic, fits(model)%bic, use_aic)
         if (criterion < best) then
            best = criterion
            index = model
         end if
      end do
   end function best_fit_index

   pure function covariance_root(covariance) result(root)
      !! Compute a clipped lower-triangular square root of a covariance matrix.
      real(dp), intent(in) :: covariance(:, :) !! Symmetric covariance matrix.
      real(dp) :: root(size(covariance, 1), size(covariance, 2))
      real(dp) :: value
      integer :: row, column, inner, count

      count = size(covariance, 1)
      root = 0.0_dp
      do row = 1, count
         do column = 1, row
            value = 0.5_dp*(covariance(row, column) + covariance(column, row))
            do inner = 1, column - 1
               value = value - root(row, inner)*root(column, inner)
            end do
            if (row == column) then
               root(row, column) = sqrt(max(value, 0.0_dp))
            else if (root(column, column) > tiny(1.0_dp)) then
               root(row, column) = value/root(column, column)
            end if
         end do
      end do
   end function covariance_root

   pure elemental real(dp) function safe_t(coefficient, variance_value) &
      result(statistic)
      !! Form a t statistic while guarding a zero estimated variance.
      real(dp), intent(in) :: coefficient !! Estimated coefficient.
      real(dp), intent(in) :: variance_value !! Estimated coefficient variance.

      statistic = coefficient/sqrt(max(variance_value, tiny(1.0_dp)))
   end function safe_t

   pure elemental real(dp) function normal_upper_tail(value) result(probability)
      !! Return the standard-normal upper-tail probability.
      real(dp), intent(in) :: value !! Standard-normal quantile.

      probability = 0.5_dp*erfc(value/sqrt(2.0_dp))
   end function normal_upper_tail

   pure subroutine insertion_sort(values)
      !! Sort a short real vector in ascending order in place.
      real(dp), intent(inout) :: values(:) !! Values to sort.
      real(dp) :: current
      integer :: index, previous

      do index = 2, size(values)
         current = values(index)
         previous = index - 1
         do while (previous >= 1)
            if (values(previous) <= current) exit
            values(previous + 1) = values(previous)
            previous = previous - 1
         end do
         values(previous + 1) = current
      end do
   end subroutine insertion_sort

   subroutine display_rugarch_bootstrap(value, unit)
      !! Display a bootstrap forecast summary.
      type(rugarch_bootstrap_t), intent(in) :: value !! Bootstrap result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch bootstrap forecast'
      write(output, '(a,i0)') 'paths: ', size(value%paths, 2)
      write(output, '(a,i0)') 'horizon: ', size(value%paths, 1)
      write(output, '(a,f8.4)') 'level: ', value%level
      write(output, '(a,i0)') 'info: ', value%info
   end subroutine display_rugarch_bootstrap

   subroutine display_rugarch_roll(value, unit)
      !! Display a rolling-forecast summary.
      type(rugarch_roll_t), intent(in) :: value !! Rolling result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch rolling forecast'
      write(output, '(a,i0)') 'forecasts: ', size(value%mean)
      write(output, '(a,i0)') 'refits: ', count(value%refitted)
      write(output, '(a,i0)') 'info: ', value%info
   end subroutine display_rugarch_roll

   subroutine display_rugarch_multifit(value, unit)
      !! Display a multi-model fit summary.
      type(rugarch_multifit_t), intent(in) :: value !! Multi-fit result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch multi-fit'
      write(output, '(a,i0)') 'models: ', size(value%fits)
      write(output, '(a,i0)') 'best AIC model: ', value%best_aic
      write(output, '(a,i0)') 'best BIC model: ', value%best_bic
   end subroutine display_rugarch_multifit

   subroutine display_rugarch_parameter_distribution(value, unit)
      !! Display a parameter-distribution summary.
      type(rugarch_parameter_distribution_t), intent(in) :: value !! Draw result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch parameter distribution'
      write(output, '(a,i0)') 'parameters: ', size(value%draws, 1)
      write(output, '(a,i0)') 'draws: ', size(value%draws, 2)
      write(output, '(a,i0)') 'info: ', value%info
   end subroutine display_rugarch_parameter_distribution

   subroutine display_rugarch_multifilter(value, unit)
      !! Display a multi-model filtering summary.
      type(rugarch_multifilter_t), intent(in) :: value !! Multi-filter result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch multi-filter'
      write(output, '(a,i0)') 'models: ', size(value%filters)
      write(output, '(a,i0)') 'info: ', value%info
   end subroutine display_rugarch_multifilter

   subroutine display_rugarch_multiforecast(value, unit)
      !! Display a multi-model forecasting summary.
      type(rugarch_multiforecast_t), intent(in) :: value !! Multi-forecast result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch multi-forecast'
      write(output, '(a,i0)') 'models: ', size(value%forecasts)
      write(output, '(a,i0)') 'info: ', value%info
   end subroutine display_rugarch_multiforecast

   subroutine display_rugarch_model_confidence(value, unit)
      !! Display a model-confidence summary.
      type(rugarch_model_confidence_t), intent(in) :: value !! Confidence set.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch model confidence set'
      write(output, '(a,i0)') 'best model: ', value%best_model
      write(output, '(a,i0)') 'included models: ', count(value%included)
      write(output, '(a,f8.4)') 'level: ', value%level
   end subroutine display_rugarch_model_confidence

   subroutine display_rugarch_news_impact(value, unit)
      !! Display a news-impact curve summary.
      type(rugarch_news_impact_t), intent(in) :: value !! News-impact result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch news-impact curve'
      write(output, '(a,i0)') 'grid points: ', size(value%shock)
      write(output, '(a,es14.6)') 'previous variance: ', value%previous_variance
   end subroutine display_rugarch_news_impact

   subroutine display_rugarch_sign_bias_test(value, unit)
      !! Display Engle-Ng sign and size-bias test results.
      type(rugarch_sign_bias_test_t), intent(in) :: value !! Test result.
      integer, intent(in), optional :: unit !! Output unit.
      integer :: output

      output = output_unit
      if (present(unit)) output = unit
      write(output, '(a)') 'rugarch sign-bias test'
      write(output, '(a,f10.4,a,f10.4)') 'sign statistic: ', &
         value%sign_statistic, '  p-value: ', value%sign_p_value
      write(output, '(a,f10.4,a,f10.4)') 'joint statistic: ', &
         value%joint_statistic, '  p-value: ', value%joint_p_value
   end subroutine display_rugarch_sign_bias_test

end module rugarch_extensions_mod
