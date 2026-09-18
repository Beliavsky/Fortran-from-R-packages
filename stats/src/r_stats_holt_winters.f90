! SPDX-License-Identifier: MIT
! SPDX-FileComment: Holt-Winters filtering and forecasting compatible with R stats conventions.
module r_stats_holt_winters
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use r_distributions, only: r_qnorm
   use r_kinds, only: dp
   use r_stats_decomposition, only: decompose
   use r_stats_types, only: holt_winters_fit_t, holt_winters_forecast_t, &
                            seasonal_decomposition_t
   implicit none
   private

   public :: holt_winters, holt_winters_auto, holt_winters_optimize
   public :: holt_winters_reduced, predict_holt_winters, predict_holt_winters_interval

contains

   pure function holt_winters(x, period, alpha, beta, gamma, initial_level, initial_trend, &
                              initial_season, multiplicative, use_trend, use_season) result(fit)
      !! Applies the R `HoltWinters` recurrences from explicitly supplied initial states.
      real(dp), intent(in) :: x(:) !! Observations, including values preceding the first fit.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      real(dp), intent(in) :: alpha !! Level smoothing parameter in the closed unit interval.
      real(dp), intent(in) :: beta !! Trend smoothing parameter in the closed unit interval.
      real(dp), intent(in) :: gamma !! Seasonal smoothing parameter in the closed unit interval.
      real(dp), intent(in) :: initial_level !! Level preceding the first fitted observation.
      real(dp), intent(in) :: initial_trend !! Trend preceding the first fitted observation.
      real(dp), intent(in) :: initial_season(:) !! Seasonal states with shape `(period)` when used.
      logical, intent(in), optional :: multiplicative !! Select multiplicative seasonality.
      logical, intent(in), optional :: use_trend !! Include a trend component; defaults to true.
      logical, intent(in), optional :: use_season !! Include a seasonal component; defaults to true.
      type(holt_winters_fit_t) :: fit
      real(dp), allocatable :: season_state(:)
      real(dp) :: baseline, level_state, trend_state, new_level, new_trend, new_season
      integer :: i, j, n_fit, start_offset
      logical :: has_season, has_trend, use_multiplicative

      use_multiplicative = .false.
      if (present(multiplicative)) use_multiplicative = multiplicative
      has_trend = .true.
      if (present(use_trend)) has_trend = use_trend
      has_season = .true.
      if (present(use_season)) has_season = use_season
      if (.not. has_season) use_multiplicative = .false.
      fit%period = period
      fit%alpha = alpha
      fit%beta = beta
      fit%gamma = gamma
      fit%has_season = has_season
      fit%has_trend = has_trend
      fit%multiplicative = use_multiplicative
      fit%initial_level = initial_level
      fit%initial_trend = merge(initial_trend, 0.0_dp, has_trend)

      start_offset = 1
      if (has_trend) start_offset = 2
      if (has_season) start_offset = period
      if (period < 1 .or. size(x) <= start_offset .or. &
          (has_season .and. size(initial_season) /= period)) then
         fit%status = 1
         return
      end if
      if (.not. valid_smoothing_parameter(alpha) .or. &
          (has_trend .and. .not. valid_smoothing_parameter(beta)) .or. &
          (has_season .and. .not. valid_smoothing_parameter(gamma)) .or. &
          .not. ieee_is_finite(initial_level) .or. &
          .not. ieee_is_finite(initial_trend) .or. .not. all(ieee_is_finite(x)) .or. &
          (has_season .and. .not. all(ieee_is_finite(initial_season)))) then
         fit%status = 2
         return
      end if
      if (use_multiplicative .and. (any(x == 0.0_dp) .or. any(initial_season == 0.0_dp))) then
         fit%status = 3
         return
      end if

      n_fit = size(x) - start_offset
      allocate (fit%fitted(n_fit), fit%residuals(n_fit), fit%level(n_fit), &
                fit%trend(n_fit), fit%seasonal(n_fit))
      if (has_season) then
         allocate (fit%final_season(period), fit%initial_season(period))
         fit%initial_season = initial_season
         season_state = initial_season
      else
         allocate (fit%final_season(0), fit%initial_season(0), season_state(1), source=0.0_dp)
      end if
      level_state = initial_level
      trend_state = merge(initial_trend, 0.0_dp, has_trend)

      do i = 1, n_fit
         j = 1
         if (has_season) j = modulo(i - 1, period) + 1
         baseline = level_state
         if (has_trend) baseline = baseline + trend_state
         fit%level(i) = level_state
         fit%trend(i) = trend_state
         fit%seasonal(i) = season_state(j)
         if (.not. has_season) then
            fit%fitted(i) = baseline
            new_level = alpha*x(start_offset + i) + (1.0_dp - alpha)*baseline
            new_season = 0.0_dp
         else if (use_multiplicative) then
            fit%fitted(i) = baseline*season_state(j)
            new_level = alpha*x(start_offset + i)/season_state(j) + &
                        (1.0_dp - alpha)*baseline
            new_season = gamma*x(start_offset + i)/new_level + &
                         (1.0_dp - gamma)*season_state(j)
         else
            fit%fitted(i) = baseline + season_state(j)
            new_level = alpha*(x(start_offset + i) - season_state(j)) + &
                        (1.0_dp - alpha)*baseline
            new_season = gamma*(x(start_offset + i) - new_level) + &
                         (1.0_dp - gamma)*season_state(j)
         end if
         new_trend = 0.0_dp
         if (has_trend) then
            new_trend = beta*(new_level - level_state) + (1.0_dp - beta)*trend_state
         end if
         fit%residuals(i) = x(start_offset + i) - fit%fitted(i)
         fit%sse = fit%sse + fit%residuals(i)**2
         level_state = new_level
         trend_state = new_trend
         if (has_season) season_state(j) = new_season
      end do

      fit%final_level = level_state
      fit%final_trend = trend_state
      if (has_season) then
         do i = 1, period
            j = modulo(n_fit + i - 1, period) + 1
            fit%final_season(i) = season_state(j)
         end do
      end if
   end function holt_winters

   pure function holt_winters_reduced(x, alpha, initial_level, beta, &
                                      initial_trend) result(fit)
      !! Fits simple exponential smoothing or Holt's linear trend method.
      real(dp), intent(in) :: x(:) !! Finite observations.
      real(dp), intent(in) :: alpha !! Level smoothing parameter in the closed unit interval.
      real(dp), intent(in), optional :: initial_level !! Initial level; defaults as in R.
      real(dp), intent(in), optional :: beta !! Optional trend smoothing parameter.
      real(dp), intent(in), optional :: initial_trend !! Initial trend; defaults as in R.
      type(holt_winters_fit_t) :: fit
      real(dp) :: beta_value, level_value, trend_value
      logical :: has_trend

      has_trend = present(beta)
      beta_value = 0.0_dp
      if (present(beta)) beta_value = beta
      level_value = 0.0_dp
      trend_value = 0.0_dp
      if (size(x) > 0) level_value = x(1)
      if (has_trend .and. size(x) > 1) then
         level_value = x(2)
         trend_value = x(2) - x(1)
      end if
      if (present(initial_level)) level_value = initial_level
      if (present(initial_trend)) trend_value = initial_trend
      fit = holt_winters(x, 1, alpha, beta_value, 0.0_dp, level_value, trend_value, &
                         [0.0_dp], use_trend=has_trend, use_season=.false.)
   end function holt_winters_reduced

   pure function holt_winters_auto(x, period, alpha, beta, gamma, multiplicative, &
                                   start_periods, use_trend) result(fit)
      !! Applies Holt-Winters filtering after reproducing R's classical initialization.
      real(dp), intent(in) :: x(:) !! Finite observations containing at least two cycles.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      real(dp), intent(in) :: alpha !! Level smoothing parameter in the closed unit interval.
      real(dp), intent(in) :: beta !! Trend smoothing parameter in the closed unit interval.
      real(dp), intent(in) :: gamma !! Seasonal smoothing parameter in the closed unit interval.
      logical, intent(in), optional :: multiplicative !! Select multiplicative seasonality.
      integer, intent(in), optional :: start_periods !! Initial cycles used; defaults to two.
      logical, intent(in), optional :: use_trend !! Include a trend component; defaults to true.
      type(holt_winters_fit_t) :: fit
      real(dp), allocatable :: initial_season(:)
      real(dp) :: initial_level, initial_trend
      integer :: initialization_status, periods
      logical :: has_trend, use_multiplicative

      periods = 2
      if (present(start_periods)) periods = start_periods
      use_multiplicative = .false.
      if (present(multiplicative)) use_multiplicative = multiplicative
      has_trend = .true.
      if (present(use_trend)) has_trend = use_trend
      call initialize_holt_winters(x, period, periods, use_multiplicative, initial_level, &
                                   initial_trend, initial_season, initialization_status)
      if (initialization_status /= 0) then
         fit%period = period
         fit%multiplicative = use_multiplicative
         fit%status = initialization_status
         return
      end if
      fit = holt_winters(x, period, alpha, beta, gamma, initial_level, initial_trend, &
                         initial_season, use_multiplicative, use_trend=has_trend)
   end function holt_winters_auto

   pure function holt_winters_optimize(x, period, multiplicative, start_periods, &
                                       initial_parameters, max_iterations, tolerance) result(fit)
      !! Estimates all smoothing parameters by bounded Nelder-Mead minimization of SSE.
      real(dp), intent(in) :: x(:) !! Finite observations containing at least two cycles.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      logical, intent(in), optional :: multiplicative !! Select multiplicative seasonality.
      integer, intent(in), optional :: start_periods !! Initial cycles used; defaults to two.
      real(dp), intent(in), optional :: initial_parameters(:) !! Initial `(alpha, beta, gamma)`.
      integer, intent(in), optional :: max_iterations !! Iteration limit; defaults to 2000.
      real(dp), intent(in), optional :: tolerance !! Simplex convergence tolerance.
      type(holt_winters_fit_t) :: fit
      real(dp), allocatable :: initial_season(:)
      real(dp) :: best_objective, candidate_objective, convergence_tolerance
      real(dp) :: initial_level, initial_trend, parameters(3), candidate(3), starts(3, 5)
      integer :: i, initialization_status, iteration_limit, iterations, local_iterations, periods
      logical :: converged, local_converged, use_multiplicative

      periods = 2
      if (present(start_periods)) periods = start_periods
      use_multiplicative = .false.
      if (present(multiplicative)) use_multiplicative = multiplicative
      iteration_limit = 2000
      if (present(max_iterations)) iteration_limit = max_iterations
      convergence_tolerance = 1.0e-9_dp
      if (present(tolerance)) convergence_tolerance = tolerance
      parameters = [0.3_dp, 0.1_dp, 0.1_dp]
      if (present(initial_parameters)) then
         if (size(initial_parameters) /= 3) then
            fit%status = 1
            return
         end if
         parameters = initial_parameters
      end if
      if (iteration_limit < 1 .or. convergence_tolerance <= 0.0_dp .or. &
          .not. all(ieee_is_finite(parameters))) then
         fit%status = 2
         return
      end if

      call initialize_holt_winters(x, period, periods, use_multiplicative, initial_level, &
                                   initial_trend, initial_season, initialization_status)
      if (initialization_status /= 0) then
         fit%period = period
         fit%multiplicative = use_multiplicative
         fit%status = initialization_status
         return
      end if
      parameters = bounded_unit(parameters)
      starts(:, 1) = parameters
      starts(:, 2) = [0.9_dp, 0.1_dp, 0.9_dp]
      starts(:, 3) = [0.5_dp, 0.0_dp, 0.9_dp]
      starts(:, 4) = [0.9_dp, 0.5_dp, 0.9_dp]
      starts(:, 5) = [0.5_dp, 0.5_dp, 0.5_dp]
      best_objective = huge(1.0_dp)
      converged = .false.
      iterations = 0
      do i = 1, size(starts, 2)
         candidate = starts(:, i)
         call optimize_smoothing_parameters(x, period, initial_level, initial_trend, &
                                            initial_season, use_multiplicative, candidate, &
                                            iteration_limit, convergence_tolerance, local_iterations, &
                                            local_converged)
         candidate_objective = smoothing_objective(x, period, initial_level, initial_trend, &
                                                   initial_season, use_multiplicative, candidate)
         if (candidate_objective < best_objective) then
            best_objective = candidate_objective
            parameters = candidate
            iterations = local_iterations
            converged = local_converged
         end if
      end do
      fit = holt_winters(x, period, parameters(1), parameters(2), parameters(3), &
                         initial_level, initial_trend, initial_season, use_multiplicative)
      fit%iterations = iterations
      fit%optimized = .true.
      if (.not. converged) fit%status = 5
   end function holt_winters_optimize

   pure subroutine initialize_holt_winters(x, period, start_periods, multiplicative, &
                                           initial_level, initial_trend, initial_season, status)
      !! Derives initial states using R's moving-average decomposition and trend regression.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      integer, intent(in) :: start_periods !! Complete cycles used for initialization.
      logical, intent(in) :: multiplicative !! Whether seasonal ratios replace differences.
      real(dp), intent(out) :: initial_level !! Regressed initial level.
      real(dp), intent(out) :: initial_trend !! Regressed trend increment.
      real(dp), allocatable, intent(out) :: initial_season(:) !! Seasonal states `(period)`.
      integer, intent(out) :: status !! Zero on success or one through four for invalid input.
      type(seasonal_decomposition_t) :: decomposition
      real(dp), allocatable :: trend(:)
      real(dp) :: denominator, mean_index, mean_trend
      integer :: half_width, i, n_trend, window

      status = 0
      initial_level = 0.0_dp
      initial_trend = 0.0_dp
      window = start_periods*period
      if (period < 2 .or. start_periods < 2 .or. window > size(x)) then
         status = 1
         return
      end if
      half_width = period/2
      n_trend = window - 2*half_width
      decomposition = decompose(x(:window), period, multiplicative)
      if (decomposition%status /= 0) then
         status = decomposition%status
         return
      end if
      trend = decomposition%trend(half_width + 1:window - half_width)
      initial_season = decomposition%figure

      mean_index = 0.5_dp*real(n_trend + 1, dp)
      mean_trend = sum(trend)/real(n_trend, dp)
      denominator = 0.0_dp
      do i = 1, n_trend
         denominator = denominator + (real(i, dp) - mean_index)**2
         initial_trend = initial_trend + (real(i, dp) - mean_index)*(trend(i) - mean_trend)
      end do
      if (denominator <= 0.0_dp) then
         status = 4
         return
      end if
      initial_trend = initial_trend/denominator
      initial_level = mean_trend - initial_trend*mean_index
   end subroutine initialize_holt_winters

   pure subroutine optimize_smoothing_parameters(x, period, initial_level, initial_trend, &
                                                 initial_season, multiplicative, parameters, &
                                                 max_iterations, tolerance, iterations, converged)
      !! Minimizes Holt-Winters SSE over the smoothing-parameter unit cube.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      real(dp), intent(in) :: initial_level !! Fixed initial level.
      real(dp), intent(in) :: initial_trend !! Fixed initial trend.
      real(dp), intent(in) :: initial_season(:) !! Fixed initial seasonal states `(period)`.
      logical, intent(in) :: multiplicative !! Whether seasonal effects are multiplicative.
      real(dp), intent(inout) :: parameters(3) !! Initial and optimized smoothing parameters.
      integer, intent(in) :: max_iterations !! Positive iteration limit.
      real(dp), intent(in) :: tolerance !! Positive convergence tolerance.
      integer, intent(out) :: iterations !! Number of completed simplex iterations.
      logical, intent(out) :: converged !! Whether the convergence test was satisfied.
      real(dp) :: centroid(3), contracted(3), expanded(3), reflected(3)
      real(dp) :: simplex(3, 4), values(4)
      real(dp) :: contracted_value, expanded_value, reflected_value
      integer :: i

      simplex(:, 1) = parameters
      do i = 1, 3
         simplex(:, i + 1) = parameters
         simplex(i, i + 1) = min(1.0_dp, parameters(i) + 0.1_dp)
         if (simplex(i, i + 1) == parameters(i)) then
            simplex(i, i + 1) = max(0.0_dp, parameters(i) - 0.1_dp)
         end if
      end do
      do i = 1, 4
         values(i) = smoothing_objective(x, period, initial_level, initial_trend, &
                                         initial_season, multiplicative, simplex(:, i))
      end do

      converged = .false.
      do iterations = 1, max_iterations
         call sort_smoothing_simplex(simplex, values)
         if (maxval(abs(simplex(:, 2:) - spread(simplex(:, 1), 2, 3))) <= tolerance .and. &
             maxval(abs(values(2:) - values(1))) <= tolerance*(1.0_dp + abs(values(1)))) then
            converged = .true.
            exit
         end if
         centroid = sum(simplex(:, :3), dim=2)/3.0_dp
         reflected = bounded_unit(2.0_dp*centroid - simplex(:, 4))
         reflected_value = smoothing_objective(x, period, initial_level, initial_trend, &
                                               initial_season, multiplicative, reflected)
         if (reflected_value < values(1)) then
            expanded = bounded_unit(centroid + 2.0_dp*(reflected - centroid))
            expanded_value = smoothing_objective(x, period, initial_level, initial_trend, &
                                                  initial_season, multiplicative, expanded)
            if (expanded_value < reflected_value) then
               simplex(:, 4) = expanded
               values(4) = expanded_value
            else
               simplex(:, 4) = reflected
               values(4) = reflected_value
            end if
         else if (reflected_value < values(3)) then
            simplex(:, 4) = reflected
            values(4) = reflected_value
         else
            if (reflected_value < values(4)) then
               contracted = bounded_unit(centroid + 0.5_dp*(reflected - centroid))
            else
               contracted = bounded_unit(centroid + 0.5_dp*(simplex(:, 4) - centroid))
            end if
            contracted_value = smoothing_objective(x, period, initial_level, initial_trend, &
                                                    initial_season, multiplicative, contracted)
            if (contracted_value < min(reflected_value, values(4))) then
               simplex(:, 4) = contracted
               values(4) = contracted_value
            else
               do i = 2, 4
                  simplex(:, i) = bounded_unit(simplex(:, 1) + &
                                               0.5_dp*(simplex(:, i) - simplex(:, 1)))
                  values(i) = smoothing_objective(x, period, initial_level, initial_trend, &
                                                  initial_season, multiplicative, simplex(:, i))
               end do
            end if
         end if
      end do
      if (.not. converged) iterations = max_iterations
      call sort_smoothing_simplex(simplex, values)
      parameters = simplex(:, 1)
   end subroutine optimize_smoothing_parameters

   pure real(dp) function smoothing_objective(x, period, initial_level, initial_trend, &
                                              initial_season, multiplicative, parameters)
      !! Returns Holt-Winters SSE for one smoothing-parameter vector.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      real(dp), intent(in) :: initial_level !! Fixed initial level.
      real(dp), intent(in) :: initial_trend !! Fixed initial trend.
      real(dp), intent(in) :: initial_season(:) !! Fixed initial seasonal states `(period)`.
      logical, intent(in) :: multiplicative !! Whether seasonal effects are multiplicative.
      real(dp), intent(in) :: parameters(3) !! Candidate `(alpha, beta, gamma)`.
      real(dp) :: fitted, level_state, new_level, new_season
      real(dp) :: new_trend, season_state(size(initial_season)), trend_state
      integer :: i, j

      smoothing_objective = 0.0_dp
      season_state = initial_season
      level_state = initial_level
      trend_state = initial_trend
      do i = 1, size(x) - period
         j = modulo(i - 1, period) + 1
         if (multiplicative) then
            fitted = (level_state + trend_state)*season_state(j)
            new_level = parameters(1)*x(period + i)/season_state(j) + &
                        (1.0_dp - parameters(1))*(level_state + trend_state)
            new_season = parameters(3)*x(period + i)/new_level + &
                         (1.0_dp - parameters(3))*season_state(j)
         else
            fitted = level_state + trend_state + season_state(j)
            new_level = parameters(1)*(x(period + i) - season_state(j)) + &
                        (1.0_dp - parameters(1))*(level_state + trend_state)
            new_season = parameters(3)*(x(period + i) - new_level) + &
                         (1.0_dp - parameters(3))*season_state(j)
         end if
         new_trend = parameters(2)*(new_level - level_state) + &
                     (1.0_dp - parameters(2))*trend_state
         smoothing_objective = smoothing_objective + (x(period + i) - fitted)**2
         level_state = new_level
         trend_state = new_trend
         season_state(j) = new_season
      end do
   end function smoothing_objective

   pure subroutine sort_smoothing_simplex(simplex, values)
      !! Sorts the four simplex vertices in ascending objective order.
      real(dp), intent(inout) :: simplex(3, 4) !! Parameter vertices arranged by column.
      real(dp), intent(inout) :: values(4) !! Objective values corresponding to vertices.
      real(dp) :: temporary(3), temporary_value
      integer :: i, j

      do i = 2, 4
         j = i
         do while (j > 1)
            if (values(j) >= values(j - 1)) exit
            temporary = simplex(:, j)
            simplex(:, j) = simplex(:, j - 1)
            simplex(:, j - 1) = temporary
            temporary_value = values(j)
            values(j) = values(j - 1)
            values(j - 1) = temporary_value
            j = j - 1
         end do
      end do
   end subroutine sort_smoothing_simplex

   pure function bounded_unit(parameters) result(bounded)
      !! Restricts smoothing parameters to the closed unit interval.
      real(dp), intent(in) :: parameters(:) !! Candidate smoothing parameters.
      real(dp) :: bounded(size(parameters))

      bounded = max(0.0_dp, min(1.0_dp, parameters))
   end function bounded_unit

   pure function predict_holt_winters(fit, n_ahead) result(forecast)
      !! Returns point forecasts from a successful Holt-Winters fit.
      type(holt_winters_fit_t), intent(in) :: fit !! Filtered states and model settings.
      integer, intent(in) :: n_ahead !! Required nonnegative forecast horizon.
      real(dp), allocatable :: forecast(:)
      real(dp) :: baseline
      integer :: i, j

      if (fit%status /= 0 .or. n_ahead < 0 .or. fit%period < 1 .or. &
          .not. allocated(fit%final_season)) then
         allocate (forecast(0))
         return
      end if
      allocate (forecast(n_ahead))
      do i = 1, n_ahead
         baseline = fit%final_level
         if (fit%has_trend) baseline = baseline + real(i, dp)*fit%final_trend
         if (.not. fit%has_season) then
            forecast(i) = baseline
         else
            j = modulo(i - 1, fit%period) + 1
         end if
         if (fit%has_season .and. fit%multiplicative) then
            forecast(i) = baseline*fit%final_season(j)
         else if (fit%has_season) then
            forecast(i) = baseline + fit%final_season(j)
         end if
      end do
   end function predict_holt_winters

   pure function predict_holt_winters_interval(fit, n_ahead, confidence_level) result(prediction)
      !! Returns normal prediction intervals using R's Holt-Winters variance recurrence.
      type(holt_winters_fit_t), intent(in) :: fit !! Filtered states and model settings.
      integer, intent(in) :: n_ahead !! Required nonnegative forecast horizon.
      real(dp), intent(in), optional :: confidence_level !! Central probability; defaults to 0.95.
      type(holt_winters_forecast_t) :: prediction
      real(dp) :: beta, denominator, gamma, mean_residual, probability
      real(dp) :: psi, residual_variance, seasonal_ratio, variance_factor, z_value
      integer :: h, j, relative_season, wrapped_season

      prediction%confidence_level = 0.95_dp
      if (present(confidence_level)) prediction%confidence_level = confidence_level
      if (fit%status /= 0 .or. n_ahead < 0 .or. prediction%confidence_level <= 0.0_dp .or. &
          prediction%confidence_level >= 1.0_dp) then
         prediction%status = 1
         return
      end if
      if (.not. allocated(fit%residuals)) then
         prediction%status = 1
         return
      end if
      if (size(fit%residuals) < 2) then
         prediction%status = 1
         return
      end if
      if (fit%multiplicative .and. (.not. fit%has_trend .or. abs(fit%final_trend) <= tiny(1.0_dp))) then
         prediction%status = 2
         return
      end if

      prediction%fit = predict_holt_winters(fit, n_ahead)
      allocate (prediction%upper(n_ahead), prediction%lower(n_ahead), &
                prediction%standard_error(n_ahead))
      mean_residual = sum(fit%residuals)/real(size(fit%residuals), dp)
      residual_variance = sum((fit%residuals - mean_residual)**2)/ &
                          real(size(fit%residuals) - 1, dp)
      probability = 0.5_dp*(1.0_dp + prediction%confidence_level)
      z_value = r_qnorm(probability)
      beta = merge(fit%beta, 0.0_dp, fit%has_trend)
      gamma = merge(fit%gamma, 0.0_dp, fit%has_season)

      do h = 1, n_ahead
         if (.not. fit%multiplicative) then
            variance_factor = 1.0_dp
            do j = 1, h - 1
               psi = fit%alpha*(1.0_dp + real(j, dp)*beta)
               if (fit%has_season .and. modulo(j, fit%period) == 0) then
                  psi = psi + gamma*(1.0_dp - fit%alpha)
               end if
               variance_factor = variance_factor + psi**2
            end do
         else
            relative_season = modulo(h - 1, fit%period) + 1
            variance_factor = 0.0_dp
            do j = 0, h - 1
               psi = fit%alpha*(1.0_dp + real(j, dp)*beta)
               if (modulo(j, fit%period) == 0) psi = psi + gamma*(1.0_dp - fit%alpha)
               wrapped_season = modulo(relative_season - j, fit%period)
               if (wrapped_season == 0) then
                  denominator = fit%final_trend
               else
                  denominator = fit%final_season(wrapped_season)
               end if
               seasonal_ratio = fit%final_season(relative_season)/denominator
               variance_factor = variance_factor + (psi*seasonal_ratio)**2
            end do
         end if
         prediction%standard_error(h) = sqrt(residual_variance*variance_factor)
         prediction%upper(h) = prediction%fit(h) + z_value*prediction%standard_error(h)
         prediction%lower(h) = prediction%fit(h) - z_value*prediction%standard_error(h)
      end do
   end function predict_holt_winters_interval

   pure elemental logical function valid_smoothing_parameter(value)
      !! Reports whether a finite smoothing parameter lies in the closed unit interval.
      real(dp), intent(in) :: value !! Candidate smoothing parameter.

      valid_smoothing_parameter = ieee_is_finite(value) .and. value >= 0.0_dp .and. value <= 1.0_dp
   end function valid_smoothing_parameter

end module r_stats_holt_winters
