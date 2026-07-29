! SPDX-License-Identifier: GPL-3.0-only
! SPDX-FileComment: Algorithms translated from the R rugarch package.
! Core variance recursions are adapted from the MIT-licensed GARCH-BFGS project.
module rugarch_mod
   !! Univariate ARMA-GARCH filtering, estimation, forecasting, and simulation.
   use kind_mod, only: dp
   use linalg_mod, only: symmetric_pseudoinverse
   use optimization_mod, only: optimization_result_t, bfgs_minimize_fd, &
      finite_difference_hessian
   use random_mod, only: set_random_seed, random_uniform, random_standard_normal, &
      random_standard_student, random_gamma
   use distribution_mod, only: standardized_log_density, random_standardized, &
      shared_distribution_has_skew => distribution_has_skew, &
      shared_distribution_has_shape => distribution_has_shape, &
      distribution_has_lambda, shared_distribution_name => distribution_name
   use special_functions_mod, only: regularized_gamma_q
   use stats_mod, only: variance, normal_quantile
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none
   private

   integer, parameter, public :: rugarch_model_sgarch = 1
   integer, parameter, public :: rugarch_model_igarch = 2
   integer, parameter, public :: rugarch_model_egarch = 3
   integer, parameter, public :: rugarch_model_gjrgarch = 4
   integer, parameter, public :: rugarch_model_aparch = 5
   integer, parameter, public :: rugarch_model_figarch = 6
   integer, parameter, public :: rugarch_model_csgarch = 7
   integer, parameter, public :: rugarch_model_realgarch = 8
   integer, parameter, public :: rugarch_model_fgarch = 9

   integer, parameter, public :: rugarch_fgarch_garch = 1
   integer, parameter, public :: rugarch_fgarch_tgarch = 2
   integer, parameter, public :: rugarch_fgarch_avgarch = 3
   integer, parameter, public :: rugarch_fgarch_ngarch = 4
   integer, parameter, public :: rugarch_fgarch_nagarch = 5
   integer, parameter, public :: rugarch_fgarch_aparch = 6
   integer, parameter, public :: rugarch_fgarch_allgarch = 7
   integer, parameter, public :: rugarch_fgarch_gjrgarch = 8

   integer, parameter, public :: rugarch_distribution_normal = 1
   integer, parameter, public :: rugarch_distribution_student = 2
   integer, parameter, public :: rugarch_distribution_ged = 3
   integer, parameter, public :: rugarch_distribution_skew_normal = 4
   integer, parameter, public :: rugarch_distribution_skew_student = 5
   integer, parameter, public :: rugarch_distribution_skew_ged = 6
   integer, parameter, public :: rugarch_distribution_johnson_su = 7
   integer, parameter, public :: rugarch_distribution_nig = 8
   integer, parameter, public :: rugarch_distribution_ghyp = 9
   integer, parameter, public :: rugarch_distribution_gh_skew_student = 10
   integer, parameter, public :: rugarch_direction_pt = 1
   integer, parameter, public :: rugarch_direction_ag = 2

   type, public :: rugarch_spec_t
      !! Model specification corresponding to the principal ugarchspec options.
      integer :: variance_model = rugarch_model_sgarch
      integer :: distribution = rugarch_distribution_normal
      integer :: arch_order = 1
      integer :: garch_order = 1
      integer :: ar_order = 0
      integer :: ma_order = 0
      integer :: fgarch_submodel = rugarch_fgarch_allgarch
      integer :: mean_regressor_count = 0
      integer :: variance_regressor_count = 0
      logical :: include_mean = .true.
      logical :: fractional_mean = .false.
      logical :: arch_in_mean = .false.
      real(dp) :: arch_in_mean_power = 1.0_dp
      integer :: truncation_lag = 1000
   end type rugarch_spec_t

   type, public :: rugarch_parameters_t
      !! Physical ARMA and conditional-variance parameters.
      real(dp) :: mean = 0.0_dp
      real(dp) :: arch_in_mean = 0.0_dp
      real(dp) :: omega = 0.0_dp
      real(dp) :: power = 2.0_dp
      real(dp) :: shape = 0.0_dp
      real(dp) :: skew = 1.0_dp
      real(dp) :: lambda = -0.5_dp
      real(dp) :: fractional = 0.0_dp
      real(dp) :: mean_fractional = 0.0_dp
      real(dp) :: fgarch_delta = 0.0_dp
      real(dp) :: fgarch_lambda = 2.0_dp
      real(dp) :: component_rho = 0.0_dp
      real(dp) :: component_phi = 0.0_dp
      real(dp) :: measurement_intercept = 0.0_dp
      real(dp) :: measurement_slope = 1.0_dp
      real(dp) :: measurement_leverage1 = 0.0_dp
      real(dp) :: measurement_leverage2 = 0.0_dp
      real(dp) :: measurement_sigma = 0.3_dp
      real(dp), allocatable :: ar(:)
      real(dp), allocatable :: ma(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: asymmetry(:)
      real(dp), allocatable :: fgarch_shift(:)
      real(dp), allocatable :: mean_regression(:)
      real(dp), allocatable :: variance_regression(:)
   end type rugarch_parameters_t

   type, public :: rugarch_filter_t
      !! Conditional mean, variance, innovation, and likelihood paths.
      real(dp), allocatable :: conditional_mean(:)
      real(dp), allocatable :: conditional_variance(:)
      real(dp), allocatable :: conditional_sigma(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: standardized_residuals(:)
      real(dp), allocatable :: log_density(:)
      real(dp), allocatable :: component_variance(:)
      real(dp), allocatable :: realized_variance(:)
      real(dp), allocatable :: measurement_residuals(:)
      real(dp), allocatable :: measurement_log_density(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      integer :: info = 0
   end type rugarch_filter_t

   type, public :: rugarch_fit_t
      !! Estimated rugarch model with inference and filtered paths.
      type(rugarch_spec_t) :: specification
      type(rugarch_parameters_t) :: parameters
      type(rugarch_filter_t) :: filtered
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: scores(:, :)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: hqic = huge(1.0_dp)
      real(dp) :: persistence = 0.0_dp
      real(dp) :: unconditional_variance = 0.0_dp
      real(dp) :: half_life = huge(1.0_dp)
      integer :: iterations = 0
      integer :: optimizer_info = 0
      integer :: info = 0
      logical :: converged = .false.
   end type rugarch_fit_t

   type, public :: rugarch_forecast_t
      !! Recursive conditional mean, variance, and sigma forecasts.
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: realized_variance(:)
      integer :: info = 0
   end type rugarch_forecast_t

   type, public :: rugarch_simulation_t
      !! Simulated observations and latent conditional quantities.
      type(rugarch_spec_t) :: specification
      type(rugarch_parameters_t) :: parameters
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: conditional_mean(:)
      real(dp), allocatable :: conditional_variance(:)
      real(dp), allocatable :: innovations(:)
      real(dp), allocatable :: realized_variance(:)
      real(dp), allocatable :: measurement_errors(:)
      integer :: burnin = 0
      integer :: info = 0
   end type rugarch_simulation_t

   type, public :: rugarch_berkowitz_test_t
      !! Berkowitz likelihood-ratio and Jarque-Bera density calibration tests.
      real(dp), allocatable :: autoregression(:)
      real(dp) :: unrestricted_log_likelihood = -huge(1.0_dp)
      real(dp) :: restricted_log_likelihood = -huge(1.0_dp)
      real(dp) :: likelihood_ratio = 0.0_dp
      real(dp) :: likelihood_ratio_p_value = 1.0_dp
      real(dp) :: jarque_bera = 0.0_dp
      real(dp) :: jarque_bera_p_value = 1.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      integer :: lags = 0
      integer :: info = 0
   end type rugarch_berkowitz_test_t

   type, public :: rugarch_direction_test_t
      !! Pesaran-Timmermann or Anatolyev-Gerko directional-accuracy test.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: directional_accuracy = 0.0_dp
      integer :: method = rugarch_direction_pt
      integer :: info = 0
   end type rugarch_direction_test_t

   type, public :: rugarch_var_test_t
      !! Kupiec unconditional and Christoffersen conditional VaR coverage tests.
      real(dp) :: unconditional_statistic = 0.0_dp
      real(dp) :: unconditional_p_value = 1.0_dp
      real(dp) :: unconditional_critical = 0.0_dp
      real(dp) :: conditional_statistic = 0.0_dp
      real(dp) :: conditional_p_value = 1.0_dp
      real(dp) :: conditional_critical = 0.0_dp
      integer :: expected_exceedances = 0
      integer :: actual_exceedances = 0
      integer :: info = 0
   end type rugarch_var_test_t

   type, public :: rugarch_es_test_t
      !! Conditional expected-shortfall excess-violation test.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: expected_exceedances = 0
      integer :: actual_exceedances = 0
      integer :: info = 0
   end type rugarch_es_test_t

   interface display
      module procedure display_rugarch_spec
      module procedure display_rugarch_parameters
      module procedure display_rugarch_filter
      module procedure display_rugarch_fit
      module procedure display_rugarch_forecast
      module procedure display_rugarch_simulation
      module procedure display_rugarch_berkowitz_test
      module procedure display_rugarch_direction_test
      module procedure display_rugarch_var_test
      module procedure display_rugarch_es_test
   end interface display

   public :: rugarch_spec, rugarch_filter, rugarch_fit
   public :: rugarch_forecast, rugarch_simulate
   public :: rugarch_coefficients, rugarch_persistence
   public :: rugarch_unconditional_variance, rugarch_half_life
   public :: rugarch_log_density, display
   public :: rugarch_berkowitz_test, rugarch_direction_test
   public :: rugarch_var_test, rugarch_es_test
   public :: random_standardized_innovation
   public :: rugarch_fgarch_exponent

contains

   pure function rugarch_spec(variance_model, distribution, arch_order, &
      garch_order, ar_order, ma_order, include_mean, arch_in_mean, &
      arch_in_mean_power, truncation_lag, fgarch_submodel, &
      mean_regressor_count, variance_regressor_count, fractional_mean) &
      result(specification)
      !! Construct the supported numerical subset of a rugarch specification.
      integer, intent(in), optional :: variance_model !! Conditional-variance model code.
      integer, intent(in), optional :: distribution !! Standardized innovation distribution code.
      integer, intent(in), optional :: arch_order !! ARCH order.
      integer, intent(in), optional :: garch_order !! GARCH order.
      integer, intent(in), optional :: ar_order !! Autoregressive mean order.
      integer, intent(in), optional :: ma_order !! Moving-average mean order.
      logical, intent(in), optional :: include_mean !! Include an unconditional mean.
      logical, intent(in), optional :: arch_in_mean !! Include sigma in the conditional mean.
      real(dp), intent(in), optional :: arch_in_mean_power !! Power of sigma used in the mean.
      integer, intent(in), optional :: truncation_lag !! FIGARCH fractional-recursion truncation.
      integer, intent(in), optional :: fgarch_submodel !! Hentschel fGARCH submodel code.
      integer, intent(in), optional :: mean_regressor_count !! Number of external mean regressors.
      integer, intent(in), optional :: variance_regressor_count !! Number of external variance regressors.
      logical, intent(in), optional :: fractional_mean !! Include an ARFIMA fractional mean term.
      type(rugarch_spec_t) :: specification

      if (present(variance_model)) specification%variance_model = variance_model
      if (present(distribution)) specification%distribution = distribution
      if (present(arch_order)) specification%arch_order = arch_order
      if (present(garch_order)) specification%garch_order = garch_order
      if (present(ar_order)) specification%ar_order = ar_order
      if (present(ma_order)) specification%ma_order = ma_order
      if (present(include_mean)) specification%include_mean = include_mean
      if (present(arch_in_mean)) specification%arch_in_mean = arch_in_mean
      if (present(arch_in_mean_power)) &
         specification%arch_in_mean_power = arch_in_mean_power
      if (present(truncation_lag)) specification%truncation_lag = truncation_lag
      if (present(fgarch_submodel)) specification%fgarch_submodel = fgarch_submodel
      if (present(mean_regressor_count)) &
         specification%mean_regressor_count = mean_regressor_count
      if (present(variance_regressor_count)) &
         specification%variance_regressor_count = variance_regressor_count
      if (present(fractional_mean)) specification%fractional_mean = fractional_mean
   end function rugarch_spec

   pure logical function valid_specification(specification) result(valid)
      !! Check model codes, distribution codes, orders, and powers.
      type(rugarch_spec_t), intent(in) :: specification !! Candidate specification.

      valid = specification%variance_model >= rugarch_model_sgarch .and. &
         specification%variance_model <= rugarch_model_fgarch .and. &
         specification%distribution >= rugarch_distribution_normal .and. &
         specification%distribution <= rugarch_distribution_gh_skew_student .and. &
         specification%arch_order >= 0 .and. &
         specification%garch_order >= 0 .and. &
         max(specification%arch_order, specification%garch_order) > 0 .and. &
         specification%ar_order >= 0 .and. specification%ma_order >= 0 .and. &
         specification%arch_in_mean_power > 0.0_dp .and. &
         specification%truncation_lag > 0 .and. &
         specification%mean_regressor_count >= 0 .and. &
         specification%variance_regressor_count >= 0 .and. &
         specification%fgarch_submodel >= rugarch_fgarch_garch .and. &
         specification%fgarch_submodel <= rugarch_fgarch_gjrgarch
      if (specification%variance_model == rugarch_model_egarch .or. &
         specification%variance_model == rugarch_model_gjrgarch .or. &
         specification%variance_model == rugarch_model_aparch .or. &
         specification%variance_model == rugarch_model_fgarch) &
         valid = valid .and. specification%arch_order > 0
      if (specification%variance_model == rugarch_model_figarch) &
         valid = valid .and. specification%arch_order == 1 .and. &
            specification%garch_order == 1
      if (specification%variance_model == rugarch_model_csgarch .or. &
         specification%variance_model == rugarch_model_realgarch) &
         valid = valid .and. specification%arch_order == 1 .and. &
            specification%garch_order == 1
   end function valid_specification

   pure function rugarch_filter(series, specification, parameters, &
      realized_variance, mean_regressors, variance_regressors) result(out)
      !! Filter an ARMA model with a supported rugarch variance recursion.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp), intent(in), optional :: realized_variance(:) !! Realized variance for realGARCH.
      real(dp), intent(in), optional :: mean_regressors(:, :) !! External regressors in the conditional mean.
      real(dp), intent(in), optional :: variance_regressors(:, :) !! External regressors in the variance recursion.
      type(rugarch_filter_t) :: out
      real(dp), allocatable :: powered_scale(:), fractional_weight(:)
      real(dp), allocatable :: mean_fractional_weight(:)
      real(dp) :: backcast, location, variance_value, scale_value
      real(dp) :: component_value, previous_component, previous_variance
      real(dp) :: previous_shock, log_variance, measurement_mean
      real(dp) :: standardized, expected_absolute
      integer :: observations, time, lag, truncation

      observations = size(series)
      if (.not. valid_specification(specification) .or. observations < 5 .or. &
         .not. all(ieee_is_finite(series)) .or. &
         .not. valid_parameters(specification, parameters)) then
         out%info = 1
         return
      end if
      if (specification%variance_model == rugarch_model_realgarch) then
         if (.not. present(realized_variance)) then
            out%info = 1
            return
         end if
         if (size(realized_variance) /= observations .or. &
            any(realized_variance <= 0.0_dp) .or. &
            .not. all(ieee_is_finite(realized_variance))) then
            out%info = 1
            return
         end if
      end if
      if (.not. valid_regressors(observations, specification, &
         mean_regressors, variance_regressors)) then
         out%info = 1
         return
      end if
      allocate(out%conditional_mean(observations), &
         out%conditional_variance(observations), &
         out%conditional_sigma(observations), out%residuals(observations), &
         out%standardized_residuals(observations), &
         out%log_density(observations), powered_scale(observations), &
         out%component_variance(observations), &
         out%realized_variance(observations), &
         out%measurement_residuals(observations), &
         out%measurement_log_density(observations))
      backcast = max(variance(series), 100.0_dp*tiny(1.0_dp))
      out%conditional_mean = 0.0_dp
      out%residuals = 0.0_dp
      out%conditional_variance = backcast
      out%conditional_sigma = sqrt(backcast)
      out%standardized_residuals = 0.0_dp
      out%log_density = 0.0_dp
      out%component_variance = 0.0_dp
      out%realized_variance = 0.0_dp
      out%measurement_residuals = 0.0_dp
      out%measurement_log_density = 0.0_dp
      if (present(realized_variance)) out%realized_variance = realized_variance
      powered_scale = backcast**(0.5_dp*parameters%power)
      previous_component = max(parameters%omega/max(1.0_dp - &
         parameters%component_rho, 1.0e-8_dp), tiny(1.0_dp))
      previous_variance = max(previous_component, backcast)
      previous_shock = backcast
      log_variance = log(backcast)
      if (specification%variance_model == rugarch_model_figarch) then
         truncation = specification%truncation_lag
         fractional_weight = figarch_weights(parameters%alpha(1), &
            parameters%fractional, parameters%beta(1), truncation)
      end if
      if (specification%fractional_mean) then
         allocate(mean_fractional_weight(min(specification%truncation_lag, &
            observations - 1)))
         mean_fractional_weight = fractional_mean_weights( &
            parameters%mean_fractional, size(mean_fractional_weight))
      end if
      expected_absolute = 0.0_dp
      if (specification%variance_model == rugarch_model_egarch) &
         expected_absolute = innovation_absolute_moment( &
         specification%distribution, parameters%shape, parameters%skew, &
         parameters%lambda, 1.0_dp)
      do time = 1, observations
         select case (specification%variance_model)
         case (rugarch_model_sgarch, rugarch_model_igarch, &
            rugarch_model_gjrgarch)
            variance_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  variance_value = variance_value + &
                     parameters%alpha(lag)*out%residuals(time - lag)**2
                  if (specification%variance_model == &
                     rugarch_model_gjrgarch .and. &
                     out%residuals(time - lag) < 0.0_dp) &
                     variance_value = variance_value + &
                        parameters%asymmetry(lag)* &
                        out%residuals(time - lag)**2
               else
                  variance_value = variance_value + &
                     parameters%alpha(lag)*backcast
                  if (specification%variance_model == &
                     rugarch_model_gjrgarch) variance_value = variance_value + &
                     0.5_dp*parameters%asymmetry(lag)*backcast
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  variance_value = variance_value + parameters%beta(lag)* &
                     out%conditional_variance(time - lag)
               else
                  variance_value = variance_value + parameters%beta(lag)*backcast
               end if
            end do
         case (rugarch_model_egarch)
            scale_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  standardized = out%standardized_residuals(time - lag)
                  scale_value = scale_value + parameters%alpha(lag)* &
                     (abs(standardized) - expected_absolute) + &
                     parameters%asymmetry(lag)*standardized
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%beta(lag)* &
                     log(max(out%conditional_variance(time - lag), &
                     tiny(1.0_dp)))
               else
                  scale_value = scale_value + parameters%beta(lag)*log(backcast)
               end if
            end do
            variance_value = exp(max(-50.0_dp, min(50.0_dp, scale_value)))
         case (rugarch_model_aparch)
            scale_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%alpha(lag)* &
                     (abs(out%residuals(time - lag)) - &
                     parameters%asymmetry(lag)*out%residuals(time - lag))** &
                     parameters%power
               else
                  scale_value = scale_value + parameters%alpha(lag)* &
                     innovation_aparch_moment(specification%distribution, &
                     parameters%shape, parameters%skew, parameters%lambda, &
                     parameters%power, &
                     parameters%asymmetry(lag))*powered_scale(1)
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%beta(lag)* &
                     powered_scale(time - lag)
               else
                  scale_value = scale_value + parameters%beta(lag)* &
                     powered_scale(1)
               end if
            end do
            powered_scale(time) = max(scale_value, tiny(1.0_dp))
            variance_value = powered_scale(time)**(2.0_dp/parameters%power)
         case (rugarch_model_figarch)
            variance_value = parameters%omega/max(1.0_dp - &
               parameters%beta(1), 1.0e-8_dp)
            if (time <= size(fractional_weight)) variance_value = &
               variance_value + sum(fractional_weight(time:))*backcast
            do lag = 1, min(time - 1, size(fractional_weight))
                variance_value = variance_value + fractional_weight(lag)* &
                   out%residuals(time - lag)**2
             end do
         case (rugarch_model_csgarch)
            component_value = parameters%omega + &
               parameters%component_rho*previous_component + &
               parameters%component_phi*(previous_shock - previous_variance)
            component_value = max(component_value, tiny(1.0_dp))
            variance_value = component_value + parameters%alpha(1)* &
               (previous_shock - previous_component) + parameters%beta(1)* &
               (previous_variance - previous_component)
            variance_value = max(variance_value, tiny(1.0_dp))
            out%component_variance(time) = component_value
         case (rugarch_model_realgarch)
            if (time > 1) log_variance = parameters%omega + &
               parameters%alpha(1)*log(realized_variance(time - 1)) + &
               parameters%beta(1)*log_variance
            variance_value = exp(max(-50.0_dp, min(50.0_dp, log_variance)))
         case (rugarch_model_fgarch)
            scale_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  standardized = out%standardized_residuals(time - lag) - &
                     parameters%fgarch_shift(lag)
                  scale_value = scale_value + parameters%alpha(lag)* &
                     (sqrt(1.0e-6_dp + standardized**2) - &
                     parameters%asymmetry(lag)*standardized)** &
                     max(parameters%fgarch_delta + fgarch_fk( &
                     specification%fgarch_submodel)*parameters%fgarch_lambda, &
                     0.01_dp)*out%conditional_sigma(time - lag)** &
                     parameters%fgarch_lambda
               else
                  scale_value = scale_value + parameters%alpha(lag)* &
                     fgarch_innovation_moment(specification, parameters, lag)* &
                     backcast**(0.5_dp*parameters%fgarch_lambda)
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%beta(lag)* &
                     out%conditional_sigma(time - lag)** &
                     parameters%fgarch_lambda
               else
                  scale_value = scale_value + parameters%beta(lag)* &
                     backcast**(0.5_dp*parameters%fgarch_lambda)
               end if
            end do
            variance_value = max(scale_value, tiny(1.0_dp))** &
               (2.0_dp/parameters%fgarch_lambda)
         end select
         if (specification%variance_regressor_count > 0) then
            if (specification%variance_model == rugarch_model_egarch .or. &
               specification%variance_model == rugarch_model_realgarch) then
               variance_value = variance_value*exp(dot_product( &
                  variance_regressors(time, :), &
                  parameters%variance_regression))
            else
               variance_value = variance_value + dot_product( &
                  variance_regressors(time, :), &
                  parameters%variance_regression)
            end if
         end if
         if (.not. ieee_is_finite(variance_value) .or. &
            variance_value <= tiny(1.0_dp)) then
            out%info = 2
            out%log_likelihood = -huge(1.0_dp)
            return
         end if
         out%conditional_variance(time) = variance_value
         out%conditional_sigma(time) = sqrt(variance_value)
         location = merge(parameters%mean, 0.0_dp, specification%include_mean)
         do lag = 1, specification%ar_order
            if (time - lag >= 1) location = location + parameters%ar(lag)* &
               (series(time - lag) - merge(parameters%mean, 0.0_dp, &
               specification%include_mean))
         end do
         do lag = 1, specification%ma_order
            if (time - lag >= 1) location = location + &
               parameters%ma(lag)*out%residuals(time - lag)
         end do
         if (specification%fractional_mean) then
            do lag = 1, min(time - 1, size(mean_fractional_weight))
               location = location - mean_fractional_weight(lag)* &
                  (series(time - lag) - merge(parameters%mean, 0.0_dp, &
                  specification%include_mean))
            end do
         end if
         if (specification%mean_regressor_count > 0) location = location + &
            dot_product(mean_regressors(time, :), parameters%mean_regression)
         if (specification%arch_in_mean) location = location + &
            parameters%arch_in_mean*out%conditional_sigma(time)** &
            specification%arch_in_mean_power
         out%conditional_mean(time) = location
         out%residuals(time) = series(time) - location
         out%standardized_residuals(time) = out%residuals(time)/ &
            out%conditional_sigma(time)
         out%log_density(time) = rugarch_log_density( &
            out%standardized_residuals(time), specification%distribution, &
            parameters%shape, parameters%skew, parameters%lambda) - &
            log(out%conditional_sigma(time))
         if (specification%variance_model == rugarch_model_realgarch) then
            measurement_mean = parameters%measurement_intercept + &
               parameters%measurement_slope*log_variance + &
               parameters%measurement_leverage1*out%standardized_residuals(time) + &
               parameters%measurement_leverage2* &
               (out%standardized_residuals(time)**2 - 1.0_dp)
            out%measurement_residuals(time) = log(realized_variance(time)) - &
               measurement_mean
            out%measurement_log_density(time) = -0.5_dp* &
               (log(2.0_dp*acos(-1.0_dp)) + &
               2.0_dp*log(parameters%measurement_sigma) + &
               (out%measurement_residuals(time)/ &
               parameters%measurement_sigma)**2)
         end if
         if (specification%variance_model == rugarch_model_csgarch) then
            previous_component = component_value
            previous_variance = variance_value
            previous_shock = out%residuals(time)**2
         end if
      end do
      out%log_likelihood = sum(out%log_density) + &
         sum(out%measurement_log_density)
   end function rugarch_filter

   pure logical function valid_regressors(observations, specification, &
      mean_regressors, variance_regressors) result(valid)
      !! Check optional external-regressor presence and dimensions.
      integer, intent(in) :: observations !! Number of observations.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      real(dp), intent(in), optional :: mean_regressors(:, :) !! Mean regressors.
      real(dp), intent(in), optional :: variance_regressors(:, :) !! Variance regressors.

      valid = specification%mean_regressor_count == 0
      if (present(mean_regressors)) valid = &
         size(mean_regressors, 1) == observations .and. &
         size(mean_regressors, 2) == specification%mean_regressor_count .and. &
         all(ieee_is_finite(mean_regressors))
      if (.not. valid) return
      valid = specification%variance_regressor_count == 0
      if (present(variance_regressors)) valid = &
         size(variance_regressors, 1) == observations .and. &
         size(variance_regressors, 2) == &
         specification%variance_regressor_count .and. &
         all(ieee_is_finite(variance_regressors))
   end function valid_regressors

   pure function fractional_mean_weights(fractional, truncation) result(weights)
      !! Return coefficients of the fractional-difference polynomial.
      real(dp), intent(in) :: fractional !! Fractional mean parameter.
      integer, intent(in) :: truncation !! Number of retained coefficients.
      real(dp) :: weights(truncation)
      integer :: lag

      if (truncation < 1) return
      weights(1) = -fractional
      do lag = 2, truncation
         weights(lag) = weights(lag - 1)* &
            (real(lag - 1, dp) - fractional)/real(lag, dp)
      end do
   end function fractional_mean_weights

   pure logical function valid_parameters(specification, parameters) &
      result(valid)
      !! Validate physical parameter dimensions and distributional constraints.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp) :: persistence

      valid = allocated(parameters%ar) .and. allocated(parameters%ma) .and. &
         allocated(parameters%alpha) .and. allocated(parameters%beta) .and. &
         allocated(parameters%asymmetry)
      if (.not. valid) return
      if (specification%mean_regressor_count > 0) valid = &
         allocated(parameters%mean_regression)
      if (valid .and. specification%variance_regressor_count > 0) valid = &
         allocated(parameters%variance_regression)
      if (valid .and. specification%variance_model == rugarch_model_fgarch) &
         valid = allocated(parameters%fgarch_shift)
      if (.not. valid) return
      if (specification%mean_regressor_count > 0) valid = &
         size(parameters%mean_regression) == specification%mean_regressor_count
      if (valid .and. specification%variance_regressor_count > 0) valid = &
         size(parameters%variance_regression) == &
         specification%variance_regressor_count
      if (valid .and. specification%variance_model == rugarch_model_fgarch) &
         valid = size(parameters%fgarch_shift) == specification%arch_order
      if (.not. valid) return
      valid = size(parameters%ar) == specification%ar_order .and. &
         size(parameters%ma) == specification%ma_order .and. &
         size(parameters%alpha) == specification%arch_order .and. &
         size(parameters%beta) == specification%garch_order .and. &
         size(parameters%asymmetry) == specification%arch_order .and. &
         all(parameters%alpha >= 0.0_dp) .and. &
         all(parameters%beta >= 0.0_dp)
      if (.not. valid) return
      if (specification%variance_model == rugarch_model_egarch .or. &
         specification%variance_model == rugarch_model_realgarch) then
         valid = ieee_is_finite(parameters%omega)
      else
         valid = parameters%omega > 0.0_dp
      end if
      if (.not. valid) return
      if (specification%distribution == rugarch_distribution_student) &
         valid = parameters%shape > 2.0_dp
      if (specification%distribution == rugarch_distribution_ged) &
         valid = parameters%shape > 0.0_dp
      if (specification%distribution == rugarch_distribution_skew_normal) &
         valid = parameters%skew > 0.0_dp
      if (specification%distribution == rugarch_distribution_skew_student) &
         valid = parameters%shape > 2.0_dp .and. parameters%skew > 0.0_dp
      if (specification%distribution == rugarch_distribution_skew_ged) &
         valid = parameters%shape > 0.0_dp .and. parameters%skew > 0.0_dp
      if (specification%distribution == rugarch_distribution_johnson_su) &
         valid = parameters%shape > 0.0_dp .and. &
            ieee_is_finite(parameters%skew)
      if (specification%distribution == rugarch_distribution_nig .or. &
         specification%distribution == rugarch_distribution_ghyp) &
         valid = parameters%shape > 0.0_dp .and. abs(parameters%skew) < 1.0_dp
      if (specification%distribution == rugarch_distribution_gh_skew_student) &
         valid = parameters%shape > 4.0_dp .and. ieee_is_finite(parameters%skew)
      if (specification%distribution == rugarch_distribution_ghyp) &
         valid = valid .and. ieee_is_finite(parameters%lambda)
      if (specification%variance_model == rugarch_model_gjrgarch) &
         valid = valid .and. all(parameters%alpha + parameters%asymmetry >= 0.0_dp)
      if (specification%variance_model == rugarch_model_aparch) &
         valid = valid .and. parameters%power > 0.0_dp .and. &
            all(abs(parameters%asymmetry) < 1.0_dp)
      if (specification%variance_model == rugarch_model_figarch) &
         valid = valid .and. parameters%fractional > 0.0_dp .and. &
            parameters%fractional < 1.0_dp .and. &
            parameters%alpha(1) <= 0.5_dp*(1.0_dp - &
            parameters%fractional) .and. &
            parameters%beta(1) <= parameters%fractional + &
             parameters%alpha(1)
      if (specification%variance_model == rugarch_model_csgarch) &
         valid = valid .and. parameters%component_rho > 0.0_dp .and. &
            parameters%component_rho < 1.0_dp .and. &
            parameters%component_phi >= 0.0_dp .and. &
            parameters%component_phi <= 0.25_dp
      if (specification%variance_model == rugarch_model_realgarch) &
         valid = valid .and. parameters%measurement_slope > 0.0_dp .and. &
            parameters%measurement_sigma > 0.0_dp
      if (specification%fractional_mean) valid = valid .and. &
         parameters%mean_fractional > -0.5_dp .and. &
         parameters%mean_fractional < 0.5_dp
      if (specification%variance_model == rugarch_model_fgarch) &
         valid = valid .and. parameters%fgarch_lambda > 0.0_dp .and. &
            parameters%fgarch_delta + parameters%fgarch_lambda > 0.0_dp .and. &
            all(abs(parameters%asymmetry) < 1.0_dp)
      persistence = rugarch_persistence(specification, parameters)
      if (specification%variance_model /= rugarch_model_egarch .and. &
         specification%variance_model /= rugarch_model_igarch .and. &
         specification%variance_model /= rugarch_model_figarch .and. &
         specification%variance_model /= rugarch_model_fgarch) &
         valid = valid .and. persistence < 1.0_dp
      if (specification%variance_model == rugarch_model_egarch) &
         valid = valid .and. sum(parameters%beta) < 1.0_dp
   end function valid_parameters

   pure elemental real(dp) function rugarch_log_density(value, distribution, &
      shape, skew, lambda) result(log_density)
      !! Evaluate a standardized rugarch innovation log density.
      real(dp), intent(in) :: value !! Standardized innovation.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Distribution shape parameter.
      real(dp), intent(in), optional :: skew !! Distribution skew parameter.
      real(dp), intent(in), optional :: lambda !! Generalized-hyperbolic lambda.
      real(dp) :: scale, skew_value, first_moment, raw_mean, raw_second
      real(dp) :: raw_scale, transformed, side_scale, inverse_shape
      real(dp) :: weight, omega, log_standardizer, normal_value
      real(dp) :: lambda_value

      skew_value = 1.0_dp
      if (present(skew)) skew_value = skew
      lambda_value = -0.5_dp
      if (present(lambda)) lambda_value = lambda
      log_density = standardized_log_density(value, distribution, shape, &
         skew_value, lambda_value)
      return

      select case (distribution)
      case (rugarch_distribution_normal)
         log_density = -0.5_dp*(log(2.0_dp*acos(-1.0_dp)) + value**2)
      case (rugarch_distribution_student)
         if (shape <= 2.0_dp) then
            log_density = -huge(1.0_dp)
         else
            log_density = log_gamma(0.5_dp*(shape + 1.0_dp)) - &
               log_gamma(0.5_dp*shape) - &
               0.5_dp*log(acos(-1.0_dp)*(shape - 2.0_dp)) - &
               0.5_dp*(shape + 1.0_dp)* &
               log(1.0_dp + value**2/(shape - 2.0_dp))
         end if
      case (rugarch_distribution_ged)
         if (shape <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            scale = exp(0.5_dp*(log_gamma(1.0_dp/shape) - &
               log_gamma(3.0_dp/shape)))
             log_density = log(shape) - log(2.0_dp) - log(scale) - &
                log_gamma(1.0_dp/shape) - (abs(value)/scale)**shape
          end if
      case (rugarch_distribution_skew_normal)
         first_moment = sqrt(2.0_dp/acos(-1.0_dp))
         call fs_standardization(skew_value, first_moment, raw_mean, raw_scale)
         transformed = value*raw_scale + raw_mean
         side_scale = merge(skew_value, 1.0_dp/skew_value, transformed >= 0.0_dp)
         log_density = log(2.0_dp/(skew_value + 1.0_dp/skew_value)) + &
            log(raw_scale) - 0.5_dp*log(2.0_dp*acos(-1.0_dp)) - &
            0.5_dp*(transformed/side_scale)**2
      case (rugarch_distribution_skew_student)
         if (shape <= 2.0_dp .or. skew_value <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            first_moment = exp(0.5_dp*log(shape - 2.0_dp) + &
               log_gamma(0.5_dp*(shape - 1.0_dp)) - &
               0.5_dp*log(acos(-1.0_dp)) - log_gamma(0.5_dp*shape))
            call fs_standardization(skew_value, first_moment, raw_mean, raw_scale)
            transformed = value*raw_scale + raw_mean
            side_scale = merge(skew_value, 1.0_dp/skew_value, &
               transformed >= 0.0_dp)
            log_density = log(2.0_dp/(skew_value + 1.0_dp/skew_value)) + &
               log(raw_scale) + log_gamma(0.5_dp*(shape + 1.0_dp)) - &
               log_gamma(0.5_dp*shape) - &
               0.5_dp*log(acos(-1.0_dp)*(shape - 2.0_dp)) - &
               0.5_dp*(shape + 1.0_dp)*log(1.0_dp + &
               (transformed/side_scale)**2/(shape - 2.0_dp))
         end if
      case (rugarch_distribution_skew_ged)
         if (shape <= 0.0_dp .or. skew_value <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            scale = exp(0.5_dp*(log_gamma(1.0_dp/shape) - &
               log_gamma(3.0_dp/shape)))
            first_moment = scale*exp(log_gamma(2.0_dp/shape) - &
               log_gamma(1.0_dp/shape))
            call fs_standardization(skew_value, first_moment, raw_mean, raw_scale)
            transformed = value*raw_scale + raw_mean
            side_scale = merge(skew_value, 1.0_dp/skew_value, &
               transformed >= 0.0_dp)
            log_density = log(2.0_dp/(skew_value + 1.0_dp/skew_value)) + &
               log(raw_scale) + log(shape) - log(2.0_dp) - log(scale) - &
               log_gamma(1.0_dp/shape) - &
               (abs(transformed/side_scale)/scale)**shape
         end if
      case (rugarch_distribution_johnson_su)
         if (shape <= 0.0_dp) then
            log_density = -huge(1.0_dp)
         else
            inverse_shape = 1.0_dp/shape
            weight = exp(min(inverse_shape**2, 50.0_dp))
            omega = -skew_value*inverse_shape
            log_standardizer = -0.5_dp*log(0.5_dp*(weight - 1.0_dp)* &
               (weight*cosh(2.0_dp*omega) + 1.0_dp))
            scale = exp(log_standardizer)
            transformed = (value - scale*sqrt(weight)*sinh(omega))/scale
            normal_value = -skew_value + asinh(transformed)/inverse_shape
            log_density = -log(scale) - log(inverse_shape) - &
               0.5_dp*log(1.0_dp + transformed**2) - &
               0.5_dp*log(2.0_dp*acos(-1.0_dp)) - 0.5_dp*normal_value**2
         end if
      case default
         log_density = -huge(1.0_dp)
      end select
   end function rugarch_log_density

   pure subroutine fs_standardization(skew, first_moment, mean, scale)
      !! Return the mean and scale of a Fernandez-Steel transformed variate.
      real(dp), intent(in) :: skew !! Positive Fernandez-Steel skew parameter.
      real(dp), intent(in) :: first_moment !! Absolute first moment of the symmetric base.
      real(dp), intent(out) :: mean !! Raw transformed mean.
      real(dp), intent(out) :: scale !! Raw transformed standard deviation.
      real(dp) :: second_moment

      mean = first_moment*(skew - 1.0_dp/skew)
      second_moment = (skew**3 + skew**(-3))/ &
         (skew + 1.0_dp/skew)
      scale = sqrt(max(second_moment - mean**2, tiny(1.0_dp)))
   end subroutine fs_standardization

   pure real(dp) function innovation_absolute_moment(distribution, shape, skew, &
      lambda, power) &
      result(moment)
      !! Return an absolute moment of a standardized symmetric innovation.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Distribution shape parameter.
      real(dp), intent(in) :: skew !! Distribution skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda.
      real(dp), intent(in) :: power !! Positive absolute-moment power.
      real(dp) :: scale

      select case (distribution)
      case (rugarch_distribution_normal)
         moment = exp(0.5_dp*power*log(2.0_dp) + &
            log_gamma(0.5_dp*(power + 1.0_dp)) - &
            0.5_dp*log(acos(-1.0_dp)))
      case (rugarch_distribution_student)
         if (shape <= max(2.0_dp, power)) then
            moment = huge(1.0_dp)
         else
            moment = exp(0.5_dp*power*log(shape - 2.0_dp) + &
               log_gamma(0.5_dp*(power + 1.0_dp)) + &
               log_gamma(0.5_dp*(shape - power)) - &
               0.5_dp*log(acos(-1.0_dp)) - log_gamma(0.5_dp*shape))
         end if
      case (rugarch_distribution_ged)
         scale = exp(0.5_dp*(log_gamma(1.0_dp/shape) - &
            log_gamma(3.0_dp/shape)))
          moment = scale**power*exp(log_gamma((power + 1.0_dp)/shape) - &
             log_gamma(1.0_dp/shape))
      case (rugarch_distribution_skew_student)
         if (shape <= max(2.0_dp, power)) then
            moment = huge(1.0_dp)
         else
            moment = numerical_innovation_moment(distribution, shape, skew, &
               lambda, power, 0.0_dp)
         end if
      case (rugarch_distribution_skew_normal, rugarch_distribution_skew_ged, &
         rugarch_distribution_johnson_su, rugarch_distribution_nig, &
         rugarch_distribution_ghyp)
         moment = numerical_innovation_moment(distribution, shape, skew, lambda, &
            power, 0.0_dp)
      case (rugarch_distribution_gh_skew_student)
         if (shape <= max(4.0_dp, power)) then
            moment = huge(1.0_dp)
         else
            moment = numerical_innovation_moment(distribution, shape, skew, &
               lambda, power, 0.0_dp)
         end if
      case default
         moment = huge(1.0_dp)
      end select
   end function innovation_absolute_moment

   pure real(dp) function innovation_aparch_moment(distribution, shape, skew, &
      lambda, power, asymmetry) result(moment)
      !! Return E[(abs(z)-gamma*z)^delta] for a symmetric innovation.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Distribution shape parameter.
      real(dp), intent(in) :: skew !! Distribution skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda.
      real(dp), intent(in) :: power !! APARCH power delta.
      real(dp), intent(in) :: asymmetry !! APARCH leverage parameter gamma.

      if (distribution <= rugarch_distribution_ged) then
         moment = 0.5_dp*((1.0_dp - asymmetry)**power + &
            (1.0_dp + asymmetry)**power)* &
            innovation_absolute_moment(distribution, shape, skew, lambda, power)
      else
         moment = numerical_innovation_moment(distribution, shape, skew, lambda, &
            power, asymmetry)
      end if
   end function innovation_aparch_moment

   pure real(dp) function numerical_innovation_moment(distribution, shape, &
      skew, lambda, power, asymmetry) result(moment)
      !! Integrate a standardized innovation power moment over the real line.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Distribution shape parameter.
      real(dp), intent(in) :: skew !! Distribution skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda.
      real(dp), intent(in) :: power !! Positive moment power.
      real(dp), intent(in) :: asymmetry !! APARCH leverage parameter.
      integer, parameter :: intervals = 512
      real(dp) :: angle, value, jacobian
      integer :: index

      moment = 0.0_dp
      do index = 1, intervals
         angle = acos(-1.0_dp)*( &
            (real(index, dp) - 0.5_dp)/real(intervals, dp) - 0.5_dp)
         value = tan(angle)
         jacobian = acos(-1.0_dp)/cos(angle)**2
         moment = moment + (abs(value) - asymmetry*value)**power* &
            exp(rugarch_log_density(value, distribution, shape, skew, lambda))* &
            jacobian
      end do
      moment = moment/real(intervals, dp)
   end function numerical_innovation_moment

   pure function rugarch_fit(series, specification, initial, max_iterations, &
      tolerance, realized_variance, mean_regressors, variance_regressors) &
      result(out)
      !! Translate ugarchfit using shared finite-difference BFGS optimization.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in), optional :: initial !! Optional physical starting values.
      integer, intent(in), optional :: max_iterations !! Maximum BFGS iterations.
      real(dp), intent(in), optional :: tolerance !! Optimizer gradient tolerance.
      real(dp), intent(in), optional :: realized_variance(:) !! Realized variance for realGARCH.
      real(dp), intent(in), optional :: mean_regressors(:, :) !! External mean regressors.
      real(dp), intent(in), optional :: variance_regressors(:, :) !! External variance regressors.
      type(rugarch_fit_t) :: out
      type(optimization_result_t) :: optimum
      type(rugarch_parameters_t) :: starting, selected
      type(rugarch_filter_t) :: filtered
      real(dp), allocatable :: raw_start(:), hessian(:, :), inverse(:, :)
      real(dp), allocatable :: jacobian(:, :), physical_covariance(:, :)
      real(dp) :: selected_tolerance
      integer :: limit, count, status, parameter

      out%specification = specification
      limit = 500
      if (present(max_iterations)) limit = max_iterations
      selected_tolerance = 1.0e-5_dp
      if (present(tolerance)) selected_tolerance = tolerance
      if (.not. valid_specification(specification) .or. size(series) < 20 .or. &
         .not. all(ieee_is_finite(series)) .or. limit < 1 .or. &
         selected_tolerance <= 0.0_dp) then
         out%info = 1
         return
      end if
      if (specification%variance_model == rugarch_model_realgarch) then
         if (.not. present(realized_variance)) then
            out%info = 1
            return
         end if
         if (size(realized_variance) /= size(series)) then
            out%info = 1
            return
         end if
      end if
      if (.not. valid_regressors(size(series), specification, &
         mean_regressors, variance_regressors)) then
         out%info = 1
         return
      end if
      if (present(initial)) then
         starting = initial
         if (.not. valid_parameters(specification, starting)) then
            out%info = 2
            return
         end if
      else
         starting = default_parameters(series, specification)
      end if
      raw_start = parameters_to_raw(specification, starting)
      optimum = bfgs_minimize_fd(objective, raw_start, limit, selected_tolerance)
      if (.not. allocated(optimum%parameters)) then
         out%info = 3
         return
      end if
      selected = raw_to_parameters(specification, optimum%parameters)
      filtered = rugarch_filter(series, specification, selected, &
         realized_variance, mean_regressors, variance_regressors)
      if (filtered%info /= 0 .or. .not. ieee_is_finite(filtered%log_likelihood)) then
         selected = starting
         filtered = rugarch_filter(series, specification, selected, &
            realized_variance, mean_regressors, variance_regressors)
         optimum%parameters = raw_start
         optimum%converged = .false.
      end if
      if (filtered%info /= 0) then
         out%info = 4
         return
      end if
      out%parameters = selected
      out%filtered = filtered
      out%coefficients = rugarch_coefficients(specification, selected)
      out%log_likelihood = filtered%log_likelihood
      out%persistence = rugarch_persistence(specification, selected)
      out%unconditional_variance = rugarch_unconditional_variance( &
         specification, selected)
      out%half_life = rugarch_half_life(specification, selected)
      out%iterations = optimum%iterations
      out%optimizer_info = optimum%info
      out%converged = optimum%converged
      out%scores = likelihood_scores(optimum%parameters)
      count = size(out%coefficients)
      out%aic = -2.0_dp*out%log_likelihood + 2.0_dp*real(count, dp)
      out%bic = -2.0_dp*out%log_likelihood + &
         log(real(size(series), dp))*real(count, dp)
      out%hqic = -2.0_dp*out%log_likelihood + &
         2.0_dp*log(log(real(size(series), dp)))*real(count, dp)
      hessian = finite_difference_hessian(objective, optimum%parameters)
      allocate(inverse(size(hessian, 1), size(hessian, 2)))
      call symmetric_pseudoinverse(hessian, inverse, status)
      if (status == 0) then
         jacobian = physical_parameter_jacobian(specification, &
            optimum%parameters)
         physical_covariance = matmul(matmul(jacobian, inverse), &
            transpose(jacobian))
         out%covariance = physical_covariance
         allocate(out%standard_errors(count))
         do parameter = 1, count
            out%standard_errors(parameter) = &
               sqrt(max(physical_covariance(parameter, parameter), 0.0_dp))
         end do
      else
         allocate(out%covariance(0, 0), out%standard_errors(0))
      end if

   contains

      pure real(dp) function objective(raw) result(value)
         !! Evaluate the negative conditional log likelihood.
         real(dp), intent(in) :: raw(:) !! Unconstrained optimizer parameters.
         type(rugarch_parameters_t) :: candidate
         type(rugarch_filter_t) :: candidate_filter

         candidate = raw_to_parameters(specification, raw)
         candidate_filter = rugarch_filter(series, specification, candidate, &
            realized_variance, mean_regressors, variance_regressors)
         if (candidate_filter%info /= 0 .or. &
            .not. ieee_is_finite(candidate_filter%log_likelihood)) then
            value = huge(1.0_dp)
         else
            value = -candidate_filter%log_likelihood
            value = value + arma_constraint_penalty(candidate)
         end if
      end function objective

      pure function likelihood_scores(raw) result(scores)
         !! Numerically differentiate observation log likelihood contributions.
         real(dp), intent(in) :: raw(:) !! Unconstrained optimizer parameters.
         real(dp), allocatable :: scores(:, :)
         real(dp) :: shifted(size(raw)), step
         type(rugarch_filter_t) :: upper_filter, lower_filter
         integer :: parameter

         allocate(scores(size(series), size(raw)))
         scores = 0.0_dp
         do parameter = 1, size(raw)
            step = epsilon(1.0_dp)**(1.0_dp/3.0_dp)* &
               max(1.0_dp, abs(raw(parameter)))
            shifted = raw
            shifted(parameter) = raw(parameter) + step
            upper_filter = rugarch_filter(series, specification, &
               raw_to_parameters(specification, shifted), realized_variance, &
               mean_regressors, variance_regressors)
            shifted(parameter) = raw(parameter) - step
            lower_filter = rugarch_filter(series, specification, &
               raw_to_parameters(specification, shifted), realized_variance, &
               mean_regressors, variance_regressors)
            if (upper_filter%info == 0 .and. lower_filter%info == 0) &
               scores(:, parameter) = (upper_filter%log_density + &
               upper_filter%measurement_log_density - &
               lower_filter%log_density - &
               lower_filter%measurement_log_density)/(2.0_dp*step)
         end do
      end function likelihood_scores

   end function rugarch_fit

   pure function default_parameters(series, specification) result(parameters)
      !! Construct stable data-scaled starting parameters.
      real(dp), intent(in) :: series(:) !! Time-series observations.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t) :: parameters
      real(dp) :: variance_value, arch_total, garch_total

      allocate(parameters%ar(specification%ar_order), &
         parameters%ma(specification%ma_order), &
         parameters%alpha(specification%arch_order), &
         parameters%beta(specification%garch_order), &
         parameters%asymmetry(specification%arch_order), &
         parameters%fgarch_shift(specification%arch_order), &
         parameters%mean_regression(specification%mean_regressor_count), &
         parameters%variance_regression( &
         specification%variance_regressor_count))
      parameters%ar = 0.0_dp
      parameters%ma = 0.0_dp
      parameters%asymmetry = 0.0_dp
      parameters%fgarch_shift = 0.0_dp
      parameters%mean_regression = 0.0_dp
      parameters%variance_regression = 0.0_dp
      parameters%mean = sum(series)/real(size(series), dp)
      if (.not. specification%include_mean) parameters%mean = 0.0_dp
      parameters%arch_in_mean = 0.0_dp
      parameters%power = 2.0_dp
      parameters%fractional = 0.0_dp
      parameters%mean_fractional = 0.10_dp
      parameters%fgarch_delta = 0.0_dp
      parameters%fgarch_lambda = 2.0_dp
      parameters%skew = 1.0_dp
      parameters%lambda = -0.5_dp
      variance_value = max(variance(series), 100.0_dp*tiny(1.0_dp))
      arch_total = merge(0.08_dp, 0.0_dp, specification%arch_order > 0)
      garch_total = merge(0.90_dp, 0.0_dp, specification%garch_order > 0)
      if (specification%garch_order == 0) arch_total = 0.90_dp
      if (specification%arch_order == 0) garch_total = 0.95_dp
      if (specification%variance_model == rugarch_model_igarch) then
         if (specification%arch_order > 0 .and. &
            specification%garch_order > 0) then
            arch_total = 0.08_dp
            garch_total = 0.92_dp
         else if (specification%arch_order > 0) then
            arch_total = 1.0_dp
         else
            garch_total = 1.0_dp
         end if
      end if
      if (specification%arch_order > 0) parameters%alpha = &
         arch_total/real(specification%arch_order, dp)
      if (specification%garch_order > 0) parameters%beta = &
         garch_total/real(specification%garch_order, dp)
      select case (specification%variance_model)
      case (rugarch_model_egarch)
         parameters%omega = (1.0_dp - sum(parameters%beta))*log(variance_value)
      case (rugarch_model_igarch)
         parameters%omega = 0.01_dp*variance_value
      case (rugarch_model_figarch)
         parameters%omega = 0.03_dp*variance_value
         parameters%fractional = 0.40_dp
         parameters%alpha = 0.10_dp
         parameters%beta = 0.40_dp
      case (rugarch_model_csgarch)
         parameters%omega = 0.05_dp*variance_value
         parameters%alpha = 0.10_dp
         parameters%beta = 0.60_dp
         parameters%component_rho = 0.95_dp
         parameters%component_phi = 0.08_dp
      case (rugarch_model_realgarch)
         parameters%alpha = 0.10_dp
         parameters%beta = 0.80_dp
         parameters%omega = 0.10_dp*log(variance_value)
         parameters%measurement_intercept = 0.0_dp
         parameters%measurement_slope = 1.0_dp
         parameters%measurement_leverage1 = -0.05_dp
         parameters%measurement_leverage2 = 0.02_dp
         parameters%measurement_sigma = 0.30_dp
      case (rugarch_model_fgarch)
         parameters%omega = max(0.02_dp*variance_value, tiny(1.0_dp))
         parameters%alpha = arch_total/real(max(1, specification%arch_order), dp)
         parameters%beta = garch_total/real(max(1, specification%garch_order), dp)
         call set_fgarch_defaults(specification%fgarch_submodel, parameters)
      case default
         parameters%omega = max(variance_value* &
            (1.0_dp - arch_total - garch_total), 0.01_dp*variance_value)
      end select
      select case (specification%distribution)
      case (rugarch_distribution_student, rugarch_distribution_skew_student)
         parameters%shape = 8.0_dp
      case (rugarch_distribution_ged, rugarch_distribution_skew_ged)
         parameters%shape = 1.5_dp
      case (rugarch_distribution_johnson_su)
         parameters%shape = 1.5_dp
         parameters%skew = 0.0_dp
      case (rugarch_distribution_nig)
         parameters%shape = 3.0_dp
         parameters%skew = 0.0_dp
      case (rugarch_distribution_ghyp)
         parameters%shape = 3.0_dp
         parameters%skew = 0.0_dp
         parameters%lambda = 1.0_dp
      case (rugarch_distribution_gh_skew_student)
         parameters%shape = 8.0_dp
         parameters%skew = 0.0_dp
      case default
         parameters%shape = 0.0_dp
      end select
   end function default_parameters

   pure subroutine set_fgarch_defaults(submodel, parameters)
      !! Apply rugarch Hentschel-family starting values for a named submodel.
      integer, intent(in) :: submodel !! Hentschel fGARCH submodel code.
      type(rugarch_parameters_t), intent(inout) :: parameters !! Parameters to initialize.

      parameters%fgarch_lambda = 2.0_dp
      parameters%fgarch_delta = 0.0_dp
      parameters%asymmetry = 0.0_dp
      parameters%fgarch_shift = 0.0_dp
      select case (submodel)
      case (rugarch_fgarch_tgarch)
         parameters%fgarch_lambda = 1.0_dp
         parameters%fgarch_delta = 1.0_dp
         parameters%asymmetry = 0.05_dp
      case (rugarch_fgarch_avgarch)
         parameters%fgarch_lambda = 1.0_dp
         parameters%fgarch_delta = 1.0_dp
         parameters%asymmetry = 0.02_dp
         parameters%fgarch_shift = 0.05_dp
      case (rugarch_fgarch_ngarch)
         parameters%fgarch_delta = 0.0_dp
      case (rugarch_fgarch_nagarch)
         parameters%fgarch_delta = 2.0_dp
         parameters%fgarch_shift = 0.05_dp
      case (rugarch_fgarch_aparch)
         parameters%fgarch_lambda = 1.0_dp
         parameters%asymmetry = 0.05_dp
      case (rugarch_fgarch_allgarch)
         parameters%asymmetry = 0.05_dp
         parameters%fgarch_shift = 0.05_dp
      case (rugarch_fgarch_gjrgarch)
         parameters%fgarch_delta = 2.0_dp
         parameters%asymmetry = 0.05_dp
      case default
         parameters%fgarch_delta = 2.0_dp
      end select
   end subroutine set_fgarch_defaults

   pure subroutine apply_fgarch_submodel(submodel, parameters)
      !! Enforce the fixed Hentschel parameters defining each fGARCH submodel.
      integer, intent(in) :: submodel !! Hentschel fGARCH submodel code.
      type(rugarch_parameters_t), intent(inout) :: parameters !! Parameters to constrain.

      select case (submodel)
      case (rugarch_fgarch_garch)
         parameters%fgarch_lambda = 2.0_dp
         parameters%fgarch_delta = 2.0_dp
         parameters%asymmetry = 0.0_dp
         parameters%fgarch_shift = 0.0_dp
      case (rugarch_fgarch_tgarch)
         parameters%fgarch_lambda = 1.0_dp
         parameters%fgarch_delta = 1.0_dp
         parameters%fgarch_shift = 0.0_dp
      case (rugarch_fgarch_avgarch)
         parameters%fgarch_lambda = 1.0_dp
         parameters%fgarch_delta = 1.0_dp
      case (rugarch_fgarch_ngarch)
         parameters%fgarch_delta = 0.0_dp
         parameters%asymmetry = 0.0_dp
         parameters%fgarch_shift = 0.0_dp
      case (rugarch_fgarch_nagarch)
         parameters%fgarch_lambda = 2.0_dp
         parameters%fgarch_delta = 2.0_dp
         parameters%asymmetry = 0.0_dp
      case (rugarch_fgarch_aparch)
         parameters%fgarch_delta = 0.0_dp
         parameters%fgarch_shift = 0.0_dp
      case (rugarch_fgarch_allgarch)
         parameters%fgarch_delta = 0.0_dp
      case (rugarch_fgarch_gjrgarch)
         parameters%fgarch_lambda = 2.0_dp
         parameters%fgarch_delta = 2.0_dp
         parameters%fgarch_shift = 0.0_dp
      end select
   end subroutine apply_fgarch_submodel

   pure elemental real(dp) function fgarch_fk(submodel) result(value)
      !! Return the Hentschel switch coupling delta to the variance power.
      integer, intent(in) :: submodel !! Hentschel fGARCH submodel code.

      select case (submodel)
      case (rugarch_fgarch_ngarch, rugarch_fgarch_aparch, &
         rugarch_fgarch_allgarch)
         value = 1.0_dp
      case default
         value = 0.0_dp
      end select
   end function fgarch_fk

   pure real(dp) function rugarch_fgarch_exponent(specification, parameters) &
      result(value)
      !! Return the submodel-specific exponent of the standardized fGARCH shock.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.

      value = max(parameters%fgarch_delta + fgarch_fk( &
         specification%fgarch_submodel)*parameters%fgarch_lambda, 0.01_dp)
   end function rugarch_fgarch_exponent

   pure real(dp) function fgarch_innovation_moment(specification, parameters, &
      lag) result(moment)
      !! Numerically evaluate the standardized fGARCH shock moment.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      integer, intent(in) :: lag !! ARCH lag whose asymmetry parameters are used.
      integer, parameter :: intervals = 512
      real(dp), parameter :: lower = -12.0_dp
      real(dp), parameter :: upper = 12.0_dp
      real(dp) :: value, width, exponent
      integer :: point

      width = (upper - lower)/real(intervals, dp)
      exponent = rugarch_fgarch_exponent(specification, parameters)
      moment = 0.0_dp
      do point = 1, intervals
         value = lower + (real(point, dp) - 0.5_dp)*width
         moment = moment + (abs(value - parameters%fgarch_shift(lag)) - &
            parameters%asymmetry(lag)*(value - &
            parameters%fgarch_shift(lag)))**exponent*exp(rugarch_log_density( &
            value, specification%distribution, parameters%shape, &
            parameters%skew, parameters%lambda))
      end do
      moment = max(moment*width, tiny(1.0_dp))
   end function fgarch_innovation_moment

   pure integer function raw_parameter_count(specification) result(count)
      !! Return the unconstrained parameter-vector length.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.

      count = merge(1, 0, specification%include_mean) + &
         specification%ar_order + specification%ma_order + &
         merge(1, 0, specification%arch_in_mean) + &
         merge(1, 0, specification%fractional_mean) + &
         specification%mean_regressor_count + &
         specification%variance_regressor_count + 1
      select case (specification%variance_model)
      case (rugarch_model_sgarch, rugarch_model_igarch)
         count = count + specification%arch_order + specification%garch_order
      case (rugarch_model_gjrgarch)
         count = count + 2*specification%arch_order + &
            specification%garch_order
      case (rugarch_model_egarch)
         count = count + 2*specification%arch_order + &
            specification%garch_order
      case (rugarch_model_aparch)
         count = count + 2*specification%arch_order + &
            specification%garch_order + 1
      case (rugarch_model_figarch)
         count = count + 3
      case (rugarch_model_csgarch)
         count = count + specification%arch_order + &
            specification%garch_order + 2
      case (rugarch_model_realgarch)
         count = count + specification%arch_order + &
            specification%garch_order + 5
      case (rugarch_model_fgarch)
         count = count + 3*specification%arch_order + &
            specification%garch_order + 2
      end select
      if (distribution_has_skew(specification%distribution)) count = count + 1
      if (distribution_has_shape(specification%distribution)) count = count + 1
      if (distribution_has_lambda(specification%distribution)) count = count + 1
   end function raw_parameter_count

   pure function raw_to_parameters(specification, raw) result(parameters)
      !! Map unconstrained optimizer values to physical rugarch parameters.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      real(dp), intent(in) :: raw(:) !! Unconstrained parameter vector.
      type(rugarch_parameters_t) :: parameters
      real(dp), allocatable :: weights(:)
      integer :: index, count, p, q

      p = specification%arch_order
      q = specification%garch_order
      allocate(parameters%ar(specification%ar_order), &
         parameters%ma(specification%ma_order), parameters%alpha(p), &
         parameters%beta(q), parameters%asymmetry(p), &
         parameters%fgarch_shift(p), &
         parameters%mean_regression(specification%mean_regressor_count), &
         parameters%variance_regression( &
         specification%variance_regressor_count))
      parameters%ar = 0.0_dp
      parameters%ma = 0.0_dp
      parameters%alpha = 0.0_dp
      parameters%beta = 0.0_dp
      parameters%asymmetry = 0.0_dp
      parameters%fgarch_shift = 0.0_dp
      parameters%mean_regression = 0.0_dp
      parameters%variance_regression = 0.0_dp
      index = 1
      if (specification%include_mean) then
         parameters%mean = raw(index)
         index = index + 1
      end if
      if (specification%ar_order > 0) then
         parameters%ar = raw(index:index + specification%ar_order - 1)
         index = index + specification%ar_order
      end if
      if (specification%ma_order > 0) then
         parameters%ma = raw(index:index + specification%ma_order - 1)
         index = index + specification%ma_order
      end if
      if (specification%arch_in_mean) then
         parameters%arch_in_mean = raw(index)
         index = index + 1
      end if
      if (specification%fractional_mean) then
         parameters%mean_fractional = 0.49_dp*tanh(raw(index))
         index = index + 1
      end if
      if (specification%mean_regressor_count > 0) then
         parameters%mean_regression = raw(index:index + &
            specification%mean_regressor_count - 1)
         index = index + specification%mean_regressor_count
      end if
      if (specification%variance_regressor_count > 0) then
         parameters%variance_regression = raw(index:index + &
            specification%variance_regressor_count - 1)
         index = index + specification%variance_regressor_count
      end if
      select case (specification%variance_model)
      case (rugarch_model_sgarch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         count = p + q
         weights = constrained_weights(raw(index:index + count - 1), &
            0.995_dp, .false.)
         if (p > 0) parameters%alpha = weights(:p)
         if (q > 0) parameters%beta = weights(p + 1:)
         index = index + count
      case (rugarch_model_igarch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         count = p + q
         weights = constrained_weights(raw(index:index + count - 1), &
            1.0_dp, .true.)
         if (p > 0) parameters%alpha = weights(:p)
         if (q > 0) parameters%beta = weights(p + 1:)
         index = index + count
      case (rugarch_model_gjrgarch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         count = 2*p + q
         weights = constrained_weights(raw(index:index + count - 1), &
            0.995_dp, .false.)
         if (p > 0) then
            parameters%alpha = 2.0_dp*weights(:p)
            parameters%asymmetry = 2.0_dp*weights(p + 1:2*p) - &
               parameters%alpha
         end if
         if (q > 0) parameters%beta = weights(2*p + 1:)
         index = index + count
      case (rugarch_model_egarch)
         parameters%omega = raw(index)
         index = index + 1
         if (p > 0) then
            parameters%alpha = raw(index:index + p - 1)
            index = index + p
            parameters%asymmetry = raw(index:index + p - 1)
            index = index + p
         end if
         if (q > 0) then
            parameters%beta = constrained_weights(raw(index:index + q - 1), &
               0.995_dp, .false.)
            index = index + q
         end if
      case (rugarch_model_aparch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         count = p + q
         weights = constrained_weights(raw(index:index + count - 1), &
            0.98_dp, .false.)
         if (p > 0) parameters%alpha = weights(:p)
         if (q > 0) parameters%beta = weights(p + 1:)
         index = index + count
         if (p > 0) then
            parameters%asymmetry = tanh(raw(index:index + p - 1))
            index = index + p
         end if
         parameters%power = 0.25_dp + 3.75_dp*logistic(raw(index))
         index = index + 1
      case (rugarch_model_figarch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         parameters%fractional = logistic(raw(index))
         index = index + 1
         parameters%alpha(1) = 0.5_dp*(1.0_dp - parameters%fractional)* &
            logistic(raw(index))
         index = index + 1
         parameters%beta(1) = (parameters%fractional + &
            parameters%alpha(1))*logistic(raw(index))
         index = index + 1
      case (rugarch_model_csgarch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         weights = constrained_weights(raw(index:index + 1), 0.98_dp, .false.)
         parameters%alpha(1) = weights(1)
         parameters%beta(1) = weights(2)
         index = index + 2
         parameters%component_rho = 0.999_dp*logistic(raw(index))
         index = index + 1
         parameters%component_phi = 0.25_dp*logistic(raw(index))
         index = index + 1
      case (rugarch_model_realgarch)
         parameters%omega = raw(index)
         index = index + 1
         weights = constrained_weights(raw(index:index + 1), 0.98_dp, .false.)
         parameters%alpha(1) = weights(1)
         parameters%beta(1) = weights(2)
         index = index + 2
         parameters%measurement_intercept = raw(index)
         index = index + 1
         parameters%measurement_slope = exp(min(raw(index), 20.0_dp))
         index = index + 1
         parameters%measurement_leverage1 = raw(index)
         index = index + 1
         parameters%measurement_leverage2 = raw(index)
         index = index + 1
         parameters%measurement_sigma = exp(min(raw(index), 20.0_dp))
         index = index + 1
      case (rugarch_model_fgarch)
         parameters%omega = exp(min(raw(index), 50.0_dp))
         index = index + 1
         count = p + q
         weights = constrained_weights(raw(index:index + count - 1), &
            0.98_dp, .false.)
         if (p > 0) parameters%alpha = weights(:p)
         if (q > 0) parameters%beta = weights(p + 1:)
         index = index + count
         parameters%asymmetry = tanh(raw(index:index + p - 1))
         index = index + p
         parameters%fgarch_shift = 10.0_dp*tanh(raw(index:index + p - 1))
         index = index + p
         parameters%fgarch_delta = 4.0_dp*logistic(raw(index))
         index = index + 1
         parameters%fgarch_lambda = 0.01_dp + &
            3.99_dp*logistic(raw(index))
         index = index + 1
         call apply_fgarch_submodel(specification%fgarch_submodel, parameters)
      end select
      if (distribution_has_skew(specification%distribution)) then
         if (specification%distribution == rugarch_distribution_johnson_su .or. &
            specification%distribution == &
            rugarch_distribution_gh_skew_student) then
            parameters%skew = raw(index)
         else if (specification%distribution == rugarch_distribution_nig .or. &
            specification%distribution == rugarch_distribution_ghyp) then
            parameters%skew = 0.999_dp*tanh(raw(index))
         else
            parameters%skew = exp(min(raw(index), 20.0_dp))
         end if
         index = index + 1
      end if
      if (specification%distribution == rugarch_distribution_student .or. &
         specification%distribution == rugarch_distribution_skew_student) then
         parameters%shape = 2.01_dp + exp(min(raw(index), 20.0_dp))
      else if (specification%distribution == rugarch_distribution_ged .or. &
         specification%distribution == rugarch_distribution_skew_ged) then
         parameters%shape = 0.20_dp + exp(min(raw(index), 20.0_dp))
      else if (specification%distribution == rugarch_distribution_johnson_su) then
         parameters%shape = 0.20_dp + 19.80_dp*logistic(raw(index))
      else if (specification%distribution == rugarch_distribution_nig .or. &
         specification%distribution == rugarch_distribution_ghyp) then
         parameters%shape = 0.10_dp + exp(min(raw(index), 20.0_dp))
      else if (specification%distribution == &
         rugarch_distribution_gh_skew_student) then
         parameters%shape = 4.01_dp + exp(min(raw(index), 20.0_dp))
      end if
      if (distribution_has_shape(specification%distribution)) index = index + 1
      if (distribution_has_lambda(specification%distribution)) &
         parameters%lambda = raw(index)
   end function raw_to_parameters

   pure function parameters_to_raw(specification, parameters) result(raw)
      !! Map physical starting values to unconstrained optimizer coordinates.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp), allocatable :: raw(:)
      real(dp), allocatable :: weights(:)
      integer :: index, count, p, q

      p = specification%arch_order
      q = specification%garch_order
      allocate(raw(raw_parameter_count(specification)))
      raw = 0.0_dp
      index = 1
      if (specification%include_mean) then
         raw(index) = parameters%mean
         index = index + 1
      end if
      if (specification%ar_order > 0) then
         raw(index:index + specification%ar_order - 1) = parameters%ar
         index = index + specification%ar_order
      end if
      if (specification%ma_order > 0) then
         raw(index:index + specification%ma_order - 1) = parameters%ma
         index = index + specification%ma_order
      end if
      if (specification%arch_in_mean) then
         raw(index) = parameters%arch_in_mean
         index = index + 1
      end if
      if (specification%fractional_mean) then
         raw(index) = atanh(max(-0.999_dp, min(0.999_dp, &
            parameters%mean_fractional/0.49_dp)))
         index = index + 1
      end if
      if (specification%mean_regressor_count > 0) then
         raw(index:index + specification%mean_regressor_count - 1) = &
            parameters%mean_regression
         index = index + specification%mean_regressor_count
      end if
      if (specification%variance_regressor_count > 0) then
         raw(index:index + specification%variance_regressor_count - 1) = &
            parameters%variance_regression
         index = index + specification%variance_regressor_count
      end if
      select case (specification%variance_model)
      case (rugarch_model_sgarch, rugarch_model_igarch)
         raw(index) = log(max(parameters%omega, tiny(1.0_dp)))
         index = index + 1
         weights = [parameters%alpha, parameters%beta]
         count = size(weights)
         if (specification%variance_model == rugarch_model_igarch) then
            raw(index:index + count - 1) = log(max(weights, 1.0e-8_dp))
         else
            raw(index:index + count - 1) = inverse_constrained_weights( &
               weights, 0.995_dp)
         end if
         index = index + count
      case (rugarch_model_gjrgarch)
         raw(index) = log(max(parameters%omega, tiny(1.0_dp)))
         index = index + 1
         weights = [0.5_dp*parameters%alpha, &
            0.5_dp*(parameters%alpha + parameters%asymmetry), parameters%beta]
         count = size(weights)
         raw(index:index + count - 1) = inverse_constrained_weights( &
            weights, 0.995_dp)
         index = index + count
      case (rugarch_model_egarch)
         raw(index) = parameters%omega
         index = index + 1
         if (p > 0) then
            raw(index:index + p - 1) = parameters%alpha
            index = index + p
            raw(index:index + p - 1) = parameters%asymmetry
            index = index + p
         end if
         if (q > 0) then
            raw(index:index + q - 1) = inverse_constrained_weights( &
               parameters%beta, 0.995_dp)
            index = index + q
         end if
      case (rugarch_model_aparch)
         raw(index) = log(max(parameters%omega, tiny(1.0_dp)))
         index = index + 1
         weights = [parameters%alpha, parameters%beta]
         count = size(weights)
         raw(index:index + count - 1) = inverse_constrained_weights( &
            weights, 0.98_dp)
         index = index + count
         if (p > 0) then
            raw(index:index + p - 1) = atanh(max(-0.999_dp, &
               min(0.999_dp, parameters%asymmetry)))
            index = index + p
         end if
         raw(index) = logit((parameters%power - 0.25_dp)/3.75_dp)
         index = index + 1
      case (rugarch_model_figarch)
         raw(index) = log(max(parameters%omega, tiny(1.0_dp)))
         index = index + 1
         raw(index) = logit(parameters%fractional)
         index = index + 1
         raw(index) = logit(parameters%alpha(1)/max(0.5_dp* &
            (1.0_dp - parameters%fractional), 1.0e-8_dp))
         index = index + 1
         raw(index) = logit(parameters%beta(1)/max(parameters%fractional + &
            parameters%alpha(1), 1.0e-8_dp))
         index = index + 1
      case (rugarch_model_csgarch)
         raw(index) = log(max(parameters%omega, tiny(1.0_dp)))
         index = index + 1
         weights = [parameters%alpha, parameters%beta]
         raw(index:index + 1) = inverse_constrained_weights(weights, 0.98_dp)
         index = index + 2
         raw(index) = logit(parameters%component_rho/0.999_dp)
         index = index + 1
         raw(index) = logit(parameters%component_phi/0.25_dp)
         index = index + 1
      case (rugarch_model_realgarch)
         raw(index) = parameters%omega
         index = index + 1
         weights = [parameters%alpha, parameters%beta]
         raw(index:index + 1) = inverse_constrained_weights(weights, 0.98_dp)
         index = index + 2
         raw(index) = parameters%measurement_intercept
         index = index + 1
         raw(index) = log(max(parameters%measurement_slope, 1.0e-8_dp))
         index = index + 1
         raw(index) = parameters%measurement_leverage1
         index = index + 1
         raw(index) = parameters%measurement_leverage2
         index = index + 1
         raw(index) = log(max(parameters%measurement_sigma, 1.0e-8_dp))
         index = index + 1
      case (rugarch_model_fgarch)
         raw(index) = log(max(parameters%omega, tiny(1.0_dp)))
         index = index + 1
         weights = [parameters%alpha, parameters%beta]
         count = size(weights)
         raw(index:index + count - 1) = inverse_constrained_weights( &
            weights, 0.98_dp)
         index = index + count
         raw(index:index + p - 1) = atanh(max(-0.999_dp, &
            min(0.999_dp, parameters%asymmetry)))
         index = index + p
         raw(index:index + p - 1) = atanh(max(-0.999_dp, &
            min(0.999_dp, parameters%fgarch_shift/10.0_dp)))
         index = index + p
         raw(index) = logit(parameters%fgarch_delta/4.0_dp)
         index = index + 1
         raw(index) = logit((parameters%fgarch_lambda - 0.01_dp)/3.99_dp)
         index = index + 1
      end select
      if (distribution_has_skew(specification%distribution)) then
         if (specification%distribution == rugarch_distribution_johnson_su .or. &
            specification%distribution == &
            rugarch_distribution_gh_skew_student) then
            raw(index) = parameters%skew
         else if (specification%distribution == rugarch_distribution_nig .or. &
            specification%distribution == rugarch_distribution_ghyp) then
            raw(index) = atanh(max(-0.999_dp, min(0.999_dp, &
               parameters%skew/0.999_dp)))
         else
            raw(index) = log(max(parameters%skew, 1.0e-8_dp))
         end if
         index = index + 1
      end if
      if (specification%distribution == rugarch_distribution_student .or. &
         specification%distribution == rugarch_distribution_skew_student) then
         raw(index) = log(max(parameters%shape - 2.01_dp, 1.0e-6_dp))
      else if (specification%distribution == rugarch_distribution_ged .or. &
         specification%distribution == rugarch_distribution_skew_ged) then
         raw(index) = log(max(parameters%shape - 0.20_dp, 1.0e-6_dp))
      else if (specification%distribution == rugarch_distribution_johnson_su) then
         raw(index) = logit((parameters%shape - 0.20_dp)/19.80_dp)
      else if (specification%distribution == rugarch_distribution_nig .or. &
         specification%distribution == rugarch_distribution_ghyp) then
         raw(index) = log(max(parameters%shape - 0.10_dp, 1.0e-6_dp))
      else if (specification%distribution == &
         rugarch_distribution_gh_skew_student) then
         raw(index) = log(max(parameters%shape - 4.01_dp, 1.0e-6_dp))
      end if
      if (distribution_has_shape(specification%distribution)) index = index + 1
      if (distribution_has_lambda(specification%distribution)) &
         raw(index) = parameters%lambda
   end function parameters_to_raw

   pure function constrained_weights(raw, scale, exact_sum) result(weights)
      !! Map unconstrained values to positive coefficients with bounded sum.
      real(dp), intent(in) :: raw(:) !! Unconstrained coefficient values.
      real(dp), intent(in) :: scale !! Target upper bound or exact sum.
      logical, intent(in) :: exact_sum !! Normalize coefficients to the exact scale.
      real(dp) :: weights(size(raw))
      real(dp) :: exponential(size(raw)), maximum, denominator

      if (size(raw) == 0) return
      maximum = max(0.0_dp, maxval(raw))
      exponential = exp(max(-50.0_dp, min(50.0_dp, raw - maximum)))
      if (exact_sum) then
         denominator = sum(exponential)
      else
         denominator = exp(-maximum) + sum(exponential)
      end if
      weights = scale*exponential/max(denominator, tiny(1.0_dp))
   end function constrained_weights

   pure function inverse_constrained_weights(weights, scale) result(raw)
      !! Return softmax coordinates for positive coefficients below a sum bound.
      real(dp), intent(in) :: weights(:) !! Positive constrained coefficients.
      real(dp), intent(in) :: scale !! Sum upper bound.
      real(dp) :: raw(size(weights))
      real(dp) :: slack

      slack = max(scale - sum(weights), 1.0e-8_dp)
      raw = log(max(weights, 1.0e-8_dp)/slack)
   end function inverse_constrained_weights

   pure elemental real(dp) function logistic(value) result(transformed)
      !! Apply a numerically stable logistic transform.
      real(dp), intent(in) :: value !! Unconstrained value.

      if (value >= 0.0_dp) then
         transformed = 1.0_dp/(1.0_dp + exp(-min(value, 50.0_dp)))
      else
         transformed = exp(max(value, -50.0_dp))/ &
            (1.0_dp + exp(max(value, -50.0_dp)))
      end if
   end function logistic

   pure elemental real(dp) function logit(probability) result(value)
      !! Apply a clipped inverse logistic transform.
      real(dp), intent(in) :: probability !! Probability value.
      real(dp) :: clipped

      clipped = max(1.0e-8_dp, min(1.0_dp - 1.0e-8_dp, probability))
      value = log(clipped/(1.0_dp - clipped))
   end function logit

   pure elemental logical function distribution_has_skew(distribution) &
      result(has_skew)
      !! Report whether an innovation distribution has a skew parameter.
      integer, intent(in) :: distribution !! Distribution code.

      has_skew = shared_distribution_has_skew(distribution)
   end function distribution_has_skew

   pure elemental logical function distribution_has_shape(distribution) &
      result(has_shape)
      !! Report whether an innovation distribution has a shape parameter.
      integer, intent(in) :: distribution !! Distribution code.

      has_shape = shared_distribution_has_shape(distribution)
   end function distribution_has_shape

   pure real(dp) function arma_constraint_penalty(parameters) result(penalty)
      !! Penalize simple AR and MA coefficient sums outside stable bounds.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp) :: excess

      penalty = 0.0_dp
      if (size(parameters%ar) > 0) then
         excess = max(0.0_dp, sum(abs(parameters%ar)) - 0.999_dp)
         penalty = penalty + 1.0e6_dp*excess**2
      end if
      if (size(parameters%ma) > 0) then
         excess = max(0.0_dp, sum(abs(parameters%ma)) - 0.999_dp)
         penalty = penalty + 1.0e6_dp*excess**2
      end if
   end function arma_constraint_penalty

   pure function rugarch_coefficients(specification, parameters) result(values)
      !! Pack physical coefficients in rugarch mean-then-variance order.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp), allocatable :: values(:)

      allocate(values(0))
      if (specification%include_mean) values = [values, parameters%mean]
      values = [values, parameters%ar, parameters%ma]
      if (specification%arch_in_mean) &
         values = [values, parameters%arch_in_mean]
      if (specification%fractional_mean) &
         values = [values, parameters%mean_fractional]
      if (specification%mean_regressor_count > 0) &
         values = [values, parameters%mean_regression]
      if (specification%variance_regressor_count > 0) &
         values = [values, parameters%variance_regression]
      values = [values, parameters%omega, parameters%alpha, parameters%beta]
      if (specification%variance_model == rugarch_model_egarch .or. &
         specification%variance_model == rugarch_model_gjrgarch .or. &
         specification%variance_model == rugarch_model_aparch) &
         values = [values, parameters%asymmetry]
      if (specification%variance_model == rugarch_model_aparch) &
         values = [values, parameters%power]
      if (specification%variance_model == rugarch_model_figarch) &
         values = [values, parameters%fractional]
      if (specification%variance_model == rugarch_model_csgarch) &
         values = [values, parameters%component_rho, parameters%component_phi]
      if (specification%variance_model == rugarch_model_realgarch) &
         values = [values, parameters%measurement_intercept, &
         parameters%measurement_slope, parameters%measurement_leverage1, &
         parameters%measurement_leverage2, parameters%measurement_sigma]
      if (specification%variance_model == rugarch_model_fgarch) &
         values = [values, parameters%asymmetry, parameters%fgarch_shift, &
         parameters%fgarch_delta, parameters%fgarch_lambda]
      if (distribution_has_skew(specification%distribution)) &
         values = [values, parameters%skew]
      if (distribution_has_shape(specification%distribution)) &
         values = [values, parameters%shape]
      if (distribution_has_lambda(specification%distribution)) &
         values = [values, parameters%lambda]
   end function rugarch_coefficients

   pure function physical_parameter_jacobian(specification, raw) &
      result(jacobian)
      !! Numerically differentiate packed physical coefficients by raw values.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      real(dp), intent(in) :: raw(:) !! Unconstrained optimizer parameters.
      real(dp), allocatable :: jacobian(:, :)
      real(dp) :: shifted(size(raw)), step
      real(dp), allocatable :: upper(:), lower(:)
      integer :: parameter

      upper = rugarch_coefficients(specification, &
         raw_to_parameters(specification, raw))
      allocate(jacobian(size(upper), size(raw)))
      do parameter = 1, size(raw)
         step = epsilon(1.0_dp)**(1.0_dp/3.0_dp)* &
            max(1.0_dp, abs(raw(parameter)))
         shifted = raw
         shifted(parameter) = raw(parameter) + step
         upper = rugarch_coefficients(specification, &
            raw_to_parameters(specification, shifted))
         shifted(parameter) = raw(parameter) - step
         lower = rugarch_coefficients(specification, &
            raw_to_parameters(specification, shifted))
         jacobian(:, parameter) = (upper - lower)/(2.0_dp*step)
      end do
   end function physical_parameter_jacobian

   pure function figarch_weights(phi, fractional, beta, truncation) &
      result(weights)
      !! Construct Chung-style FIGARCH fractional ARCH weights.
      real(dp), intent(in) :: phi !! Short-memory fractional AR coefficient.
      real(dp), intent(in) :: fractional !! Fractional differencing parameter.
      real(dp), intent(in) :: beta !! GARCH coefficient.
      integer, intent(in) :: truncation !! Positive truncation lag.
      real(dp) :: weights(truncation)
      real(dp) :: previous_delta, current_delta
      integer :: lag

      if (truncation < 1) return
      weights(1) = phi - beta + fractional
      previous_delta = fractional
      do lag = 2, truncation
         current_delta = (real(lag - 1, dp) - fractional)/ &
            real(lag, dp)*previous_delta
         weights(lag) = beta*weights(lag - 1) + current_delta - &
            phi*previous_delta
         previous_delta = current_delta
      end do
   end function figarch_weights

   pure real(dp) function rugarch_persistence(specification, parameters) &
      result(value)
      !! Return the model-specific conditional-scale persistence.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      integer :: lag

      select case (specification%variance_model)
      case (rugarch_model_sgarch)
         value = sum(parameters%alpha) + sum(parameters%beta)
      case (rugarch_model_igarch)
         value = 1.0_dp
      case (rugarch_model_gjrgarch)
         value = sum(parameters%alpha + 0.5_dp*parameters%asymmetry) + &
            sum(parameters%beta)
      case (rugarch_model_egarch)
         value = sum(parameters%beta)
      case (rugarch_model_aparch)
         value = sum(parameters%beta)
         do lag = 1, size(parameters%alpha)
             value = value + parameters%alpha(lag)*innovation_aparch_moment( &
                specification%distribution, parameters%shape, &
                parameters%skew, parameters%lambda, parameters%power, &
                parameters%asymmetry(lag))
         end do
      case (rugarch_model_figarch)
         value = 1.0_dp
      case (rugarch_model_csgarch)
         value = parameters%component_rho
      case (rugarch_model_realgarch)
         value = parameters%beta(1) + &
            parameters%alpha(1)*parameters%measurement_slope
      case (rugarch_model_fgarch)
         value = sum(parameters%beta)
         do lag = 1, size(parameters%alpha)
            value = value + parameters%alpha(lag)* &
               fgarch_innovation_moment(specification, parameters, lag)
         end do
      case default
         value = huge(1.0_dp)
      end select
   end function rugarch_persistence

   pure real(dp) function rugarch_unconditional_variance(specification, &
      parameters) result(value)
      !! Return the model's unconditional variance when it is finite.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp) :: persistence

      persistence = rugarch_persistence(specification, parameters)
      if (specification%variance_model == rugarch_model_igarch .or. &
         persistence >= 1.0_dp) then
         value = huge(1.0_dp)
      else if (specification%variance_model == rugarch_model_egarch) then
         value = exp(parameters%omega/max(1.0_dp - persistence, &
            tiny(1.0_dp)))
      else if (specification%variance_model == rugarch_model_aparch) then
         value = (parameters%omega/max(1.0_dp - persistence, &
            tiny(1.0_dp)))**(2.0_dp/parameters%power)
      else if (specification%variance_model == rugarch_model_fgarch) then
         value = (parameters%omega/max(1.0_dp - persistence, &
            tiny(1.0_dp)))**(2.0_dp/parameters%fgarch_lambda)
      else if (specification%variance_model == rugarch_model_realgarch) then
         value = exp((parameters%omega + parameters%alpha(1)* &
            parameters%measurement_intercept)/max(1.0_dp - persistence, &
            tiny(1.0_dp)))
      else
         value = parameters%omega/max(1.0_dp - persistence, tiny(1.0_dp))
      end if
   end function rugarch_unconditional_variance

   pure real(dp) function rugarch_half_life(specification, parameters) &
      result(value)
      !! Return the shock half-life implied by conditional-scale persistence.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      real(dp) :: persistence

      persistence = rugarch_persistence(specification, parameters)
      if (persistence <= 0.0_dp) then
         value = 0.0_dp
      else if (persistence >= 1.0_dp) then
         value = huge(1.0_dp)
      else
         value = log(0.5_dp)/log(persistence)
      end if
   end function rugarch_half_life

   pure function rugarch_forecast(fit, horizon, mean_regressors, &
      variance_regressors) result(out)
      !! Translate ugarchforecast analytic recursive mean and variance forecasts.
      type(rugarch_fit_t), intent(in) :: fit !! Fitted rugarch model.
      integer, intent(in) :: horizon !! Positive forecast horizon.
      real(dp), intent(in), optional :: mean_regressors(:, :) !! Future mean regressors.
      real(dp), intent(in), optional :: variance_regressors(:, :) !! Future variance regressors.
      type(rugarch_forecast_t) :: out
      real(dp), allocatable :: observed(:), residuals(:), variances(:)
      real(dp), allocatable :: standardized(:), powered_scale(:)
      real(dp), allocatable :: fractional_weight(:)
      real(dp), allocatable :: mean_fractional_weight(:)
      real(dp), allocatable :: component_path(:), realized_path(:)
      real(dp) :: location, variance_value, scale_value, expected_absolute
      real(dp) :: component_value, previous_shock
      integer :: observations, total, time, lag, step

      if (fit%info /= 0 .or. horizon < 1 .or. &
         .not. allocated(fit%filtered%residuals)) then
         out%info = 1
         return
      end if
      if (.not. valid_regressors(horizon, fit%specification, &
         mean_regressors, variance_regressors)) then
         out%info = 1
         return
      end if
      observations = size(fit%filtered%residuals)
      total = observations + horizon
      allocate(observed(total), residuals(total), variances(total), &
         standardized(total), powered_scale(total), component_path(total), &
         realized_path(total))
      observed(:observations) = fit%filtered%conditional_mean + &
         fit%filtered%residuals
      residuals(:observations) = fit%filtered%residuals
      variances(:observations) = fit%filtered%conditional_variance
      standardized(:observations) = fit%filtered%standardized_residuals
      powered_scale(:observations) = variances(:observations)** &
         (0.5_dp*fit%parameters%power)
      component_path = 0.0_dp
      realized_path = 0.0_dp
      component_path(:observations) = fit%filtered%component_variance
      realized_path(:observations) = fit%filtered%realized_variance
      if (fit%specification%variance_model == rugarch_model_figarch) &
         fractional_weight = figarch_weights(fit%parameters%alpha(1), &
            fit%parameters%fractional, fit%parameters%beta(1), &
            fit%specification%truncation_lag)
      if (fit%specification%fractional_mean) then
         allocate(mean_fractional_weight(min(fit%specification%truncation_lag, &
            total - 1)))
         mean_fractional_weight = fractional_mean_weights( &
            fit%parameters%mean_fractional, size(mean_fractional_weight))
      end if
      residuals(observations + 1:) = 0.0_dp
      standardized(observations + 1:) = 0.0_dp
      allocate(out%mean(horizon), out%variance(horizon), out%sigma(horizon), &
         out%realized_variance(horizon))
      out%realized_variance = 0.0_dp
      expected_absolute = 0.0_dp
      if (fit%specification%variance_model == rugarch_model_egarch) &
         expected_absolute = innovation_absolute_moment( &
         fit%specification%distribution, fit%parameters%shape, &
         fit%parameters%skew, fit%parameters%lambda, 1.0_dp)
      do step = 1, horizon
         time = observations + step
         select case (fit%specification%variance_model)
         case (rugarch_model_sgarch, rugarch_model_igarch, &
            rugarch_model_gjrgarch)
            variance_value = fit%parameters%omega
            do lag = 1, fit%specification%arch_order
               if (time - lag <= observations) then
                  variance_value = variance_value + fit%parameters%alpha(lag)* &
                     residuals(time - lag)**2
                  if (fit%specification%variance_model == &
                     rugarch_model_gjrgarch .and. residuals(time - lag) < 0.0_dp) &
                     variance_value = variance_value + &
                        fit%parameters%asymmetry(lag)*residuals(time - lag)**2
               else
                  variance_value = variance_value + &
                     (fit%parameters%alpha(lag) + merge(0.5_dp* &
                     fit%parameters%asymmetry(lag), 0.0_dp, &
                     fit%specification%variance_model == &
                     rugarch_model_gjrgarch))*variances(time - lag)
               end if
            end do
            do lag = 1, fit%specification%garch_order
               variance_value = variance_value + fit%parameters%beta(lag)* &
                  variances(time - lag)
            end do
         case (rugarch_model_egarch)
            scale_value = fit%parameters%omega
            do lag = 1, fit%specification%arch_order
               if (time - lag <= observations) scale_value = scale_value + &
                  fit%parameters%alpha(lag)*(abs(standardized(time - lag)) - &
                  expected_absolute) + fit%parameters%asymmetry(lag)* &
                  standardized(time - lag)
            end do
            do lag = 1, fit%specification%garch_order
               scale_value = scale_value + fit%parameters%beta(lag)* &
                  log(max(variances(time - lag), tiny(1.0_dp)))
            end do
            variance_value = exp(max(-50.0_dp, min(50.0_dp, scale_value)))
         case (rugarch_model_aparch)
            scale_value = fit%parameters%omega
            do lag = 1, fit%specification%arch_order
               if (time - lag <= observations) then
                  scale_value = scale_value + fit%parameters%alpha(lag)* &
                     (abs(residuals(time - lag)) - &
                     fit%parameters%asymmetry(lag)*residuals(time - lag))** &
                     fit%parameters%power
               else
                  scale_value = scale_value + fit%parameters%alpha(lag)* &
                      innovation_aparch_moment(fit%specification%distribution, &
                      fit%parameters%shape, fit%parameters%skew, &
                      fit%parameters%lambda, fit%parameters%power, &
                     fit%parameters%asymmetry(lag))*powered_scale(time - lag)
               end if
            end do
            do lag = 1, fit%specification%garch_order
               scale_value = scale_value + fit%parameters%beta(lag)* &
                  powered_scale(time - lag)
            end do
            powered_scale(time) = max(scale_value, tiny(1.0_dp))
            variance_value = powered_scale(time)**(2.0_dp/fit%parameters%power)
         case (rugarch_model_figarch)
            variance_value = fit%parameters%omega/max(1.0_dp - &
               fit%parameters%beta(1), 1.0e-8_dp)
            if (time <= size(fractional_weight)) variance_value = &
               variance_value + sum(fractional_weight(time:))*variances(1)
            do lag = 1, min(time - 1, size(fractional_weight))
               if (time - lag <= observations) then
                  variance_value = variance_value + fractional_weight(lag)* &
                     residuals(time - lag)**2
               else
                  variance_value = variance_value + fractional_weight(lag)* &
                     variances(time - lag)
                end if
             end do
         case (rugarch_model_csgarch)
            previous_shock = merge(residuals(time - 1)**2, &
               variances(time - 1), step == 1)
            component_value = fit%parameters%omega + &
               fit%parameters%component_rho*component_path(time - 1) + &
               fit%parameters%component_phi*(previous_shock - &
               variances(time - 1))
            component_path(time) = max(component_value, tiny(1.0_dp))
            variance_value = component_path(time) + fit%parameters%alpha(1)* &
               (previous_shock - component_path(time - 1)) + &
               fit%parameters%beta(1)*(variances(time - 1) - &
               component_path(time - 1))
         case (rugarch_model_realgarch)
            scale_value = fit%parameters%omega + fit%parameters%alpha(1)* &
               log(max(realized_path(time - 1), tiny(1.0_dp))) + &
               fit%parameters%beta(1)*log(max(variances(time - 1), &
               tiny(1.0_dp)))
            variance_value = exp(max(-50.0_dp, min(50.0_dp, scale_value)))
            realized_path(time) = exp(fit%parameters%measurement_intercept + &
               fit%parameters%measurement_slope*log(variance_value))
         case (rugarch_model_fgarch)
            scale_value = fit%parameters%omega
            do lag = 1, fit%specification%arch_order
               scale_value = scale_value + fit%parameters%alpha(lag)* &
                  fgarch_innovation_moment(fit%specification, &
                  fit%parameters, lag)* &
                  variances(time - lag)** &
                  (0.5_dp*fit%parameters%fgarch_lambda)
            end do
            do lag = 1, fit%specification%garch_order
               scale_value = scale_value + fit%parameters%beta(lag)* &
                  variances(time - lag)** &
                  (0.5_dp*fit%parameters%fgarch_lambda)
            end do
            variance_value = max(scale_value, tiny(1.0_dp))** &
               (2.0_dp/fit%parameters%fgarch_lambda)
          end select
         if (fit%specification%variance_regressor_count > 0) then
            if (fit%specification%variance_model == rugarch_model_egarch .or. &
               fit%specification%variance_model == rugarch_model_realgarch) then
               variance_value = variance_value*exp(dot_product( &
                  variance_regressors(step, :), &
                  fit%parameters%variance_regression))
            else
               variance_value = variance_value + dot_product( &
                  variance_regressors(step, :), &
                  fit%parameters%variance_regression)
            end if
         end if
         variances(time) = max(variance_value, tiny(1.0_dp))
         location = merge(fit%parameters%mean, 0.0_dp, &
            fit%specification%include_mean)
         do lag = 1, fit%specification%ar_order
            location = location + fit%parameters%ar(lag)* &
               (observed(time - lag) - merge(fit%parameters%mean, 0.0_dp, &
               fit%specification%include_mean))
         end do
         do lag = 1, fit%specification%ma_order
            location = location + fit%parameters%ma(lag)*residuals(time - lag)
         end do
         if (fit%specification%fractional_mean) then
            do lag = 1, min(time - 1, size(mean_fractional_weight))
               location = location - mean_fractional_weight(lag)* &
                  (observed(time - lag) - merge(fit%parameters%mean, &
                  0.0_dp, fit%specification%include_mean))
            end do
         end if
         if (fit%specification%mean_regressor_count > 0) location = location + &
            dot_product(mean_regressors(step, :), &
            fit%parameters%mean_regression)
         if (fit%specification%arch_in_mean) location = location + &
            fit%parameters%arch_in_mean*sqrt(variances(time))** &
            fit%specification%arch_in_mean_power
         observed(time) = location
         out%mean(step) = location
         out%variance(step) = variances(time)
         out%sigma(step) = sqrt(variances(time))
         out%realized_variance(step) = realized_path(time)
      end do
   end function rugarch_forecast

   function rugarch_simulate(specification, parameters, observations, burnin, &
      seed) result(out)
      !! Translate ugarchsim using the library's shared random-number stream.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t), intent(in) :: parameters !! Physical parameters.
      integer, intent(in) :: observations !! Number of retained observations.
      integer, intent(in), optional :: burnin !! Number of discarded initial observations.
      integer, intent(in), optional :: seed !! Shared random-number seed.
      type(rugarch_simulation_t) :: out
      real(dp), allocatable :: series(:), means(:), variances(:), residuals(:)
      real(dp), allocatable :: standardized(:), powered_scale(:)
      real(dp), allocatable :: fractional_weight(:)
      real(dp), allocatable :: mean_fractional_weight(:)
      real(dp), allocatable :: component_path(:), realized_path(:)
      real(dp), allocatable :: measurement_errors(:)
      real(dp) :: initial_variance, location, variance_value, scale_value
      real(dp) :: expected_absolute, previous_shock, component_value
      real(dp) :: previous_component
      integer :: discarded, total, time, lag, first

      discarded = 500
      if (present(burnin)) discarded = burnin
      out%specification = specification
      out%parameters = parameters
      out%burnin = discarded
      if (.not. valid_specification(specification) .or. &
         .not. valid_parameters(specification, parameters) .or. &
         observations < 1 .or. discarded < 0) then
         out%info = 1
         return
      end if
      if (present(seed)) call set_random_seed(seed)
      total = observations + discarded
      allocate(series(total), means(total), variances(total), residuals(total), &
         standardized(total), powered_scale(total), component_path(total), &
         realized_path(total), measurement_errors(total))
      series = 0.0_dp
      means = 0.0_dp
      residuals = 0.0_dp
      standardized = 0.0_dp
      initial_variance = rugarch_unconditional_variance(specification, parameters)
      if (.not. ieee_is_finite(initial_variance) .or. &
         initial_variance >= 0.5_dp*huge(1.0_dp)) initial_variance = &
         max(parameters%omega, 1.0_dp)
      variances = initial_variance
      component_path = max(parameters%omega/max(1.0_dp - &
         parameters%component_rho, 1.0e-8_dp), tiny(1.0_dp))
      realized_path = initial_variance
      measurement_errors = 0.0_dp
      powered_scale = initial_variance**(0.5_dp*parameters%power)
      if (specification%variance_model == rugarch_model_figarch) &
         fractional_weight = figarch_weights(parameters%alpha(1), &
            parameters%fractional, parameters%beta(1), &
            specification%truncation_lag)
      if (specification%fractional_mean) then
         allocate(mean_fractional_weight(min(specification%truncation_lag, &
            total - 1)))
         mean_fractional_weight = fractional_mean_weights( &
            parameters%mean_fractional, size(mean_fractional_weight))
      end if
      expected_absolute = 0.0_dp
      if (specification%variance_model == rugarch_model_egarch) &
         expected_absolute = innovation_absolute_moment( &
         specification%distribution, parameters%shape, parameters%skew, &
         parameters%lambda, 1.0_dp)
      do time = 1, total
         select case (specification%variance_model)
         case (rugarch_model_sgarch, rugarch_model_igarch, &
            rugarch_model_gjrgarch)
            variance_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  variance_value = variance_value + parameters%alpha(lag)* &
                     residuals(time - lag)**2
                  if (specification%variance_model == &
                     rugarch_model_gjrgarch .and. residuals(time - lag) < 0.0_dp) &
                     variance_value = variance_value + &
                        parameters%asymmetry(lag)*residuals(time - lag)**2
               else
                  variance_value = variance_value + &
                     (parameters%alpha(lag) + merge(0.5_dp* &
                     parameters%asymmetry(lag), 0.0_dp, &
                     specification%variance_model == rugarch_model_gjrgarch))* &
                     initial_variance
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  variance_value = variance_value + parameters%beta(lag)* &
                     variances(time - lag)
               else
                  variance_value = variance_value + parameters%beta(lag)* &
                     initial_variance
               end if
            end do
         case (rugarch_model_egarch)
            scale_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) scale_value = scale_value + &
                  parameters%alpha(lag)*(abs(standardized(time - lag)) - &
                  expected_absolute) + parameters%asymmetry(lag)* &
                  standardized(time - lag)
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%beta(lag)* &
                     log(max(variances(time - lag), tiny(1.0_dp)))
               else
                  scale_value = scale_value + parameters%beta(lag)* &
                     log(initial_variance)
               end if
            end do
            variance_value = exp(max(-50.0_dp, min(50.0_dp, scale_value)))
         case (rugarch_model_aparch)
            scale_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%alpha(lag)* &
                     (abs(residuals(time - lag)) - parameters%asymmetry(lag)* &
                     residuals(time - lag))**parameters%power
               else
                  scale_value = scale_value + parameters%alpha(lag)* &
                      innovation_aparch_moment(specification%distribution, &
                      parameters%shape, parameters%skew, parameters%lambda, &
                      parameters%power, &
                     parameters%asymmetry(lag))*powered_scale(1)
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%beta(lag)* &
                     powered_scale(time - lag)
               else
                  scale_value = scale_value + parameters%beta(lag)* &
                     powered_scale(1)
               end if
            end do
            powered_scale(time) = max(scale_value, tiny(1.0_dp))
            variance_value = powered_scale(time)**(2.0_dp/parameters%power)
         case (rugarch_model_figarch)
            variance_value = parameters%omega/max(1.0_dp - &
               parameters%beta(1), 1.0e-8_dp)
            if (time <= size(fractional_weight)) variance_value = &
               variance_value + sum(fractional_weight(time:))*initial_variance
            do lag = 1, min(time - 1, size(fractional_weight))
                variance_value = variance_value + fractional_weight(lag)* &
                   residuals(time - lag)**2
             end do
         case (rugarch_model_csgarch)
            if (time > 1) then
               previous_shock = residuals(time - 1)**2
               previous_component = component_path(time - 1)
               scale_value = variances(time - 1)
            else
               previous_shock = initial_variance
               previous_component = component_path(1)
               scale_value = initial_variance
            end if
            component_value = parameters%omega + &
               parameters%component_rho*previous_component + &
               parameters%component_phi*(previous_shock - scale_value)
            component_path(time) = max(component_value, tiny(1.0_dp))
            variance_value = component_path(time) + parameters%alpha(1)* &
               (previous_shock - previous_component) + parameters%beta(1)* &
               (scale_value - previous_component)
         case (rugarch_model_realgarch)
            if (time == 1) then
               scale_value = log(initial_variance)
            else
               scale_value = parameters%omega + parameters%alpha(1)* &
                  log(realized_path(time - 1)) + parameters%beta(1)* &
                  log(variances(time - 1))
            end if
            variance_value = exp(max(-50.0_dp, min(50.0_dp, scale_value)))
         case (rugarch_model_fgarch)
            scale_value = parameters%omega
            do lag = 1, specification%arch_order
               if (time - lag >= 1) then
                  previous_shock = standardized(time - lag) - &
                     parameters%fgarch_shift(lag)
                  scale_value = scale_value + parameters%alpha(lag)* &
                     (sqrt(1.0e-6_dp + previous_shock**2) - &
                     parameters%asymmetry(lag)*previous_shock)** &
                     max(parameters%fgarch_delta + fgarch_fk( &
                     specification%fgarch_submodel)*parameters%fgarch_lambda, &
                     0.01_dp)*sqrt(variances(time - lag))** &
                     parameters%fgarch_lambda
               else
                  scale_value = scale_value + parameters%alpha(lag)* &
                     fgarch_innovation_moment(specification, parameters, lag)* &
                     initial_variance**(0.5_dp*parameters%fgarch_lambda)
               end if
            end do
            do lag = 1, specification%garch_order
               if (time - lag >= 1) then
                  scale_value = scale_value + parameters%beta(lag)* &
                     sqrt(variances(time - lag))**parameters%fgarch_lambda
               else
                  scale_value = scale_value + parameters%beta(lag)* &
                     initial_variance**(0.5_dp*parameters%fgarch_lambda)
               end if
            end do
            variance_value = max(scale_value, tiny(1.0_dp))** &
               (2.0_dp/parameters%fgarch_lambda)
         end select
         variances(time) = max(variance_value, tiny(1.0_dp))
         location = merge(parameters%mean, 0.0_dp, specification%include_mean)
         do lag = 1, specification%ar_order
            if (time - lag >= 1) location = location + parameters%ar(lag)* &
               (series(time - lag) - merge(parameters%mean, 0.0_dp, &
               specification%include_mean))
         end do
         do lag = 1, specification%ma_order
            if (time - lag >= 1) location = location + parameters%ma(lag)* &
               residuals(time - lag)
         end do
         if (specification%fractional_mean) then
            do lag = 1, min(time - 1, size(mean_fractional_weight))
               location = location - mean_fractional_weight(lag)* &
                  (series(time - lag) - merge(parameters%mean, 0.0_dp, &
                  specification%include_mean))
            end do
         end if
         if (specification%arch_in_mean) location = location + &
            parameters%arch_in_mean*sqrt(variances(time))** &
            specification%arch_in_mean_power
          standardized(time) = random_standardized_innovation( &
             specification%distribution, parameters%shape, parameters%skew, &
             parameters%lambda)
         residuals(time) = sqrt(variances(time))*standardized(time)
         means(time) = location
         series(time) = location + residuals(time)
         if (specification%variance_model == rugarch_model_realgarch) then
            measurement_errors(time) = parameters%measurement_sigma* &
               random_standard_normal()
            realized_path(time) = exp(parameters%measurement_intercept + &
               parameters%measurement_slope*log(variances(time)) + &
               parameters%measurement_leverage1*standardized(time) + &
               parameters%measurement_leverage2*(standardized(time)**2 - &
               1.0_dp) + measurement_errors(time))
         end if
      end do
      first = discarded + 1
      out%series = series(first:)
      out%conditional_mean = means(first:)
      out%conditional_variance = variances(first:)
      out%innovations = residuals(first:)
      out%realized_variance = realized_path(first:)
      out%measurement_errors = measurement_errors(first:)
   end function rugarch_simulate

   real(dp) function random_standardized_innovation(distribution, shape, skew, &
      lambda) &
      result(value)
      !! Draw one standardized innovation from a supported distribution.
      integer, intent(in) :: distribution !! Distribution code.
      real(dp), intent(in) :: shape !! Distribution shape parameter.
      real(dp), intent(in) :: skew !! Distribution skew parameter.
      real(dp), intent(in) :: lambda !! Generalized-hyperbolic lambda.
      real(dp) :: magnitude, scale, symmetric, first_moment
      real(dp) :: raw_mean, raw_scale, weight, inverse_shape, omega

      value = random_standardized(distribution, shape, skew, lambda)
      return
      select case (distribution)
      case (rugarch_distribution_normal)
         value = random_standard_normal()
      case (rugarch_distribution_student)
         value = random_standard_student(shape)
      case (rugarch_distribution_ged)
         scale = exp(0.5_dp*(log_gamma(1.0_dp/shape) - &
            log_gamma(3.0_dp/shape)))
         magnitude = scale*random_gamma(1.0_dp/shape)**(1.0_dp/shape)
         value = merge(magnitude, -magnitude, random_uniform() >= 0.5_dp)
      case (rugarch_distribution_skew_normal)
         symmetric = random_standard_normal()
         first_moment = sqrt(2.0_dp/acos(-1.0_dp))
         call fs_standardization(skew, first_moment, raw_mean, raw_scale)
         weight = skew/(skew + 1.0_dp/skew)
         value = merge(abs(symmetric)*skew, -abs(symmetric)/skew, &
            random_uniform() < weight)
         value = (value - raw_mean)/raw_scale
      case (rugarch_distribution_skew_student)
         symmetric = random_standard_student(shape)
         first_moment = exp(0.5_dp*log(shape - 2.0_dp) + &
            log_gamma(0.5_dp*(shape - 1.0_dp)) - &
            0.5_dp*log(acos(-1.0_dp)) - log_gamma(0.5_dp*shape))
         call fs_standardization(skew, first_moment, raw_mean, raw_scale)
         weight = skew/(skew + 1.0_dp/skew)
         value = merge(abs(symmetric)*skew, -abs(symmetric)/skew, &
            random_uniform() < weight)
         value = (value - raw_mean)/raw_scale
      case (rugarch_distribution_skew_ged)
         scale = exp(0.5_dp*(log_gamma(1.0_dp/shape) - &
            log_gamma(3.0_dp/shape)))
         magnitude = scale*random_gamma(1.0_dp/shape)**(1.0_dp/shape)
         symmetric = merge(magnitude, -magnitude, random_uniform() >= 0.5_dp)
         first_moment = scale*exp(log_gamma(2.0_dp/shape) - &
            log_gamma(1.0_dp/shape))
         call fs_standardization(skew, first_moment, raw_mean, raw_scale)
         weight = skew/(skew + 1.0_dp/skew)
         value = merge(abs(symmetric)*skew, -abs(symmetric)/skew, &
            random_uniform() < weight)
         value = (value - raw_mean)/raw_scale
      case (rugarch_distribution_johnson_su)
         inverse_shape = 1.0_dp/shape
         weight = exp(min(inverse_shape**2, 50.0_dp))
         omega = -skew*inverse_shape
         scale = sqrt(1.0_dp/(0.5_dp*(weight - 1.0_dp)* &
            (weight*cosh(2.0_dp*omega) + 1.0_dp)))
         value = scale*sqrt(weight)*sinh(omega) + &
            scale*sinh(inverse_shape*(random_standard_normal() + skew))
      case default
         value = 0.0_dp
      end select
   end function random_standardized_innovation

   pure function rugarch_berkowitz_test(data, lags) result(out)
      !! Translate BerkowitzTest for standardized density-transform observations.
      real(dp), intent(in) :: data(:) !! Normalized probability-integral-transform observations.
      integer, intent(in), optional :: lags !! Positive autoregressive test order.
      type(rugarch_berkowitz_test_t) :: out
      real(dp), allocatable :: centered(:), design(:, :), response(:)
      real(dp), allocatable :: cross_product(:, :), inverse(:, :), residuals(:)
      real(dp) :: second_moment, third_moment, fourth_moment
      integer :: selected_lags, observations, row, lag, status

      selected_lags = 1
      if (present(lags)) selected_lags = lags
      observations = size(data)
      out%lags = selected_lags
      if (selected_lags < 1 .or. observations <= selected_lags + 3 .or. &
         .not. all(ieee_is_finite(data))) then
         out%info = 1
         return
      end if
      out%mean = sum(data)/real(observations, dp)
      centered = data - out%mean
      allocate(response(observations - selected_lags), &
         design(observations - selected_lags, selected_lags))
      response = centered(selected_lags + 1:)
      do row = 1, size(response)
         do lag = 1, selected_lags
            design(row, lag) = centered(selected_lags + row - lag)
         end do
      end do
      cross_product = matmul(transpose(design), design)
      allocate(inverse(selected_lags, selected_lags))
      call symmetric_pseudoinverse(cross_product, inverse, status)
      if (status /= 0) then
         out%info = 2
         return
      end if
      out%autoregression = matmul(inverse, &
         matmul(transpose(design), response))
      residuals = response - matmul(design, out%autoregression)
      out%sigma = sqrt(max(sum(residuals**2)/ &
         real(max(1, size(residuals) - selected_lags), dp), tiny(1.0_dp)))
      out%unrestricted_log_likelihood = -0.5_dp*real(size(residuals), dp)* &
         log(2.0_dp*acos(-1.0_dp)*out%sigma**2) - &
         0.5_dp*sum((residuals/out%sigma)**2)
      out%restricted_log_likelihood = -0.5_dp*real(size(response), dp)* &
         log(2.0_dp*acos(-1.0_dp)) - 0.5_dp*sum(data(selected_lags + 1:)**2)
      out%likelihood_ratio = max(0.0_dp, 2.0_dp* &
         (out%unrestricted_log_likelihood - out%restricted_log_likelihood))
      out%likelihood_ratio_p_value = regularized_gamma_q( &
         0.5_dp*real(selected_lags + 2, dp), &
         0.5_dp*out%likelihood_ratio)
      second_moment = sum(centered**2)/real(observations, dp)
      third_moment = sum(centered**3)/real(observations, dp)
      fourth_moment = sum(centered**4)/real(observations, dp)
      if (second_moment > tiny(1.0_dp)) then
         out%jarque_bera = real(observations, dp)* &
            ((third_moment/second_moment**1.5_dp)**2/6.0_dp + &
            (fourth_moment/second_moment**2 - 3.0_dp)**2/24.0_dp)
         out%jarque_bera_p_value = regularized_gamma_q(1.0_dp, &
            0.5_dp*out%jarque_bera)
      end if
   end function rugarch_berkowitz_test

   pure function rugarch_direction_test(forecast, actual, method) result(out)
      !! Translate DACTest directional predictability statistics.
      real(dp), intent(in) :: forecast(:) !! Point forecasts.
      real(dp), intent(in) :: actual(:) !! Realized observations.
      integer, intent(in), optional :: method !! PT or AG test-method code.
      type(rugarch_direction_test_t) :: out
      real(dp) :: forecast_positive, actual_positive, correct, independent
      real(dp) :: correct_variance, independent_variance, denominator
      real(dp) :: signed_payoff, benchmark, probability_positive
      integer :: selected_method, observations

      selected_method = rugarch_direction_pt
      if (present(method)) selected_method = method
      out%method = selected_method
      observations = size(actual)
      if (size(forecast) /= observations .or. observations < 2 .or. &
         selected_method < rugarch_direction_pt .or. &
         selected_method > rugarch_direction_ag .or. &
         .not. all(ieee_is_finite(forecast)) .or. &
         .not. all(ieee_is_finite(actual))) then
         out%info = 1
         return
      end if
      if (selected_method == rugarch_direction_pt) then
         forecast_positive = real(count(forecast > 0.0_dp), dp)/ &
            real(observations, dp)
         actual_positive = real(count(actual > 0.0_dp), dp)/ &
            real(observations, dp)
         correct = real(count(forecast*actual > 0.0_dp), dp)/ &
            real(observations, dp)
         independent = forecast_positive*actual_positive + &
            (1.0_dp - forecast_positive)*(1.0_dp - actual_positive)
         correct_variance = independent*(1.0_dp - independent)/ &
            real(observations, dp)
         independent_variance = (2.0_dp*forecast_positive - 1.0_dp)**2* &
            actual_positive*(1.0_dp - actual_positive)/ &
            real(observations, dp) + &
            (2.0_dp*actual_positive - 1.0_dp)**2*forecast_positive* &
            (1.0_dp - forecast_positive)/real(observations, dp) + &
            4.0_dp*forecast_positive*actual_positive* &
            (1.0_dp - forecast_positive)*(1.0_dp - actual_positive)/ &
            real(observations, dp)**2
         denominator = sqrt(max(correct_variance - independent_variance, &
            tiny(1.0_dp)))
         out%statistic = (correct - independent)/denominator
         out%directional_accuracy = correct
      else
         signed_payoff = sum(signed_unit(forecast)*actual)/ &
            real(observations, dp)
         benchmark = sum(signed_unit(forecast))/real(observations, dp)* &
            sum(actual)/real(observations, dp)
         probability_positive = 0.5_dp*(1.0_dp + &
            sum(signed_unit(forecast))/real(observations, dp))
         denominator = sqrt(max(4.0_dp/real(observations, dp)**2* &
            probability_positive*(1.0_dp - probability_positive)* &
            sum((actual - sum(actual)/real(observations, dp))**2), &
            tiny(1.0_dp)))
         out%statistic = (signed_payoff - benchmark)/denominator
         out%directional_accuracy = real(count(signed_unit(forecast)* &
            actual > 0.0_dp), dp)/real(observations, dp)
      end if
      out%p_value = 0.5_dp*erfc(out%statistic/sqrt(2.0_dp))
   end function rugarch_direction_test

   pure elemental real(dp) function signed_unit(value) result(direction)
      !! Return minus one, zero, or one according to a real value's sign.
      real(dp), intent(in) :: value !! Input value.

      if (value > 0.0_dp) then
         direction = 1.0_dp
      else if (value < 0.0_dp) then
         direction = -1.0_dp
      else
         direction = 0.0_dp
      end if
   end function signed_unit

   pure function rugarch_var_test(alpha, actual, value_at_risk, confidence) &
      result(out)
      !! Translate VaRTest coverage and exceedance-independence tests.
      real(dp), intent(in) :: alpha !! Lower-tail VaR probability.
      real(dp), intent(in) :: actual(:) !! Realized observations.
      real(dp), intent(in) :: value_at_risk(:) !! VaR forecasts.
      real(dp), intent(in), optional :: confidence !! Test confidence level.
      type(rugarch_var_test_t) :: out
      logical, allocatable :: exceedance(:)
      real(dp) :: selected_confidence, observed_probability
      real(dp) :: probability_01, probability_11, null_log_likelihood
      real(dp) :: independent_log_likelihood, markov_log_likelihood
      integer :: observations, n00, n01, n10, n11, time

      selected_confidence = 0.95_dp
      if (present(confidence)) selected_confidence = confidence
      observations = size(actual)
      if (size(value_at_risk) /= observations .or. observations < 2 .or. &
         alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. &
         selected_confidence <= 0.0_dp .or. selected_confidence >= 1.0_dp .or. &
         .not. all(ieee_is_finite(actual)) .or. &
         .not. all(ieee_is_finite(value_at_risk))) then
         out%info = 1
         return
      end if
      exceedance = actual < value_at_risk
      out%expected_exceedances = floor(alpha*real(observations, dp))
      out%actual_exceedances = count(exceedance)
      observed_probability = real(out%actual_exceedances, dp)/ &
         real(observations, dp)
      null_log_likelihood = xlogy(observations - out%actual_exceedances, &
         1.0_dp - alpha) + xlogy(out%actual_exceedances, alpha)
      independent_log_likelihood = xlogy(observations - &
         out%actual_exceedances, 1.0_dp - observed_probability) + &
         xlogy(out%actual_exceedances, observed_probability)
      out%unconditional_statistic = max(0.0_dp, 2.0_dp* &
         (independent_log_likelihood - null_log_likelihood))
      out%unconditional_p_value = regularized_gamma_q(0.5_dp, &
         0.5_dp*out%unconditional_statistic)
      n00 = 0
      n01 = 0
      n10 = 0
      n11 = 0
      do time = 2, observations
         if (.not. exceedance(time - 1) .and. .not. exceedance(time)) n00 = n00 + 1
         if (.not. exceedance(time - 1) .and. exceedance(time)) n01 = n01 + 1
         if (exceedance(time - 1) .and. .not. exceedance(time)) n10 = n10 + 1
         if (exceedance(time - 1) .and. exceedance(time)) n11 = n11 + 1
      end do
      probability_01 = real(n01, dp)/real(max(1, n00 + n01), dp)
      probability_11 = real(n11, dp)/real(max(1, n10 + n11), dp)
      markov_log_likelihood = xlogy(n00, 1.0_dp - probability_01) + &
         xlogy(n01, probability_01) + xlogy(n10, 1.0_dp - probability_11) + &
         xlogy(n11, probability_11)
      independent_log_likelihood = xlogy(n00 + n10, &
         1.0_dp - observed_probability) + xlogy(n01 + n11, &
         observed_probability)
      out%conditional_statistic = out%unconditional_statistic + &
         max(0.0_dp, 2.0_dp*(markov_log_likelihood - &
         independent_log_likelihood))
      out%conditional_p_value = regularized_gamma_q(1.0_dp, &
         0.5_dp*out%conditional_statistic)
      out%unconditional_critical = chi_square_quantile_approx( &
         selected_confidence, 1)
      out%conditional_critical = chi_square_quantile_approx( &
         selected_confidence, 2)
   end function rugarch_var_test

   pure function rugarch_es_test(alpha, actual, expected_shortfall, &
      value_at_risk) result(out)
      !! Translate ESTest's one-sided excess-violation mean test.
      real(dp), intent(in) :: alpha !! Lower-tail ES probability.
      real(dp), intent(in) :: actual(:) !! Realized observations.
      real(dp), intent(in) :: expected_shortfall(:) !! Expected-shortfall forecasts.
      real(dp), intent(in) :: value_at_risk(:) !! VaR forecasts.
      type(rugarch_es_test_t) :: out
      real(dp), allocatable :: excess(:)
      logical, allocatable :: exceedance(:)
      real(dp) :: mean_excess, standard_error
      integer :: observations

      observations = size(actual)
      if (size(expected_shortfall) /= observations .or. &
         size(value_at_risk) /= observations .or. observations < 2 .or. &
         alpha <= 0.0_dp .or. alpha >= 1.0_dp .or. &
         .not. all(ieee_is_finite(actual)) .or. &
         .not. all(ieee_is_finite(expected_shortfall)) .or. &
         .not. all(ieee_is_finite(value_at_risk))) then
         out%info = 1
         return
      end if
      exceedance = actual < value_at_risk
      out%expected_exceedances = floor(alpha*real(observations, dp))
      out%actual_exceedances = count(exceedance)
      if (out%actual_exceedances < 2) then
         out%info = 2
         return
      end if
      excess = pack(expected_shortfall - actual, exceedance)
      mean_excess = sum(excess)/real(size(excess), dp)
      standard_error = sqrt(max(sum((excess - mean_excess)**2)/ &
         real(size(excess), dp), tiny(1.0_dp))/ &
         real(size(excess) - 1, dp))
      out%statistic = mean_excess/standard_error
      out%p_value = 0.5_dp*erfc(out%statistic/sqrt(2.0_dp))
   end function rugarch_es_test

   pure elemental real(dp) function xlogy(count_value, probability) result(value)
      !! Return count times log probability with the zero-count convention.
      integer, intent(in) :: count_value !! Nonnegative count.
      real(dp), intent(in) :: probability !! Probability value.

      if (count_value == 0) then
         value = 0.0_dp
      else
         value = real(count_value, dp)*log(max(probability, tiny(1.0_dp)))
      end if
   end function xlogy

   pure elemental real(dp) function chi_square_quantile_approx(probability, &
      degrees_of_freedom) result(value)
      !! Approximate a chi-square quantile using Wilson-Hilferty transformation.
      real(dp), intent(in) :: probability !! Cumulative probability.
      integer, intent(in) :: degrees_of_freedom !! Positive degrees of freedom.
      real(dp) :: degrees

      degrees = real(max(1, degrees_of_freedom), dp)
      value = degrees*(1.0_dp - 2.0_dp/(9.0_dp*degrees) + &
         normal_quantile(probability)*sqrt(2.0_dp/(9.0_dp*degrees)))**3
   end function chi_square_quantile_approx

   pure function rugarch_model_name(model) result(name)
      !! Return the rugarch variance-model label.
      integer, intent(in) :: model !! Variance-model code.
      character(len=12) :: name

      select case (model)
      case (rugarch_model_sgarch)
         name = 'sGARCH'
      case (rugarch_model_igarch)
         name = 'iGARCH'
      case (rugarch_model_egarch)
         name = 'eGARCH'
      case (rugarch_model_gjrgarch)
         name = 'gjrGARCH'
      case (rugarch_model_aparch)
         name = 'apARCH'
      case (rugarch_model_figarch)
         name = 'fiGARCH'
      case (rugarch_model_csgarch)
         name = 'csGARCH'
      case (rugarch_model_realgarch)
         name = 'realGARCH'
      case (rugarch_model_fgarch)
         name = 'fGARCH'
      case default
         name = 'unknown'
      end select
   end function rugarch_model_name

   pure function rugarch_distribution_name(distribution) result(name)
      !! Return the rugarch innovation-distribution label.
      integer, intent(in) :: distribution !! Distribution code.
      character(len=12) :: name

      name = shared_distribution_name(distribution)
      return
      select case (distribution)
      case (rugarch_distribution_normal)
         name = 'norm'
      case (rugarch_distribution_student)
         name = 'std'
      case (rugarch_distribution_ged)
         name = 'ged'
      case (rugarch_distribution_skew_normal)
         name = 'snorm'
      case (rugarch_distribution_skew_student)
         name = 'sstd'
      case (rugarch_distribution_skew_ged)
         name = 'sged'
      case (rugarch_distribution_johnson_su)
         name = 'jsu'
      case default
         name = 'unknown'
      end select
   end function rugarch_distribution_name

   subroutine display_rugarch_spec(value, unit)
      !! Display a rugarch model specification.
      type(rugarch_spec_t), intent(in) :: value !! Model specification.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'rugarch specification'
      write(destination, '(a,a)') '  variance model: ', &
         trim(rugarch_model_name(value%variance_model))
      write(destination, '(a,a)') '  distribution: ', &
         trim(rugarch_distribution_name(value%distribution))
      write(destination, '(a,2(1x,i0))') '  variance order:', &
         value%arch_order, value%garch_order
      write(destination, '(a,2(1x,i0))') '  ARMA order:', &
         value%ar_order, value%ma_order
      write(destination, '(a,l1)') '  include mean: ', value%include_mean
      write(destination, '(a,l1)') '  ARCH in mean: ', value%arch_in_mean
      write(destination, '(a,l1)') '  fractional mean: ', value%fractional_mean
      write(destination, '(a,i0)') '  mean regressors: ', &
         value%mean_regressor_count
      write(destination, '(a,i0)') '  variance regressors: ', &
         value%variance_regressor_count
      if (value%variance_model == rugarch_model_fgarch) &
         write(destination, '(a,i0)') '  fGARCH submodel: ', &
         value%fgarch_submodel
      if (value%variance_model == rugarch_model_figarch) &
         write(destination, '(a,i0)') '  truncation lag: ', value%truncation_lag
   end subroutine display_rugarch_spec

   subroutine display_rugarch_parameters(value, unit)
      !! Display physical rugarch parameters.
      type(rugarch_parameters_t), intent(in) :: value !! Physical parameters.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'rugarch parameters'
      write(destination, '(a,es14.6)') '  mean: ', value%mean
      write(destination, '(a,es14.6)') '  ARCH-in-mean: ', value%arch_in_mean
      write(destination, '(a,es14.6)') '  omega: ', value%omega
      if (allocated(value%ar)) &
         write(destination, '(a,*(1x,es14.6))') '  AR:', value%ar
      if (allocated(value%ma)) &
         write(destination, '(a,*(1x,es14.6))') '  MA:', value%ma
      if (allocated(value%alpha)) &
         write(destination, '(a,*(1x,es14.6))') '  alpha:', value%alpha
      if (allocated(value%beta)) &
         write(destination, '(a,*(1x,es14.6))') '  beta:', value%beta
      if (allocated(value%asymmetry)) &
         write(destination, '(a,*(1x,es14.6))') '  asymmetry:', &
            value%asymmetry
      write(destination, '(a,es14.6)') '  power: ', value%power
      write(destination, '(a,es14.6)') '  skew: ', value%skew
      write(destination, '(a,es14.6)') '  shape: ', value%shape
      write(destination, '(a,es14.6)') '  lambda: ', value%lambda
      write(destination, '(a,es14.6)') '  component rho: ', &
         value%component_rho
      write(destination, '(a,es14.6)') '  component phi: ', &
         value%component_phi
      write(destination, '(a,es14.6)') '  measurement intercept: ', &
         value%measurement_intercept
      write(destination, '(a,es14.6)') '  measurement slope: ', &
         value%measurement_slope
      write(destination, '(a,es14.6)') '  measurement leverage 1: ', &
         value%measurement_leverage1
      write(destination, '(a,es14.6)') '  measurement leverage 2: ', &
         value%measurement_leverage2
      write(destination, '(a,es14.6)') '  measurement sigma: ', &
         value%measurement_sigma
      write(destination, '(a,es14.6)') '  fractional d: ', value%fractional
      write(destination, '(a,es14.6)') '  mean fractional d: ', &
         value%mean_fractional
      write(destination, '(a,es14.6)') '  fGARCH delta: ', value%fgarch_delta
      write(destination, '(a,es14.6)') '  fGARCH lambda: ', value%fgarch_lambda
      if (allocated(value%fgarch_shift)) &
         write(destination, '(a,*(1x,es14.6))') '  fGARCH shift:', &
         value%fgarch_shift
      if (allocated(value%mean_regression)) &
         write(destination, '(a,*(1x,es14.6))') '  mean regression:', &
         value%mean_regression
      if (allocated(value%variance_regression)) &
         write(destination, '(a,*(1x,es14.6))') '  variance regression:', &
         value%variance_regression
   end subroutine display_rugarch_parameters

   subroutine display_rugarch_filter(value, unit, print_obs)
      !! Display a filtered model and optionally its observation paths.
      type(rugarch_filter_t), intent(in) :: value !! Filter result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print conditional path values.
      integer :: destination, time
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'rugarch filter'
      write(destination, '(a,es14.6)') '  log likelihood: ', &
         value%log_likelihood
      write(destination, '(a,i0)') '  status: ', value%info
      if (show_observations .and. allocated(value%conditional_mean)) then
         write(destination, '(a)') '  time mean sigma residual standardized'
         do time = 1, size(value%conditional_mean)
            write(destination, '(i8,4(1x,es14.6))') time, &
               value%conditional_mean(time), value%conditional_sigma(time), &
               value%residuals(time), value%standardized_residuals(time)
         end do
      end if
   end subroutine display_rugarch_filter

   subroutine display_rugarch_fit(value, unit, print_obs)
      !! Display estimated rugarch parameters and optional filtered paths.
      type(rugarch_fit_t), intent(in) :: value !! Fit result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print filtered observation paths.
      integer :: destination
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'rugarch fit'
      write(destination, '(a,a)') '  variance model: ', &
         trim(rugarch_model_name(value%specification%variance_model))
      write(destination, '(a,a)') '  distribution: ', &
         trim(rugarch_distribution_name(value%specification%distribution))
      write(destination, '(a,es14.6)') '  log likelihood: ', value%log_likelihood
      write(destination, '(a,3(1x,es14.6))') '  AIC BIC HQIC:', &
         value%aic, value%bic, value%hqic
      write(destination, '(a,es14.6)') '  persistence: ', value%persistence
      write(destination, '(a,l1)') '  converged: ', value%converged
      write(destination, '(a,i0)') '  status: ', value%info
      call display_rugarch_parameters(value%parameters, destination)
      if (show_observations) call display_rugarch_filter(value%filtered, &
         destination, .true.)
   end subroutine display_rugarch_fit

   subroutine display_rugarch_forecast(value, unit)
      !! Display recursive rugarch forecasts.
      type(rugarch_forecast_t), intent(in) :: value !! Forecast result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination, horizon

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'rugarch forecast'
      if (allocated(value%mean)) then
         write(destination, '(a)') '  horizon mean sigma variance'
         do horizon = 1, size(value%mean)
            write(destination, '(i8,3(1x,es14.6))') horizon, &
               value%mean(horizon), value%sigma(horizon), value%variance(horizon)
         end do
      end if
      write(destination, '(a,i0)') '  status: ', value%info
   end subroutine display_rugarch_forecast

   subroutine display_rugarch_simulation(value, unit, print_obs)
      !! Display simulation settings and optionally simulated observations.
      type(rugarch_simulation_t), intent(in) :: value !! Simulation result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      logical, intent(in), optional :: print_obs !! Print simulated observation paths.
      integer :: destination, time
      logical :: show_observations

      destination = output_unit
      if (present(unit)) destination = unit
      show_observations = .false.
      if (present(print_obs)) show_observations = print_obs
      write(destination, '(a)') 'rugarch simulation'
      write(destination, '(a,a)') '  variance model: ', &
         trim(rugarch_model_name(value%specification%variance_model))
      write(destination, '(a,i0)') '  burn-in: ', value%burnin
      if (allocated(value%series)) write(destination, '(a,i0)') &
         '  observations: ', size(value%series)
      write(destination, '(a,i0)') '  status: ', value%info
      if (show_observations .and. allocated(value%series)) then
         write(destination, '(a)') '  time observation mean sigma innovation'
         do time = 1, size(value%series)
            write(destination, '(i8,4(1x,es14.6))') time, value%series(time), &
               value%conditional_mean(time), &
               sqrt(value%conditional_variance(time)), value%innovations(time)
         end do
      end if
   end subroutine display_rugarch_simulation

   subroutine display_rugarch_berkowitz_test(value, unit)
      !! Display Berkowitz and Jarque-Bera calibration statistics.
      type(rugarch_berkowitz_test_t), intent(in) :: value !! Berkowitz test result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Berkowitz density calibration test'
      write(destination, '(a,es14.6,a,es14.6)') '  LR statistic: ', &
         value%likelihood_ratio, ', p-value: ', value%likelihood_ratio_p_value
      write(destination, '(a,es14.6,a,es14.6)') '  Jarque-Bera: ', &
         value%jarque_bera, ', p-value: ', value%jarque_bera_p_value
      write(destination, '(a,es14.6)') '  mean: ', value%mean
      write(destination, '(a,es14.6)') '  sigma: ', value%sigma
      if (allocated(value%autoregression)) &
         write(destination, '(a,*(1x,es14.6))') '  AR:', value%autoregression
      write(destination, '(a,i0)') '  status: ', value%info
   end subroutine display_rugarch_berkowitz_test

   subroutine display_rugarch_direction_test(value, unit)
      !! Display a directional-accuracy test result.
      type(rugarch_direction_test_t), intent(in) :: value !! Directional test result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Directional accuracy test'
      write(destination, '(a,i0)') '  method: ', value%method
      write(destination, '(a,es14.6)') '  statistic: ', value%statistic
      write(destination, '(a,es14.6)') '  p-value: ', value%p_value
      write(destination, '(a,es14.6)') '  directional accuracy: ', &
         value%directional_accuracy
      write(destination, '(a,i0)') '  status: ', value%info
   end subroutine display_rugarch_direction_test

   subroutine display_rugarch_var_test(value, unit)
      !! Display unconditional and conditional VaR coverage tests.
      type(rugarch_var_test_t), intent(in) :: value !! VaR test result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Value-at-Risk coverage test'
      write(destination, '(a,2(1x,i0))') '  expected and actual exceedances:', &
         value%expected_exceedances, value%actual_exceedances
      write(destination, '(a,es14.6,a,es14.6)') '  unconditional LR: ', &
         value%unconditional_statistic, ', p-value: ', &
         value%unconditional_p_value
      write(destination, '(a,es14.6,a,es14.6)') '  conditional LR: ', &
         value%conditional_statistic, ', p-value: ', value%conditional_p_value
      write(destination, '(a,i0)') '  status: ', value%info
   end subroutine display_rugarch_var_test

   subroutine display_rugarch_es_test(value, unit)
      !! Display a conditional expected-shortfall test.
      type(rugarch_es_test_t), intent(in) :: value !! Expected-shortfall test result.
      integer, intent(in), optional :: unit !! Output unit; defaults to standard output.
      integer :: destination

      destination = output_unit
      if (present(unit)) destination = unit
      write(destination, '(a)') 'Expected-shortfall test'
      write(destination, '(a,2(1x,i0))') '  expected and actual exceedances:', &
         value%expected_exceedances, value%actual_exceedances
      write(destination, '(a,es14.6)') '  statistic: ', value%statistic
      write(destination, '(a,es14.6)') '  p-value: ', value%p_value
      write(destination, '(a,i0)') '  status: ', value%info
   end subroutine display_rugarch_es_test

end module rugarch_mod
