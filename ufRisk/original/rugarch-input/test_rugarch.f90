! SPDX-License-Identifier: MIT
! SPDX-FileComment: Tests for algorithms translated from the R rugarch package.
program test_rugarch
   use kind_mod, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use rugarch_mod, only: rugarch_spec_t, rugarch_parameters_t, &
      rugarch_filter_t, rugarch_fit_t, rugarch_forecast_t, &
      rugarch_simulation_t, rugarch_berkowitz_test_t, &
      rugarch_direction_test_t, rugarch_var_test_t, rugarch_es_test_t, &
      rugarch_model_sgarch, rugarch_model_igarch, &
      rugarch_model_egarch, rugarch_model_gjrgarch, rugarch_model_aparch, &
      rugarch_model_figarch, rugarch_model_csgarch, rugarch_model_realgarch, &
      rugarch_model_fgarch, &
      rugarch_distribution_normal, rugarch_distribution_student, &
      rugarch_distribution_ged, rugarch_distribution_skew_normal, &
      rugarch_distribution_skew_student, rugarch_distribution_skew_ged, &
      rugarch_distribution_johnson_su, rugarch_distribution_nig, &
      rugarch_distribution_ghyp, rugarch_distribution_gh_skew_student, &
      rugarch_spec, rugarch_filter, rugarch_fit, &
      rugarch_forecast, rugarch_simulate, rugarch_log_density, &
      rugarch_berkowitz_test, rugarch_direction_test, rugarch_var_test, &
      rugarch_es_test
   use rugarch_extensions_mod, only: rugarch_bootstrap_t, &
      rugarch_roll_t, rugarch_multifit_t, &
      rugarch_multifilter_t, rugarch_multiforecast_t, &
      rugarch_parameter_distribution_t, rugarch_model_confidence_t, &
      rugarch_news_impact_t, rugarch_sign_bias_test_t, &
      rugarch_bootstrap_forecast, rugarch_parameter_distribution, &
      rugarch_model_confidence, rugarch_news_impact, rugarch_sign_bias_test, &
      rugarch_roll, rugarch_multifit
   use rugarch_extensions_mod, only: rugarch_multifilter, &
      rugarch_multiforecast
   use rugarch_diagnostics_mod, only: rugarch_nyblom_test_t, &
      rugarch_gof_test_t, rugarch_var_duration_test_t, rugarch_gmm_test_t, &
      rugarch_hong_li_test_t, rugarch_nyblom_test, rugarch_gof_test, &
      rugarch_var_duration_test, rugarch_gmm_test, rugarch_hong_li_test
   use distribution_mod, only: standardized_cdf
   implicit none
   integer, parameter :: model_count = 8
   integer, parameter :: distribution_count = 10
   integer :: models(model_count), distributions(model_count), model
   integer :: distribution_codes(distribution_count), distribution
   type(rugarch_spec_t) :: specification
   type(rugarch_parameters_t) :: parameters
   type(rugarch_simulation_t) :: simulation
   type(rugarch_filter_t) :: filtered
   type(rugarch_fit_t) :: fitted
   type(rugarch_forecast_t) :: forecast
   type(rugarch_berkowitz_test_t) :: berkowitz
   type(rugarch_direction_test_t) :: direction
   type(rugarch_var_test_t) :: value_at_risk_test
   type(rugarch_es_test_t) :: expected_shortfall_test
   type(rugarch_bootstrap_t) :: bootstrap
   type(rugarch_roll_t) :: rolling
   type(rugarch_multifit_t) :: multiple_fit
   type(rugarch_multifilter_t) :: multiple_filter
   type(rugarch_multiforecast_t) :: multiple_forecast
   type(rugarch_parameter_distribution_t) :: parameter_distribution
   type(rugarch_model_confidence_t) :: confidence_set
   type(rugarch_news_impact_t) :: news_impact
   type(rugarch_sign_bias_test_t) :: sign_bias
   type(rugarch_nyblom_test_t) :: nyblom
   type(rugarch_gof_test_t) :: goodness_of_fit
   type(rugarch_var_duration_test_t) :: duration_test
   type(rugarch_gmm_test_t) :: gmm_test
   type(rugarch_hong_li_test_t) :: hong_li_test
   real(dp) :: value_at_risk(300), expected_shortfall(300)
   real(dp) :: losses(20, 2), shocks(5)
   type(rugarch_spec_t) :: specifications(2)
   real(dp) :: mean_regressors(180, 1), variance_regressors(180, 1)
   real(dp) :: probability_transform(300)

   models = [rugarch_model_sgarch, rugarch_model_igarch, &
      rugarch_model_egarch, rugarch_model_gjrgarch, rugarch_model_aparch, &
      rugarch_model_figarch, rugarch_model_csgarch, rugarch_model_fgarch]
   distributions = [rugarch_distribution_normal, rugarch_distribution_student, &
      rugarch_distribution_ged, rugarch_distribution_normal, &
      rugarch_distribution_student, rugarch_distribution_normal, &
      rugarch_distribution_skew_normal, rugarch_distribution_normal]
   do model = 1, model_count
      specification = rugarch_spec(variance_model=models(model), &
         distribution=distributions(model), include_mean=.false.)
      parameters = test_parameters(specification)
      simulation = rugarch_simulate(specification, parameters, 240, &
         burnin=100, seed=700 + model)
      if (simulation%info /= 0) error stop 'rugarch simulation failed'
      if (any(.not. ieee_is_finite(simulation%series))) &
         error stop 'nonfinite rugarch simulation'
      filtered = rugarch_filter(simulation%series, specification, parameters)
      if (filtered%info /= 0 .or. &
         .not. ieee_is_finite(filtered%log_likelihood)) &
         error stop 'rugarch filter failed'
      if (any(filtered%conditional_variance <= 0.0_dp)) &
         error stop 'nonpositive rugarch variance'
   end do

   distribution_codes = [rugarch_distribution_normal, &
      rugarch_distribution_student, rugarch_distribution_ged, &
      rugarch_distribution_skew_normal, rugarch_distribution_skew_student, &
      rugarch_distribution_skew_ged, rugarch_distribution_johnson_su, &
      rugarch_distribution_nig, rugarch_distribution_ghyp, &
      rugarch_distribution_gh_skew_student]
   do distribution = 1, distribution_count
      specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
         distribution=distribution_codes(distribution), include_mean=.false.)
      parameters = test_parameters(specification)
      simulation = rugarch_simulate(specification, parameters, 300, &
         burnin=150, seed=900 + distribution)
      filtered = rugarch_filter(simulation%series, specification, parameters)
      if (simulation%info /= 0 .or. filtered%info /= 0 .or. &
         .not. ieee_is_finite(filtered%log_likelihood)) &
         error stop 'rugarch innovation distribution failed'
      if (abs(sum(simulation%innovations/ &
         sqrt(simulation%conditional_variance))/ &
         real(size(simulation%innovations), dp)) > 0.20_dp) &
         error stop 'rugarch innovation centering failed'
   end do

   specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
      distribution=rugarch_distribution_skew_normal, include_mean=.false.)
   parameters = test_parameters(specification)
   simulation = rugarch_simulate(specification, parameters, 400, &
      burnin=200, seed=1776)
   fitted = rugarch_fit(simulation%series, specification, max_iterations=160)
   if (fitted%info /= 0 .or. fitted%parameters%skew <= 0.0_dp .or. &
      .not. ieee_is_finite(fitted%log_likelihood)) &
      error stop 'skewed rugarch fit failed'

   specification = rugarch_spec(variance_model=rugarch_model_realgarch, &
      distribution=rugarch_distribution_normal, include_mean=.false.)
   parameters = test_parameters(specification)
   simulation = rugarch_simulate(specification, parameters, 320, &
      burnin=200, seed=1948)
   filtered = rugarch_filter(simulation%series, specification, parameters, &
      simulation%realized_variance)
   if (simulation%info /= 0 .or. filtered%info /= 0 .or. &
      any(filtered%realized_variance <= 0.0_dp) .or. &
      .not. ieee_is_finite(filtered%log_likelihood)) &
      error stop 'realized rugarch filter failed'
   fitted = rugarch_fit(simulation%series, specification, initial=parameters, &
      max_iterations=100, realized_variance=simulation%realized_variance)
   if (fitted%info /= 0 .or. fitted%parameters%measurement_sigma <= 0.0_dp) &
      error stop 'realized rugarch fit failed'
   forecast = rugarch_forecast(fitted, 5)
   if (forecast%info /= 0 .or. any(forecast%realized_variance <= 0.0_dp)) &
      error stop 'realized rugarch forecast failed'

   specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
      distribution=rugarch_distribution_normal, ar_order=1, ma_order=1, &
      include_mean=.true., arch_in_mean=.true.)
   parameters = test_parameters(specification)
   parameters%mean = 0.04_dp
   parameters%ar = 0.20_dp
   parameters%ma = -0.10_dp
   parameters%arch_in_mean = 0.03_dp
   simulation = rugarch_simulate(specification, parameters, 180, &
      burnin=100, seed=812)
   filtered = rugarch_filter(simulation%series, specification, parameters)
   if (simulation%info /= 0 .or. filtered%info /= 0) &
      error stop 'ARMA ARCH-in-mean path failed'
   specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
      distribution=rugarch_distribution_normal, include_mean=.true., &
      fractional_mean=.true., mean_regressor_count=1, &
      variance_regressor_count=1)
   parameters = test_parameters(specification)
   parameters%mean_fractional = 0.20_dp
   parameters%mean_regression = 0.03_dp
   parameters%variance_regression = 0.002_dp
   mean_regressors(:, 1) = [(sin(0.05_dp*real(model, dp)), model=1, 180)]
   variance_regressors(:, 1) = [(cos(0.05_dp*real(model, dp)), model=1, 180)]
   filtered = rugarch_filter(simulation%series, specification, parameters, &
      mean_regressors=mean_regressors, &
      variance_regressors=variance_regressors)
   if (filtered%info /= 0 .or. &
      .not. ieee_is_finite(filtered%log_likelihood)) &
      error stop 'ARFIMA regressor path failed'

   specification = rugarch_spec(variance_model=rugarch_model_sgarch, &
      distribution=rugarch_distribution_normal, include_mean=.false.)
   parameters = test_parameters(specification)
   simulation = rugarch_simulate(specification, parameters, 300, &
      burnin=200, seed=42)
   fitted = rugarch_fit(simulation%series, specification, max_iterations=120)
   if (fitted%info /= 0) error stop 'rugarch fit failed'
   if (fitted%parameters%omega <= 0.0_dp .or. &
      sum(fitted%parameters%alpha) + sum(fitted%parameters%beta) >= 1.0_dp) &
      error stop 'invalid fitted rugarch parameters'
   forecast = rugarch_forecast(fitted, 6)
   if (forecast%info /= 0 .or. size(forecast%mean) /= 6 .or. &
      any(forecast%variance <= 0.0_dp)) error stop 'rugarch forecast failed'
   if (.not. ieee_is_finite(rugarch_log_density(0.2_dp, &
      rugarch_distribution_student, 8.0_dp))) &
      error stop 'rugarch Student density failed'
   if (.not. ieee_is_finite(rugarch_log_density(0.2_dp, &
      rugarch_distribution_ged, 1.5_dp))) &
      error stop 'rugarch GED density failed'
   if (abs(rugarch_log_density(0.2_dp, rugarch_distribution_normal, 0.0_dp) - &
      rugarch_log_density(0.2_dp, rugarch_distribution_skew_normal, 0.0_dp, &
      1.0_dp)) > 1.0e-12_dp) error stop 'skew-normal nesting failed'
   if (abs(rugarch_log_density(0.2_dp, rugarch_distribution_student, 8.0_dp) - &
      rugarch_log_density(0.2_dp, rugarch_distribution_skew_student, 8.0_dp, &
      1.0_dp)) > 1.0e-12_dp) error stop 'skew-Student nesting failed'
   if (abs(rugarch_log_density(0.2_dp, rugarch_distribution_ged, 1.5_dp) - &
      rugarch_log_density(0.2_dp, rugarch_distribution_skew_ged, 1.5_dp, &
      1.0_dp)) > 1.0e-12_dp) error stop 'skew-GED nesting failed'
   berkowitz = rugarch_berkowitz_test(fitted%filtered%standardized_residuals)
   if (berkowitz%info /= 0 .or. &
      .not. ieee_is_finite(berkowitz%likelihood_ratio)) &
      error stop 'rugarch Berkowitz test failed'
   direction = rugarch_direction_test(forecast=fitted%filtered%conditional_mean, &
      actual=simulation%series)
   if (direction%info /= 0 .or. direction%p_value < 0.0_dp .or. &
      direction%p_value > 1.0_dp) error stop 'rugarch direction test failed'
   value_at_risk = fitted%filtered%conditional_mean - &
      1.645_dp*fitted%filtered%conditional_sigma
   expected_shortfall = fitted%filtered%conditional_mean - &
      2.063_dp*fitted%filtered%conditional_sigma
   value_at_risk_test = rugarch_var_test(0.05_dp, simulation%series, &
      value_at_risk)
   if (value_at_risk_test%info /= 0) error stop 'rugarch VaR test failed'
   expected_shortfall_test = rugarch_es_test(0.05_dp, simulation%series, &
      expected_shortfall, value_at_risk)
   if (expected_shortfall_test%info /= 0) &
      error stop 'rugarch ES test failed'
   bootstrap = rugarch_bootstrap_forecast(fitted, 4, 100, seed=912)
   if (bootstrap%info /= 0 .or. any(bootstrap%upper < bootstrap%lower)) &
      error stop 'rugarch bootstrap failed'
   parameter_distribution = rugarch_parameter_distribution(fitted, 40, seed=91)
   if (parameter_distribution%info /= 0 .or. &
      size(parameter_distribution%draws, 2) /= 40) &
      error stop 'rugarch parameter distribution failed'
   losses(:, 1) = 1.0_dp
   losses(:, 2) = [(1.0_dp + 0.01_dp*real(model, dp), model=1, 20)]
   confidence_set = rugarch_model_confidence(losses)
   if (confidence_set%info /= 0 .or. confidence_set%best_model /= 1) &
      error stop 'rugarch model confidence failed'
   shocks = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   news_impact = rugarch_news_impact(fitted%specification, &
      fitted%parameters, shocks)
   if (news_impact%info /= 0 .or. any(news_impact%variance <= 0.0_dp)) &
      error stop 'rugarch news impact failed'
   sign_bias = rugarch_sign_bias_test(fitted%filtered%standardized_residuals)
   if (sign_bias%info /= 0 .or. sign_bias%joint_p_value < 0.0_dp .or. &
      sign_bias%joint_p_value > 1.0_dp) error stop 'rugarch sign bias failed'
   specifications(1) = specification
   specifications(2) = rugarch_spec(variance_model=rugarch_model_sgarch, &
      distribution=rugarch_distribution_student, include_mean=.false.)
   multiple_fit = rugarch_multifit(simulation%series, specifications, &
      max_iterations=30)
   if (multiple_fit%info /= 0 .or. multiple_fit%best_aic < 1) &
      error stop 'rugarch multifit failed'
   multiple_filter = rugarch_multifilter(simulation%series, specifications, &
      multiple_fit%fits%parameters)
   if (multiple_filter%info /= 0) error stop 'rugarch multifilter failed'
   multiple_forecast = rugarch_multiforecast(multiple_fit, 3)
   if (multiple_forecast%info /= 0 .or. &
      size(multiple_forecast%forecasts) /= 2) &
      error stop 'rugarch multiforecast failed'
   rolling = rugarch_roll(simulation%series, specification, 280, &
      refit_every=10, max_iterations=30)
   if (rolling%info /= 0 .or. size(rolling%mean) /= 20) &
      error stop 'rugarch rolling fit failed'
   nyblom = rugarch_nyblom_test(fitted)
   if (nyblom%info /= 0 .or. .not. ieee_is_finite( &
      nyblom%joint_statistic)) error stop 'rugarch Nyblom test failed'
   goodness_of_fit = rugarch_gof_test( &
      fitted%filtered%standardized_residuals, rugarch_distribution_normal, &
      [10, 20])
   if (goodness_of_fit%info /= 0 .or. &
      any(goodness_of_fit%p_value < 0.0_dp) .or. &
      any(goodness_of_fit%p_value > 1.0_dp)) &
      error stop 'rugarch goodness-of-fit test failed'
   duration_test = rugarch_var_duration_test(0.05_dp, simulation%series, &
      value_at_risk)
   if (duration_test%info /= 0 .or. duration_test%shape <= 0.0_dp) &
      error stop 'rugarch VaR duration test failed'
   gmm_test = rugarch_gmm_test(fitted%filtered%standardized_residuals, lags=2)
   if (gmm_test%info /= 0 .or. gmm_test%joint_p_value < 0.0_dp .or. &
      gmm_test%joint_p_value > 1.0_dp) error stop 'rugarch GMM test failed'
   probability_transform = standardized_cdf( &
      fitted%filtered%standardized_residuals, rugarch_distribution_normal, &
      0.0_dp, 1.0_dp, -0.5_dp)
   hong_li_test = rugarch_hong_li_test(probability_transform, lags=2)
   if (hong_li_test%info /= 0 .or. &
      any(.not. ieee_is_finite(hong_li_test%statistic))) &
      error stop 'rugarch Hong-Li test failed'

   print '(a)', 'rugarch tests passed'

contains

   pure function test_parameters(specification) result(parameters)
      !! Construct stable parameters for one supported test model.
      type(rugarch_spec_t), intent(in) :: specification !! Model specification.
      type(rugarch_parameters_t) :: parameters

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
      parameters%alpha = 0.07_dp
      parameters%beta = 0.90_dp
      parameters%asymmetry = 0.0_dp
      parameters%fgarch_shift = 0.0_dp
      parameters%mean_regression = 0.0_dp
      parameters%variance_regression = 0.0_dp
      parameters%omega = 0.03_dp
      parameters%power = 2.0_dp
      parameters%skew = 1.0_dp
      if (specification%variance_model == rugarch_model_igarch) then
         parameters%alpha = 0.08_dp
         parameters%beta = 0.92_dp
         parameters%omega = 0.01_dp
      else if (specification%variance_model == rugarch_model_egarch) then
         parameters%omega = 0.0_dp
         parameters%alpha = 0.12_dp
         parameters%asymmetry = -0.08_dp
         parameters%beta = 0.95_dp
      else if (specification%variance_model == rugarch_model_gjrgarch) then
         parameters%alpha = 0.05_dp
         parameters%asymmetry = 0.08_dp
         parameters%beta = 0.88_dp
         parameters%omega = 0.03_dp
      else if (specification%variance_model == rugarch_model_aparch) then
         parameters%alpha = 0.06_dp
         parameters%asymmetry = 0.15_dp
         parameters%beta = 0.86_dp
         parameters%omega = 0.05_dp
         parameters%power = 1.7_dp
      else if (specification%variance_model == rugarch_model_figarch) then
         parameters%omega = 0.03_dp
         parameters%alpha = 0.10_dp
         parameters%beta = 0.40_dp
         parameters%fractional = 0.40_dp
      else if (specification%variance_model == rugarch_model_csgarch) then
         parameters%omega = 0.03_dp
         parameters%alpha = 0.10_dp
         parameters%beta = 0.60_dp
         parameters%component_rho = 0.95_dp
         parameters%component_phi = 0.08_dp
      else if (specification%variance_model == rugarch_model_realgarch) then
         parameters%omega = -0.10_dp
         parameters%alpha = 0.10_dp
         parameters%beta = 0.80_dp
         parameters%measurement_intercept = 0.0_dp
         parameters%measurement_slope = 1.0_dp
         parameters%measurement_leverage1 = -0.05_dp
         parameters%measurement_leverage2 = 0.02_dp
         parameters%measurement_sigma = 0.30_dp
      else if (specification%variance_model == rugarch_model_fgarch) then
         parameters%omega = 0.03_dp
         parameters%alpha = 0.06_dp
         parameters%beta = 0.85_dp
         parameters%asymmetry = 0.10_dp
         parameters%fgarch_shift = 0.05_dp
         parameters%fgarch_delta = 0.0_dp
         parameters%fgarch_lambda = 2.0_dp
      end if
      if (specification%distribution == rugarch_distribution_student .or. &
         specification%distribution == rugarch_distribution_skew_student) &
         parameters%shape = 8.0_dp
      if (specification%distribution == rugarch_distribution_ged .or. &
         specification%distribution == rugarch_distribution_skew_ged) &
         parameters%shape = 1.5_dp
      if (specification%distribution >= rugarch_distribution_skew_normal .and. &
         specification%distribution <= rugarch_distribution_skew_ged) &
         parameters%skew = 1.35_dp
      if (specification%distribution == rugarch_distribution_johnson_su) then
         parameters%shape = 1.5_dp
         parameters%skew = -0.35_dp
      end if
      if (specification%distribution == rugarch_distribution_nig) then
         parameters%shape = 3.0_dp
         parameters%skew = -0.20_dp
         parameters%lambda = -0.5_dp
      else if (specification%distribution == rugarch_distribution_ghyp) then
         parameters%shape = 3.0_dp
         parameters%skew = 0.20_dp
         parameters%lambda = 1.0_dp
      else if (specification%distribution == &
         rugarch_distribution_gh_skew_student) then
         parameters%shape = 8.0_dp
         parameters%skew = 0.30_dp
      end if
   end function test_parameters

end program test_rugarch
