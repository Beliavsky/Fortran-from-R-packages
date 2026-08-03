! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm
  use multiatsm_kinds, only : dp, pi
  use multiatsm_types, only : var_model, varx_country_model, gvar_model, jll_model, &
    affine_loadings, atsm_likelihood_result, response_result, variance_decomposition_result, &
    forecast_result, bootstrap_result, VARX_UNCONSTRAINED, VARX_SPANNED_RESTRICTED, &
    VARX_FACTOR_RESTRICTED
  use multiatsm_pca, only : pca_weights_one_country, spanned_factors, pca_variance_explained
  use multiatsm_var, only : fit_var, fit_restricted_ols, build_star_factors, fit_varx_system, &
    build_gvar, fit_gvar, transition_matrix_year, transition_matrix_mean
  use multiatsm_jll, only : fit_jll, jll_feedback_restrictions, jll_cholesky_mask
  use multiatsm_affine, only : affine_yield_loadings, multicountry_affine_loadings, &
    pricing_factor_loadings, rotate_latent_to_observed, estimate_long_run_short_rate, &
    build_yield_intercepts
  use multiatsm_likelihood, only : gaussian_log_density, atsm_log_likelihood, yield_error_variance
  use multiatsm_outputs, only : fitted_yields, forecast_yields, forecast_rmse, impulse_responses, &
    generalized_impulse_responses, forecast_error_variance_decomposition, generalized_fevd, &
    expected_short_rate_component, term_premium, forward_rates
  use multiatsm_optimization, only : optimization_result, numerical_gradient, numerical_jacobian, &
    bfgs_minimize, nelder_mead_minimize, stabilize_transition, lower_factor_to_psd, &
    psd_to_lower_parameters, block_diagonal_psd, scale_from_jacobian
  use multiatsm_bootstrap, only : BOOTSTRAP_IID, BOOTSTRAP_WILD, BOOTSTRAP_BLOCK, &
    resample_residuals, simulate_var, bootstrap_var, percentile_bounds
  use multiatsm_bias, only : bias_correct_var, shrink_transition
  use multiatsm_random, only : set_random_seed, random_normal, random_integer
  implicit none
  public
end module multiatsm
