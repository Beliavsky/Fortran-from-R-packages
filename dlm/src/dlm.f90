! SPDX-License-Identifier: GPL-2.0-or-later
module dlm
    use dlm_types, only : dp
    use dlm_types, only : dlm_model, dlm_filter_result, dlm_smooth_result, dlm_forecast_result
    use dlm_types, only : dlm_gibbs_dig_result
    use dlm_types, only : dlm_success, dlm_invalid_shape, dlm_invalid_model, dlm_linalg_failure, dlm_invalid_argument
    use dlm_types, only : dlm_sampling_failure
    use dlm_types, only : dlm_model_is_valid, dlm_model_is_time_varying, dlm_matrices_at
    use dlm_models, only : ar_trans_pars, block_diag2, dlm_add, dlm_mod_arma, dlm_mod_poly
    use dlm_models, only : dlm_mod_reg, dlm_mod_seas, dlm_mod_trig, dlm_sum
    use dlm_kalman, only : dlm_filter, dlm_forecast, dlm_ll, dlm_residuals, dlm_smooth, dlm_svd2var
    use dlm_random, only : dlm_bsample, dlm_gibbs_dig, dlm_random_model, rwishart, seed_dlm_rng
    use dlm_mcmc, only : erg_mean, mcmc_mean, mcmc_sd
    use dlm_bounds, only : convex_bounds, dlm_indicator
    use dlm_arms, only : arms, dlm_log_density
    use dlm_mle_mod, only : dlm_mle, dlm_mle_result, dlm_model_builder
    implicit none
    private

    public :: dp
    public :: dlm_model, dlm_filter_result, dlm_smooth_result, dlm_forecast_result
    public :: dlm_gibbs_dig_result
    public :: dlm_success, dlm_invalid_shape, dlm_invalid_model, dlm_linalg_failure, dlm_invalid_argument
    public :: dlm_sampling_failure
    public :: dlm_model_is_valid, dlm_model_is_time_varying, dlm_matrices_at
    public :: ar_trans_pars, block_diag2, dlm_add, dlm_mod_arma, dlm_mod_poly
    public :: dlm_mod_reg, dlm_mod_seas, dlm_mod_trig, dlm_sum
    public :: dlm_filter, dlm_forecast, dlm_ll, dlm_residuals, dlm_smooth, dlm_svd2var
    public :: dlm_bsample, dlm_gibbs_dig, dlm_random_model, rwishart, seed_dlm_rng
    public :: erg_mean, mcmc_mean, mcmc_sd
    public :: convex_bounds, dlm_indicator
    public :: arms, dlm_log_density
    public :: dlm_mle, dlm_mle_result, dlm_model_builder

end module dlm
