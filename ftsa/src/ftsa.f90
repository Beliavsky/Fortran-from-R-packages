! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa
   use r_kinds, only : dp
   use ftsa_types, only : ftsm_fit, hurst_result
   use ftsa_statistics, only : ftsa_ok, ftsa_invalid_input, ftsa_unsupported
   use ftsa_statistics, only : ftsa_extract_time, ftsa_extract_x, ftsa_residuals, ftsa_difference
   use ftsa_statistics, only : ftsa_mean_coordinate, ftsa_median_coordinate, ftsa_variance_coordinate
   use ftsa_statistics, only : ftsa_sd_coordinate, ftsa_quantiles, ftsa_sd_vector, ftsa_sd_columns
   use ftsa_statistics, only : ftsa_centre_mean, ftsa_centre_variance
   use ftsa_statistics, only : ftsa_error_metric, ftsa_error_all, ftsa_bootstrap_mean
   use ftsa_covariance, only : ftsa_facf, ftsa_long_run_covariance
   use ftsa_factor, only : ftsa_er_gr, ftsa_h_estimate_rs
   use ftsa_simulation, only : ftsa_sim_farma, ftsa_stop_time_sim_data
   use ftsa_models, only : ftsa_ftsm, ftsa_forecast_ftsm, ftsa_ftsm_iterative_forecasts
   use ftsa_extended, only : ftsa_dynamic_flr, ftsa_dynupdate, ftsa_forecast_fplsr, ftsa_fplsr
   use ftsa_extended, only : ftsa_ftsm_weight_select, ftsa_isfe, ftsa_maf_multivariate
   use ftsa_extended, only : ftsa_pcscore_bootstrap_data, ftsa_simpls, ftsa_stop_time_detect
   use ftsa_remaining, only : ftsa_coda_bayes_nw, ftsa_coda_fpca, ftsa_farforecast
   use ftsa_remaining, only : ftsa_horta_ziegelmann_fpca, ftsa_lqdt_fpca, ftsa_mfdm, ftsa_skew_t_fun, ftsa_t_stationary
   use ftsa_gaev, only : ftsa_gaev_forecast
   implicit none
   private

   public :: dp, ftsm_fit, hurst_result
   public :: ftsa_ok, ftsa_invalid_input, ftsa_unsupported
   public :: ftsa_extract_time, ftsa_extract_x, ftsa_residuals, ftsa_difference
   public :: ftsa_mean_coordinate, ftsa_median_coordinate, ftsa_variance_coordinate
   public :: ftsa_sd_coordinate, ftsa_quantiles, ftsa_sd_vector, ftsa_sd_columns
   public :: ftsa_centre_mean, ftsa_centre_variance
   public :: ftsa_error_metric, ftsa_error_all, ftsa_bootstrap_mean
   public :: ftsa_facf, ftsa_long_run_covariance, ftsa_er_gr, ftsa_h_estimate_rs
   public :: ftsa_sim_farma, ftsa_stop_time_sim_data
   public :: ftsa_ftsm, ftsa_forecast_ftsm, ftsa_ftsm_iterative_forecasts
   public :: ftsa_dynamic_flr, ftsa_dynupdate, ftsa_forecast_fplsr, ftsa_fplsr
   public :: ftsa_ftsm_weight_select, ftsa_isfe, ftsa_maf_multivariate
   public :: ftsa_pcscore_bootstrap_data, ftsa_simpls, ftsa_stop_time_detect
   public :: ftsa_coda_bayes_nw, ftsa_coda_fpca, ftsa_farforecast
   public :: ftsa_horta_ziegelmann_fpca, ftsa_lqdt_fpca, ftsa_mfdm, ftsa_skew_t_fun, ftsa_t_stationary
   public :: ftsa_gaev_forecast

end module ftsa
