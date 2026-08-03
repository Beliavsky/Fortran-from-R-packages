! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module sandwich
   use sandwich_kinds, only : dp
   use sandwich_status, only : SANDWICH_SUCCESS, SANDWICH_INVALID_ARGUMENT, &
      SANDWICH_DIMENSION_MISMATCH, SANDWICH_SINGULAR_MATRIX, &
      SANDWICH_INSUFFICIENT_DATA, SANDWICH_NUMERICAL_FAILURE, SANDWICH_UNSUPPORTED
   use sandwich_linalg, only : identity_matrix, inverse_matrix, solve_linear, &
      covariance_matrix, symmetric_eigen_jacobi, symmetric_matrix_power, project_psd
   use sandwich_regression, only : ols_model, fit_ols, ols_scores, ols_bread, ols_hatvalues
   use sandwich_core, only : meat, sandwich_covariance, vcov_opg, bread_from_information
   use sandwich_hc, only : hc_weights, meat_hc, vcov_hc
   use sandwich_kernels, only : kernel_weight, kernel_weights
   use sandwich_auxiliary, only : pava_result, pava_blocks, pava_fitted, autocorrelation, isoacf
   use sandwich_hac, only : hac_diagnostics, meat_hac, vcov_hac, bandwidth_andrews, &
      andrews_weights, bandwidth_newey_west, newey_west_weights, lumley_weights, &
      long_run_variance, prewhite_var
   use sandwich_cluster, only : meat_cluster, vcov_cluster
   use sandwich_panel, only : PANEL_LAG_MAX, PANEL_LAG_NW1987, PANEL_LAG_NW1994, &
      meat_panel_longitudinal, vcov_panel_longitudinal, meat_panel_corrected, &
      vcov_panel_corrected
   use sandwich_bootstrap, only : bootstrap_covariance, jackknife_covariance, &
      vcov_bootstrap_ols, set_bootstrap_seed
   implicit none
   public
end module sandwich
