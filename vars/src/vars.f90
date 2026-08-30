! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars
   use r_kinds, only : dp
   use vars_types
   use vars_regression, only : fit_var, var_select, restrict_var_manual, restrict_var_ser
   use vars_regression, only : var_loglik, rebuild_var_statistics
   use vars_dynamics, only : phi_from_a, psi_from_a_sigma, var_roots
   use vars_dynamics, only : forecast_covariance, forecast_var, impulse_response, fevd_var
   use vars_dynamics, only : structural_impulse_response, structural_fevd
   use vars_dynamics, only : bq_identification, vec2var_coefficients, residual_covariance_unbiased
   use vars_diagnostics, only : arch_test_univariate, arch_test_multivariate
   use vars_diagnostics, only : jarque_bera_univariate, jarque_bera_multivariate
   use vars_diagnostics, only : portmanteau_tests, bg_serial_tests
   use vars_causality, only : granger_causality, instantaneous_causality
   use vars_structural, only : svar_negloglik, svar_fit_scoring, structural_impact
   use vars_structural, only : svec_long_run_matrix, svec_fit_scoring
   use vars_bootstrap, only : residual_bootstrap_path, bootstrap_irf_indices
   implicit none
   public
   public :: dp
end module vars
