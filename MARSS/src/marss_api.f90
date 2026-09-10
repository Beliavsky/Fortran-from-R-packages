! SPDX-License-Identifier: GPL-2.0-only
module marss_api
   use marss_kinds, only : dp, i8
   use marss_types, only : marss_model, marss_model_spec, marss_dfa_spec, marss_kf_result, marss_fit_result
   use marss_types, only : marss_hessian_result
   use marss_types, only : marss_constraint_block, marss_constraints
   use marss_types, only : marss_simulation, marss_cv_result, marss_innov_boot_result, marss_residual_result, marss_hatyt_result
   use marss_utils, only : zscore, ldiag
   use marss_model_ops, only : marss_inits, marss_model_valid
   use marss_kemcheck_mod, only : marss_kemcheck
   use marss_kalman, only : marss_kf, marss_kfss
   use marss_kfas_bridge, only : marss_kfas, marss_kfas_loglik, marss_to_kfas_model
   use marss_simulation_mod, only : marss_simulate
   use marss_em, only : marss_kem, marss_fit, marss
   use marss_analysis, only : marss_aic, marss_hatyt, marss_residuals
   use marss_analysis, only : marss_vectorizeparam, marss_unvectorizeparam
   use marss_analysis, only : marss_hessian, marss_fisher_i, marss_param_cis
   use marss_analysis, only : marss_parameter_count
   use marss_bootstrap, only : marss_boot, marss_boot_hessian, marss_bootstrap_aic, marss_bootstrap_param_cis
   use marss_innovations, only : marss_innovations_boot
   use marss_optim_mod, only : marss_optim
   use marss_cv_mod, only : marss_cv
   use marss_constraints_mod, only : marss_constraints_valid, marss_constraints_parameter_count
   use marss_constraints_mod, only : marss_constraint_from_labels, marss_constraint_from_affine, &
      marss_constraint_from_entries
   use marss_constraints_mod, only : marss_free_parameter_names, marss_vectorized_parameter_names
   use marss_constraints_mod, only : marss_reorder_free_parameters, marss_constraints_set_start_named, &
      marss_constraints_update_start_named
   use marss_constraints_mod, only : marss_vectorize_free, marss_unvectorize_free
   use marss_constraints_mod, only : marss_constraint_start_vector, marss_constraints_set_start
   use marss_constraints_mod, only : marss_apply_constraints, marss_optim_linear
   use marss_constraints_mod, only : marss_hessian_linear, marss_fisher_i_linear, marss_param_cis_linear
   use marss_constrained_em, only : marss_kem_linear
   use marss_residuals_full, only : marss_residuals_smoothed, marss_residuals_filtered, marss_residuals_one_step
   use marss_residuals_full, only : marss_residuals_harvey
   use marss_harvey_info, only : marss_fisher_i_harvey_linear, marss_param_cis_harvey_linear
   use marss_hessian_summary_mod, only : marss_hessian_summary_linear
   use marss_dfa_mod, only : marss_dfa_build, marss_dfa_model, marss_dfa_fit, marss_dfa_fit_spec
   use marss_builder, only : marss_build
   use marss_workflow, only : marss_from_data
   use marss_initialization, only : marss_inits_linear, marss_inits_named
   implicit none
   private
   public :: dp
   public :: i8
   public :: marss_model
   public :: marss_model_spec
   public :: marss_dfa_spec
   public :: marss_kf_result
   public :: marss_fit_result
   public :: marss_hessian_result
   public :: marss_constraint_block
   public :: marss_constraints
   public :: marss_simulation
   public :: marss_cv_result
   public :: marss_innov_boot_result
   public :: marss_residual_result
   public :: marss_hatyt_result
   public :: zscore
   public :: ldiag
   public :: marss_inits
   public :: marss_inits_linear
   public :: marss_inits_named
   public :: marss_kemcheck
   public :: marss_model_valid
   public :: marss_kf
   public :: marss_kfss
   public :: marss_kfas
   public :: marss_kfas_loglik
   public :: marss_to_kfas_model
   public :: marss_simulate
   public :: marss_kem
   public :: marss_fit
   public :: marss
   public :: marss_aic
   public :: marss_hatyt
   public :: marss_residuals
   public :: marss_vectorizeparam
   public :: marss_unvectorizeparam
   public :: marss_hessian
   public :: marss_fisher_i
   public :: marss_param_cis
   public :: marss_parameter_count
   public :: marss_boot
   public :: marss_boot_hessian
   public :: marss_bootstrap_aic
   public :: marss_bootstrap_param_cis
   public :: marss_innovations_boot
   public :: marss_optim
   public :: marss_cv
   public :: marss_constraints_valid
   public :: marss_constraints_parameter_count
   public :: marss_constraint_from_labels
   public :: marss_constraint_from_affine
   public :: marss_constraint_from_entries
   public :: marss_free_parameter_names
   public :: marss_vectorized_parameter_names
   public :: marss_reorder_free_parameters
   public :: marss_constraints_set_start_named
   public :: marss_constraints_update_start_named
   public :: marss_vectorize_free
   public :: marss_unvectorize_free
   public :: marss_constraint_start_vector
   public :: marss_constraints_set_start
   public :: marss_apply_constraints
   public :: marss_optim_linear
   public :: marss_hessian_linear
   public :: marss_fisher_i_linear
   public :: marss_param_cis_linear
   public :: marss_kem_linear
   public :: marss_residuals_smoothed
   public :: marss_residuals_filtered
   public :: marss_residuals_one_step
   public :: marss_residuals_harvey
   public :: marss_fisher_i_harvey_linear
   public :: marss_param_cis_harvey_linear
   public :: marss_hessian_summary_linear
   public :: marss_dfa_build
   public :: marss_dfa_model
   public :: marss_dfa_fit
   public :: marss_dfa_fit_spec
   public :: marss_build
   public :: marss_from_data
end module marss_api
