! SPDX-License-Identifier: MIT
module bayesianou
  use bayesianou_kinds
  use bayesianou_types
  use bayesianou_utils, only : zscore_train, compute_common_factor, orthogonalize_series, &
                               ou_level_spec, weighted_com_statistics, align_columns_indices
  use bayesianou_model, only : fit_nested_core => fit_ou_nested, &
                               fit_single_core => fit_ou_nonlinear_tmg, &
                               simulate_ou_nested, ou_log_likelihood, ou_mean_increment, &
                               extract_posterior_summary, build_beta_tmg_table, &
                               summarize_sv_sigmas, drift_decomposition_grid, &
                               build_accounting_block, extract_mu_trajectory
  use bayesianou_diagnostics, only : compute_fit_diagnostics, evaluate_oos, &
                                     evaluate_oos_nested, psis_loo, compare_models_loo, &
                                     kappa_stability_evidence, extract_convergence_evidence, &
                                     count_divergences, validate_ou_fit
  use bayesianou_mi, only : fit_ou_nested_mi, rubin_combine, pack_pooled_draws
  use bayesianou_geometry, only : ou_geom_target_type, ou_geom_metric_type, &
                                  ou_geom_hmc_result, ou_geom_target, &
                                  ou_geom_metric_euclidean, ou_geom_metric_riemannian, &
                                  ou_geom_hmc, ou_geom_mass, ou_geom_ebfmi
  implicit none
  private

  public :: dp, pi, status_ok, status_bad_input, status_not_converged, status_singular
  public :: ou_input, ou_options, ou_priors, ou_summary, ou_draws, ou_fit_result
  public :: ou_level_flags, ou_level_spec_type, zscore_result, loo_result, oos_metric
  public :: ou_diagnostics, stability_result, rubin_result, ou_mi_result
  public :: level_canonical, level_both_full, level_both_lean, level_n1_lean
  public :: zscore_train, compute_common_factor, orthogonalize_series, ou_level_spec
  public :: weighted_com_statistics, align_columns_indices
  public :: fit_ou_nested, fit_ou_nonlinear_tmg, fit_ou_nested_mi, simulate_ou_nested
  public :: ou_log_likelihood, ou_mean_increment, extract_posterior_summary
  public :: build_beta_tmg_table, summarize_sv_sigmas, drift_decomposition_grid
  public :: build_accounting_block, extract_mu_trajectory
  public :: compute_fit_diagnostics, evaluate_oos, evaluate_oos_nested, psis_loo
  public :: compare_models_loo, kappa_stability_evidence, extract_convergence_evidence
  public :: count_divergences, validate_ou_fit, rubin_combine, pack_pooled_draws
  public :: ou_geom_target_type, ou_geom_metric_type, ou_geom_hmc_result
  public :: ou_geom_target, ou_geom_metric_euclidean, ou_geom_metric_riemannian
  public :: ou_geom_hmc, ou_geom_mass, ou_geom_ebfmi

contains

  subroutine fit_ou_nested(input,options,result)
    type(ou_input),intent(in),target::input
    type(ou_options),intent(in)::options
    type(ou_fit_result),intent(out),target::result
    call fit_nested_core(input,options,result)
    if(result%status==status_ok)call compute_fit_diagnostics(input,result)
  end subroutine fit_ou_nested

  subroutine fit_ou_nonlinear_tmg(input,options,result)
    type(ou_input),intent(in),target::input
    type(ou_options),intent(in)::options
    type(ou_fit_result),intent(out),target::result
    type(ou_options)::opt
    opt=options;opt%n_levels=1
    call fit_single_core(input,opt,result)
    if(result%status==status_ok)call compute_fit_diagnostics(input,result)
  end subroutine fit_ou_nonlinear_tmg

end module bayesianou
