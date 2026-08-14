! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust
  use mclust_kinds, only : dp
  use mclust_types, only : em_control, mclust_fit, mclust_selection
  use mclust_models, only : fit_model, mstep_model, initialize_responsibilities, &
                            model_supported, default_model_names
  use mclust_selection_mod, only : n_var_params, n_mclust_params, bic_value, icl_value, &
                                   mclust_bic, mclust_select, pick_bic
  use mclust_math, only : mixture_posterior, dmvnorm, logsumexp, softmax, map_z, &
                          adjusted_rand_index, brier_score, count_values, &
                          covariance_weighted, random_orthogonal_matrix
  use mclust_simulation, only : simulate_mixture, simulate_fit, sample_multivariate_normal
  use mclust_density, only : component_log_density, mixture_log_density, cdf_mclust_1d, &
                             quantile_mclust_1d, hdr_levels
  use mclust_classification, only : mclust_da_fit, fit_mclust_da, predict_mclust_da, class_error_rate
  use mclust_impute, only : impute_data
  use mclust_dr, only : mclust_dr_fit, fit_mclust_dr, project_mclust_dr
  use mclust_combine, only : cluster_combination, clust_combi, entropy_z, apply_combination
  use mclust_weighted, only : fit_model_weighted
  use mclust_hierarchical, only : hc_result, hc_fit, hclass, hc_responsibilities
  use mclust_ssc, only : fit_model_ssc, mclust_ssc_select
  use mclust_bootstrap, only : bootstrap_lrt_result, parameter_bootstrap_result, &
                               bootstrap_lrt, mclust_parameter_bootstrap
  use mclust_utilities, only : unmap_classes, majority_vote, match_clusters
  use mclust_crimcoords, only : crimcoords_fit, fit_crimcoords
  implicit none
  private
  public :: dp, em_control, mclust_fit, mclust_selection
  public :: fit_model, mstep_model, estep_model, initialize_responsibilities
  public :: model_supported, default_model_names
  public :: n_var_params, n_mclust_params, bic_value, icl_value, mclust_bic, mclust_select, pick_bic
  public :: dmvnorm, logsumexp, softmax, map_z, adjusted_rand_index, brier_score, count_values
  public :: covariance_weighted, random_orthogonal_matrix
  public :: unmap_classes, majority_vote, match_clusters
  public :: crimcoords_fit, fit_crimcoords
  public :: simulate_mixture, simulate_fit, sample_multivariate_normal
  public :: component_log_density, mixture_log_density, cdf_mclust_1d, quantile_mclust_1d, hdr_levels
  public :: mclust_da_fit, fit_mclust_da, predict_mclust_da, class_error_rate
  public :: impute_data
  public :: mclust_dr_fit, fit_mclust_dr, project_mclust_dr
  public :: cluster_combination, clust_combi, entropy_z, apply_combination
  public :: fit_model_weighted
  public :: hc_result, hc_fit, hclass, hc_responsibilities
  public :: fit_model_ssc, mclust_ssc_select
  public :: bootstrap_lrt_result, parameter_bootstrap_result, bootstrap_lrt, mclust_parameter_bootstrap

contains

  subroutine estep_model(x,fit,z,loglik,log_density,status)
    real(dp),intent(in)::x(:,:)
    type(mclust_fit),intent(in)::fit
    real(dp),allocatable,intent(out)::z(:,:)
    real(dp),intent(out)::loglik
    real(dp),allocatable,intent(out),optional::log_density(:)
    integer,intent(out),optional::status
    real(dp),allocatable::ld(:)
    integer::info
    allocate(z(size(x,1),fit%g),ld(size(x,1)))
    call mixture_posterior(x,fit%pro,fit%mean,fit%sigma,z,ld,info)
    if(info==0) then; loglik=sum(ld); else; loglik=-huge(1.0_dp); end if
    if(present(log_density)) call move_alloc(ld,log_density)
    if(present(status)) status=info
  end subroutine estep_model
end module mclust
