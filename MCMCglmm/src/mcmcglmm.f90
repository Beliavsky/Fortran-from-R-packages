! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm
    use r_kinds, only : dp
    use mcmcglmm_rng, only : rng_state, rng_seed, rng_uniform, rng_normal, rng_exponential, rng_gamma, rng_chisq, rng_poisson
    use mcmcglmm_math, only : normal_logpdf, normal_cdf, normal_quantile, logistic, log1pexp, logsumexp, &
        poisson_logpmf, &
        geometric_logpmf, exponential_logpdf, weibull_logpdf, family_loglik
    use mcmcglmm_matrix, only : covariance_to_correlation, block_diagonal, commutation_matrix, kronecker_product, &
        mvn_log_density, sample_mvn_covariance, sample_mvn_precision, symmetric_eigenvalues
    use mcmcglmm_distributions, only : truncated_normal_sample, conditional_mvn_parameters, conditional_mvn_log_density, &
        truncated_conditional_mvn_sample, truncated_conditional_mvn_log_probability, inverse_wishart_sample, &
        riw_mcmcglmm, riw_mcmcglmm_conditioned, pkk_probability
    use mcmcglmm_covariance, only : conditioned_covariance_update, correlation_structure_update, &
        correlation_submatrix_update, identity_direct_sum_update, covariance_update_dispatch
    use mcmcglmm_pedigree, only : pedigree_relationship, pedigree_inverse, prune_pedigree_mask, breeding_values_pedigree
    use mcmcglmm_posterior, only : posterior_correlations, posterior_inverses, posterior_eigenvalues, &
        posterior_modes, ante_parameters
    use mcmcglmm_utilities, only : triangle_to_matrix, matrix_to_triangle, uniform_central_moment, central_moment_tensor, &
        krzanowski_compare, symmetrizer_matrix, normal_moment_matrix
    use mcmcglmm_sampler, only : gaussian_mcmc_result, ordinal_mcmc_result, gaussian_mixed_mcmc, &
        ordinal_probit_mixed_mcmc
    use mcmcglmm_families, only : binomial_logpmf, poisson_cdf, student_t_logpdf, noncentral_t_logpdf, &
        ordinal_probit_loglik, scalar_family_loglik, two_part_family_loglik, multinomial_log_kernel, &
        ztmb_log_kernel, ztmultinomial_log_kernel
    use mcmcglmm_family_sampler, only : family_mcmc_result, multivariate_family_mcmc_result, &
        threshold_cutpoint_mcmc_result, ordinal_native_mcmc_result, latent_family_mixed_mcmc, two_part_mixed_mcmc, &
        multinomial_family_mixed_mcmc, heterogeneous_family_mixed_mcmc, threshold_cutpoint_mixed_mcmc, &
        ordinal_native_mixed_mcmc
    use mcmcglmm_ante, only : ante_covariance_samples
    use mcmcglmm_phylo, only : phylogenetic_precision, breeding_values_phylo
    use mcmcglmm_design, only : path_matrix, sir_matrix, multiple_membership_design, gelman_prior_design, &
        random_effect_covariance, d_divergence_mc
    use mcmcglmm_spline, only : spline_lrtp
    use mcmcglmm_parameter_expansion, only : parameter_expansion_conditional, apply_parameter_expansion, &
        expanded_covariance
    use mcmcglmm_engine_features, only : optimal_acceptance_ratio, adaptive_mh_observe, adaptive_mh_decay, &
        adaptive_mh_finalize, binary_slice_liability_update, theta_scale_conditional, &
        categorical_measurement_error_update, structural_matrix, structural_transform, structural_gaussian_loglik, &
        structural_parameter_random_walk_update
    use mcmcglmm_theta_sampler, only : theta_scale_mcmc_result, theta_scale_multivariate_conditional, &
        theta_scale_gaussian_mixed_mcmc
    use mcmcglmm_structural_sampler, only : structural_gaussian_mcmc_result, structural_gaussian_multi_term_mcmc
    use mcmcglmm_px_sampler, only : gaussian_px_mcmc_result, gaussian_parameter_expanded_mcmc, &
        px_alpha_conditional
    use mcmcglmm_multiterm_sampler, only : multi_term_gaussian_mcmc_result, multi_term_gaussian_mixed_mcmc, &
        multi_term_coefficient_conditional, multi_term_coefficient_conditional_sparse
    use mcmcglmm_unified_sampler, only : unified_family_mcmc_result, heterogeneous_multi_term_mixed_mcmc
    use mcmcglmm_grouped_multiterm_sampler, only : grouped_multi_term_mcmc_result, &
        two_part_multi_term_mixed_mcmc, multinomial_multi_term_mixed_mcmc
    use mcmcglmm_joint_gr, only : joint_gr_decompose, joint_gr_compose, joint_gr_covariance_update
    use mcmcglmm_covu_sampler, only : covu_gaussian_mcmc_result, covu_coefficient_conditional, &
        covu_gaussian_loglik, covu_gaussian_mixed_mcmc
    use mcmcglmm_multiterm_px_sampler, only : multi_term_gaussian_px_mcmc_result, multi_term_family_px_mcmc_result, &
        multi_term_px_alpha_conditional, multi_term_gaussian_parameter_expanded_mcmc, &
        heterogeneous_multi_term_parameter_expanded_mcmc
    use mcmcglmm_simulation, only : simulate_scalar_response, simulate_ordinal_response, &
        simulate_threshold_response, simulate_two_part_response, simulate_multinomial_response, &
        simulate_multi_term_gaussian_latent
    use mcmcglmm_prediction, only : multi_term_build_v, posterior_linear_predictor, scalar_response_expectation
    use mcmcglmm_orchestrator, only : mcmcglmm_engine_scalar, mcmcglmm_engine_scalar_px, &
        mcmcglmm_engine_two_process, mcmcglmm_engine_multinomial, mcmcglmm_engine_ordinal, &
        mcmcglmm_engine_threshold, mcmcglmm_engine_theta_scale, mcmcglmm_engine_structural, &
        mcmcglmm_engine_covu, mcmcglmm_orchestrator_invalid_engine, &
        mcmcglmm_orchestrator_missing_input, mcmcglmm_orchestrator_invalid_control, mcmcglmm_numeric_model, &
        mcmcglmm_numeric_prior, mcmcglmm_control, mcmcglmm_numeric_result, mcmcglmm_fit_numeric, &
        mcmcglmm_validate_numeric
    use mcmcglmm_sparse, only : mcmcglmm_sparse_matrix, sparse_crossproduct, sparse_from_coo, &
        sparse_from_dense, sparse_is_initialized, sparse_matmul_matrix, sparse_matmul_vector, sparse_to_dense, &
        sparse_stacked_crossproduct, sparse_transpose_matmul_matrix, sparse_transpose_matmul_sparse, sparse_validate
    use mcmcglmm_sparse_factorization, only : sample_mvn_sparse_precision, sparse_cholesky_analysis, &
        sparse_cholesky_analyze, sparse_cholesky_factor, sparse_cholesky_factor_analyzed, sparse_cholesky_solve, &
        sparse_cholesky_transpose_solve, sparse_precision_cache, sparse_reverse_cuthill_mckee, &
        sparse_symmetric_permute
    implicit none
    private

    public :: dp
    public :: rng_state
    public :: rng_seed
    public :: rng_uniform
    public :: rng_normal
    public :: rng_exponential
    public :: rng_gamma
    public :: rng_chisq
    public :: rng_poisson
    public :: normal_logpdf
    public :: normal_cdf
    public :: normal_quantile
    public :: logistic
    public :: log1pexp
    public :: logsumexp
    public :: poisson_logpmf
    public :: geometric_logpmf
    public :: exponential_logpdf
    public :: weibull_logpdf
    public :: family_loglik
    public :: covariance_to_correlation
    public :: block_diagonal
    public :: commutation_matrix
    public :: kronecker_product
    public :: mvn_log_density
    public :: sample_mvn_covariance
    public :: sample_mvn_precision
    public :: symmetric_eigenvalues
    public :: truncated_normal_sample
    public :: conditional_mvn_parameters
    public :: conditional_mvn_log_density
    public :: truncated_conditional_mvn_sample
    public :: truncated_conditional_mvn_log_probability
    public :: inverse_wishart_sample
    public :: riw_mcmcglmm
    public :: riw_mcmcglmm_conditioned
    public :: pkk_probability
    public :: conditioned_covariance_update
    public :: correlation_structure_update
    public :: correlation_submatrix_update
    public :: identity_direct_sum_update
    public :: covariance_update_dispatch
    public :: pedigree_relationship
    public :: pedigree_inverse
    public :: prune_pedigree_mask
    public :: breeding_values_pedigree
    public :: posterior_correlations
    public :: posterior_inverses
    public :: posterior_eigenvalues
    public :: posterior_modes
    public :: ante_parameters
    public :: triangle_to_matrix
    public :: matrix_to_triangle
    public :: uniform_central_moment
    public :: central_moment_tensor
    public :: krzanowski_compare
    public :: symmetrizer_matrix
    public :: normal_moment_matrix
    public :: gaussian_mcmc_result
    public :: ordinal_mcmc_result
    public :: gaussian_mixed_mcmc
    public :: ordinal_probit_mixed_mcmc
    public :: binomial_logpmf
    public :: poisson_cdf
    public :: student_t_logpdf
    public :: noncentral_t_logpdf
    public :: ordinal_probit_loglik
    public :: scalar_family_loglik
    public :: two_part_family_loglik
    public :: multinomial_log_kernel
    public :: ztmb_log_kernel
    public :: ztmultinomial_log_kernel
    public :: family_mcmc_result
    public :: multivariate_family_mcmc_result
    public :: threshold_cutpoint_mcmc_result
    public :: ordinal_native_mcmc_result
    public :: latent_family_mixed_mcmc
    public :: two_part_mixed_mcmc
    public :: multinomial_family_mixed_mcmc
    public :: heterogeneous_family_mixed_mcmc
    public :: threshold_cutpoint_mixed_mcmc
    public :: ordinal_native_mixed_mcmc
    public :: ante_covariance_samples
    public :: phylogenetic_precision
    public :: breeding_values_phylo
    public :: path_matrix
    public :: sir_matrix
    public :: multiple_membership_design
    public :: gelman_prior_design
    public :: random_effect_covariance
    public :: d_divergence_mc
    public :: spline_lrtp
    public :: parameter_expansion_conditional
    public :: apply_parameter_expansion
    public :: expanded_covariance
    public :: optimal_acceptance_ratio
    public :: adaptive_mh_observe
    public :: adaptive_mh_decay
    public :: adaptive_mh_finalize
    public :: binary_slice_liability_update
    public :: theta_scale_conditional
    public :: categorical_measurement_error_update
    public :: structural_matrix
    public :: structural_transform
    public :: structural_gaussian_loglik
    public :: structural_parameter_random_walk_update
    public :: theta_scale_mcmc_result
    public :: theta_scale_multivariate_conditional
    public :: theta_scale_gaussian_mixed_mcmc
    public :: structural_gaussian_mcmc_result
    public :: structural_gaussian_multi_term_mcmc
    public :: gaussian_px_mcmc_result
    public :: gaussian_parameter_expanded_mcmc
    public :: px_alpha_conditional
    public :: multi_term_gaussian_mcmc_result
    public :: multi_term_gaussian_mixed_mcmc
    public :: multi_term_coefficient_conditional
    public :: multi_term_coefficient_conditional_sparse
    public :: unified_family_mcmc_result
    public :: heterogeneous_multi_term_mixed_mcmc
    public :: grouped_multi_term_mcmc_result
    public :: two_part_multi_term_mixed_mcmc
    public :: multinomial_multi_term_mixed_mcmc
    public :: joint_gr_decompose
    public :: joint_gr_compose
    public :: joint_gr_covariance_update
    public :: covu_gaussian_mcmc_result
    public :: covu_coefficient_conditional
    public :: covu_gaussian_loglik
    public :: covu_gaussian_mixed_mcmc
    public :: multi_term_gaussian_px_mcmc_result
    public :: multi_term_family_px_mcmc_result
    public :: multi_term_px_alpha_conditional
    public :: multi_term_gaussian_parameter_expanded_mcmc
    public :: heterogeneous_multi_term_parameter_expanded_mcmc
    public :: simulate_scalar_response
    public :: simulate_ordinal_response
    public :: simulate_threshold_response
    public :: simulate_two_part_response
    public :: simulate_multinomial_response
    public :: simulate_multi_term_gaussian_latent
    public :: multi_term_build_v
    public :: posterior_linear_predictor
    public :: scalar_response_expectation
    public :: mcmcglmm_engine_scalar
    public :: mcmcglmm_engine_scalar_px
    public :: mcmcglmm_engine_two_process
    public :: mcmcglmm_engine_multinomial
    public :: mcmcglmm_engine_ordinal
    public :: mcmcglmm_engine_threshold
    public :: mcmcglmm_engine_theta_scale
    public :: mcmcglmm_engine_structural
    public :: mcmcglmm_engine_covu
    public :: mcmcglmm_orchestrator_invalid_engine
    public :: mcmcglmm_orchestrator_missing_input
    public :: mcmcglmm_orchestrator_invalid_control
    public :: mcmcglmm_numeric_model
    public :: mcmcglmm_numeric_prior
    public :: mcmcglmm_control
    public :: mcmcglmm_numeric_result
    public :: mcmcglmm_fit_numeric
    public :: mcmcglmm_validate_numeric
    public :: mcmcglmm_sparse_matrix
    public :: sparse_crossproduct
    public :: sparse_from_coo
    public :: sparse_from_dense
    public :: sparse_is_initialized
    public :: sparse_matmul_matrix
    public :: sparse_matmul_vector
    public :: sparse_stacked_crossproduct
    public :: sparse_to_dense
    public :: sparse_transpose_matmul_matrix
    public :: sparse_transpose_matmul_sparse
    public :: sparse_validate
    public :: sample_mvn_sparse_precision
    public :: sparse_cholesky_analysis
    public :: sparse_cholesky_analyze
    public :: sparse_cholesky_factor
    public :: sparse_cholesky_factor_analyzed
    public :: sparse_cholesky_solve
    public :: sparse_cholesky_transpose_solve
    public :: sparse_reverse_cuthill_mckee
    public :: sparse_symmetric_permute
    public :: sparse_precision_cache

end module mcmcglmm
