module lme4
   use lme4_kinds, only : dp, i8
   use lme4_types, only : random_term_t, covariance_block_t, lmm_control_t, &
      glmm_control_t, nlmm_control_t, lmm_result_t, glmm_result_t, &
      nlmm_result_t, bootstrap_result_t, profile_result_t, lm_list_result_t, &
      influence_result_t, gh_rule_t, family_binomial, family_poisson, &
      family_gamma, family_inverse_gaussian, family_negative_binomial, &
      covariance_unstructured, covariance_diagonal, &
      covariance_compound_symmetry, covariance_ar1
   use lme4_covariance, only : sdcor2cov, cov2sdcor, matrix_to_lower, &
      lower_to_matrix, relative_factor_to_covariance, &
      covariance_to_relative_factor, covariance_pca, build_random_design, &
      term_covariance_from_eta
   use lme4_quadrature, only : gh_rule, gh_integrate
   use lme4_family, only : family_spec_t, gaussian_identity_family, &
      binomial_probit_family, binomial_cloglog_family, &
      quasipoisson_log_family, logistic_scalar, normal_cdf_scalar, &
      normal_pdf_scalar, inverse_normal_cdf
   use lme4_lmm, only : fit_lmm, predict_lmm, random_effects_for_term
   use lme4_lmm_pls, only : fit_lmm_pls
   use lme4_glmm, only : fit_glmm, fit_glmer_nb, predict_glmm
   use lme4_custom_glmm, only : fit_glmm_custom, predict_glmm_custom
   use lme4_aghq, only : fit_glmm_aghq
   use lme4_aghq_nd, only : fit_glmm_aghq_multidimensional
   use lme4_nlmm, only : nonlinear_mean_function, fit_nlmm, predict_nlmm, simulate_nlmm
   use lme4_simulation, only : simulate_lmm, simulate_glmm, set_random_seed
   use lme4_diagnostics, only : is_singular, re_pca, &
      standardized_residuals_lmm, pearson_residuals_glmm
   use lme4_inference, only : wald_confint_lmm, wald_confint_glmm, &
      profile_confint_lmm_beta, profile_confint_glmm_beta, &
      parametric_bootstrap_lmm, parametric_bootstrap_glmm, &
      bootstrap_percentile_confint, influence_lmm_groups, &
      influence_glmm_groups, likelihood_ratio_test
   use lme4_grouped, only : fit_lm_list, predict_lm_list
   implicit none
   public
end module lme4
