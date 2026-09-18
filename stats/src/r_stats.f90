! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public facade for the modern Fortran translation of R stats functionality.
module r_stats
   use fastcluster_types, only: hclust_result_t => hclust_result
   use r_kinds, only: dp
   use r_stats_diagnostics, only: box_test
   use r_stats_scale_tests, only: fligner_test
   use r_stats_multiple_testing, only: p_adjust
   use r_stats_anova, only: anova_lm, one_way_anova, oneway_test, oneway_test_result_t
   use r_stats_additional_tests, only: bartlett_test, binom_test, friedman_test, &
                                       mcnemar_test, poisson_test, var_test
   use r_stats_arima, only: arima_css, arima_ml, predict_arima, predict_arima_interval
   use r_stats_autoregression, only: ar, ar_burg, ar_mle, ar_ols, ar_yw, spec_ar
   use r_stats_clustering, only: cutree, dist, dist_matrix, hclust, hclust_condensed, &
                                 kmeans, kmeans_matrix, kmeans_vector
   use r_stats_decomposition, only: decompose
   use r_stats_descriptive, only: cor, cov, cov2cor, iqr, mad, median, quantile, sd, var, &
                                  weighted_mean
   use r_stats_empirical, only: bw_nrd, bw_nrd0, density, density_kernel_epanechnikov, &
                                density_kernel_gaussian, density_kernel_rectangular, &
                                density_kernel_triangular, ecdf, histogram, predict_ecdf
   use r_stats_glm, only: glm_binomial_fit, glm_family_binomial, glm_family_gamma, &
                          glm_family_gaussian, glm_family_poisson, glm_family_quasibinomial, &
                          glm_family_quasipoisson, glm_fit, glm_gamma_fit, glm_gaussian_fit, &
                          glm_inverse_link, glm_link_cauchit, glm_link_cloglog, glm_link_identity, &
                          glm_link_inverse, glm_link_log, glm_link_logit, glm_link_probit, &
                          glm_link_sqrt, glm_link_value, glm_mu_eta, glm_pearson_resid, &
                          glm_poisson_fit, glm_predict_response, glm_quasibinomial_fit, &
                          glm_quasipoisson_fit
   use r_stats_holt_winters, only: holt_winters, holt_winters_auto, holt_winters_optimize, &
                                   holt_winters_reduced, predict_holt_winters, &
                                   predict_holt_winters_interval
   use r_stats_interpolation, only: approx, approxfun, interpolation_constant, &
                                    interpolation_linear, interpolation_rule_constant, &
                                    interpolation_rule_missing, isoreg, predict_approxfun
   use r_stats_manova, only: manova_lm, mlm_fit
   use r_stats_multivariate, only: cancor, mahalanobis, prcomp, scale, scale_columns, &
                                   scale_matrix, scale_vector
   use r_stats_nls, only: nls_basis_interface, nls_confint_profile, nls_fit, &
                          nls_fit_plinear, nls_jacobian_interface, nls_max_iterations, &
                          nls_model_interface
   use r_stats_numerical, only: integrate, numerical_invalid_input, &
                                numerical_iteration_limit, optimize, &
                                root_extend_both, root_extend_lower, root_extend_none, &
                                root_extend_upper, scalar_function_interface, uniroot
   use r_stats_regression, only: lm_aic, lm_coef, lm_confint, lm_cooks_distance, lm_fit, &
                                 lm_fit_general, lm_predict_general, lm_r_squared_general, predict_lm
   use r_stats_rank_tests, only: fisher_test, kruskal_test, ks_test, ks_test_normal, &
                                 ks_test_two_sample, wilcox_test, wilcox_test_two_sample
   use r_stats_self_start, only: ss_asymp, ss_asymp_jacobian, ss_asymp_model, &
                                 ss_gompertz, ss_gompertz_jacobian, ss_gompertz_model, &
                                 ss_logis, ss_logis_jacobian, ss_logis_model, &
                                 ss_micmen, ss_micmen_jacobian, ss_micmen_model, &
                                 ss_weibull, ss_weibull_jacobian, ss_weibull_model
   use r_stats_smoothing, only: ksmooth, loess_fit, lowess, predict_loess, runmed
   use r_stats_smoothing_spline, only: predict_smooth_spline, smooth_spline
   use r_stats_spline_interpolation, only: predict_splinefun, spline, spline_method_fmm, &
                                           spline_method_natural, splinefun
   use r_stats_spectral, only: spec_pgram, spectrum
   use r_stats_stl, only: stl
   use r_stats_time_series, only: acf, arima_sim, arma_acf, ccf, filter, filter_linear, &
                                  filter_recursive, pacf
   use r_stats_tests, only: chisq_test, cor_test, prop_test, t_test, t_test_one, &
                            t_test_p_value, t_test_p_value_one, t_test_p_value_two, t_test_two
   use r_stats_types, only: acf_result_t, anova_comparison_t, ar_fit_t, arima_fit_t, &
                            arima_forecast_t, cancor_result_t, &
                            chisq_test_result_t, cor_test_result_t, &
                            cubic_spline_t, density_result_t, ecdf_t, fisher_test_result_t, &
                            exact_count_test_result_t, variance_test_result_t, &
                            glm_fit_t, holt_winters_fit_t, holt_winters_forecast_t, kmeans_result_t, &
                            histogram_t, interpolation_t, isoreg_t, &
                            integrate_result_t, optimize_result_t, uniroot_result_t, &
                            kruskal_test_result_t, &
                            ks_test_result_t, lm_fit_t, loess_fit_t, nls_fit_t, &
                            manova_result_t, mlm_fit_t, multivariate_spectrum_result_t, &
                            one_way_anova_t, prcomp_fit_t, &
                            prop_test_result_t, scale_result_t, &
                            seasonal_decomposition_t, smooth_spline_fit_t, smooth_xy_t, &
                            spectrum_result_t, stl_result_t, &
                            t_test_result_t, &
                            wilcox_test_result_t
   implicit none
   private

   public :: p_adjust
   public :: box_test
   public :: oneway_test, oneway_test_result_t
   public :: fligner_test
   public :: mcnemar_test

   public :: acf, acf_result_t, anova_comparison_t, anova_lm
   public :: ar, ar_burg, ar_fit_t, ar_mle, ar_ols, ar_yw
   public :: arima_css, arima_fit_t, arima_forecast_t, arima_ml, arima_sim, arma_acf
   public :: ccf, predict_arima, predict_arima_interval
   public :: cancor, cancor_result_t
   public :: approx, approxfun
   public :: bw_nrd, bw_nrd0
   public :: bartlett_test, binom_test
   public :: chisq_test, chisq_test_result_t
   public :: cor_test, cor_test_result_t
   public :: cor, cov, cov2cor
   public :: cubic_spline_t
   public :: density, density_result_t, ecdf, ecdf_t
   public :: density_kernel_epanechnikov, density_kernel_gaussian
   public :: density_kernel_rectangular, density_kernel_triangular
   public :: cutree, decompose, dist, dist_matrix, dp, filter, filter_linear, filter_recursive, pacf
   public :: fisher_test, fisher_test_result_t, glm_binomial_fit, glm_fit, glm_fit_t, glm_gamma_fit
   public :: exact_count_test_result_t, friedman_test
   public :: glm_gaussian_fit, glm_pearson_resid, glm_poisson_fit
   public :: glm_quasibinomial_fit, glm_quasipoisson_fit
   public :: glm_family_binomial, glm_family_gamma, glm_family_gaussian, glm_family_poisson
   public :: glm_family_quasibinomial, glm_family_quasipoisson
   public :: glm_inverse_link, glm_link_value, glm_mu_eta
   public :: glm_link_cauchit, glm_link_cloglog, glm_link_identity, glm_link_inverse
   public :: glm_link_log, glm_link_logit, glm_link_probit, glm_link_sqrt
   public :: iqr, mad, median, quantile, sd, var, weighted_mean
   public :: histogram, histogram_t, predict_ecdf
   public :: interpolation_constant, interpolation_linear, interpolation_rule_constant
   public :: interpolation_rule_missing, interpolation_t, isoreg, isoreg_t, predict_approxfun
   public :: integrate, integrate_result_t, optimize, optimize_result_t
   public :: numerical_invalid_input, numerical_iteration_limit, scalar_function_interface
   public :: root_extend_both, root_extend_lower, root_extend_none, root_extend_upper
   public :: uniroot, uniroot_result_t
   public :: predict_splinefun, spline, spline_method_fmm, spline_method_natural, splinefun
   public :: hclust, hclust_condensed, hclust_result_t, kmeans, kmeans_matrix, &
             kmeans_result_t, kmeans_vector, kruskal_test, kruskal_test_result_t, ks_test, ks_test_normal, &
             ks_test_result_t, ks_test_two_sample, &
             glm_predict_response, lm_aic, lm_coef, lm_confint, lm_cooks_distance, &
             lm_fit, lm_fit_general, lm_fit_t, lm_predict_general, &
             lm_r_squared_general, mahalanobis, prcomp, prcomp_fit_t, predict_lm, prop_test, &
             prop_test_result_t, scale, scale_columns, scale_matrix, scale_result_t, scale_vector, &
             smooth_xy_t
   public :: holt_winters, holt_winters_auto, holt_winters_fit_t, holt_winters_forecast_t
   public :: holt_winters_optimize, holt_winters_reduced, predict_holt_winters
   public :: predict_holt_winters_interval
   public :: multivariate_spectrum_result_t, spec_ar, spec_pgram, spectrum, spectrum_result_t
   public :: manova_lm, manova_result_t, mlm_fit, mlm_fit_t
   public :: one_way_anova, one_way_anova_t
   public :: nls_basis_interface, nls_confint_profile, nls_fit, nls_fit_plinear, nls_fit_t
   public :: nls_jacobian_interface, nls_max_iterations, nls_model_interface
   public :: seasonal_decomposition_t
   public :: ss_asymp, ss_asymp_jacobian, ss_asymp_model
   public :: ss_gompertz, ss_gompertz_jacobian, ss_gompertz_model
   public :: ss_logis, ss_logis_jacobian, ss_logis_model
   public :: ss_micmen, ss_micmen_jacobian, ss_micmen_model
   public :: ss_weibull, ss_weibull_jacobian, ss_weibull_model
   public :: stl, stl_result_t
   public :: t_test, t_test_one, t_test_p_value, &
             t_test_p_value_one, t_test_p_value_two, t_test_result_t, t_test_two, &
             wilcox_test, wilcox_test_result_t, wilcox_test_two_sample
   public :: poisson_test, var_test, variance_test_result_t
   public :: ksmooth, loess_fit, loess_fit_t, lowess, predict_loess, runmed
   public :: predict_smooth_spline, smooth_spline, smooth_spline_fit_t

end module r_stats
