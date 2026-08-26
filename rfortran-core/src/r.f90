! SPDX-License-Identifier: MIT
! Compatibility facade. Prefer importing directly from the defining module.
module r_mod
   use r_descriptive, only : r_correlation, r_count_nonmissing, r_covariance, r_mean, r_sd, r_variance
   use r_descriptive, only : r_weighted_correlation, r_weighted_covariance, r_weighted_mean
   use r_descriptive, only : r_weighted_sd, r_weighted_variance
   use r_distributions, only : r_dnorm, r_pnorm, r_qnorm
   use r_kinds, only : dp, i64, r_pi
   use r_missing, only : r_is_finite, r_is_na
   use r_ordering, only : r_order, r_sort_values_in_place
   use r_quantiles, only : r_qrule_hf1, r_qrule_hf2, r_qrule_hf3, r_qrule_hf4
   use r_quantiles, only : r_qrule_hf5, r_qrule_hf6, r_qrule_hf7, r_qrule_hf8, r_qrule_hf9
   use r_quantiles, only : r_qrule_math, r_qrule_school, r_qrule_shahvaish
   use r_quantiles, only : r_median, r_quantile_type7, r_weighted_quantile_ecdf
   use r_quantiles, only : r_weighted_quantile_frequency_type7, r_weighted_quantile_isotone
   use r_quantiles, only : r_weighted_quantile_linear_cdf, r_weighted_quantile_survey
   use r_rolling, only : r_roll_correlation_right, r_roll_correlation_valid, r_roll_covariance_right
   use r_rolling, only : r_roll_covariance_valid, r_roll_mean_right, r_roll_mean_valid
   use r_rolling, only : r_roll_max_right, r_roll_max_valid, r_roll_min_right, r_roll_min_valid
   use r_rolling, only : r_roll_sd_right, r_roll_sd_valid, r_roll_sum_right, r_roll_sum_valid
   use r_rolling, only : r_roll_variance_right, r_roll_variance_valid
   use r_robust, only : r_mad
   use r_special, only : r_digamma, r_log_beta, r_log_choose, r_log_factorial, r_trigamma
   use r_special, only : r_regularized_beta, r_regularized_gamma_p, r_regularized_gamma_q
   use r_status, only : r_invalid_input, r_no_data, r_ok, r_singular
   use r_sorting, only : r_average_ranks, r_sort
   use r_stability, only : r_log_mean_exp, r_log_sum_exp
   use r_time_series, only : r_autocorrelation, r_autocovariance, r_cross_correlation, r_cross_covariance
   use r_transforms, only : r_expm1, r_log1mexp, r_log1p, r_log1pexp, r_logistic, r_logit
   use r_vectors, only : r_difference
   implicit none
   private

   public :: dp, i64, r_pi
   public :: r_ok, r_invalid_input, r_no_data, r_singular
   public :: r_is_na, r_is_finite
   public :: r_order, r_sort_values_in_place
   public :: r_correlation, r_count_nonmissing, r_covariance, r_mean, r_variance, r_sd
   public :: r_weighted_correlation, r_weighted_covariance, r_weighted_mean
   public :: r_weighted_sd, r_weighted_variance
   public :: r_dnorm, r_pnorm, r_qnorm
   public :: r_sort, r_median, r_mad, r_quantile_type7, r_average_ranks
   public :: r_weighted_quantile_ecdf, r_weighted_quantile_linear_cdf
   public :: r_weighted_quantile_frequency_type7, r_weighted_quantile_isotone
   public :: r_weighted_quantile_survey
   public :: r_qrule_math, r_qrule_school, r_qrule_shahvaish
   public :: r_qrule_hf1, r_qrule_hf2, r_qrule_hf3, r_qrule_hf4, r_qrule_hf5
   public :: r_qrule_hf6, r_qrule_hf7, r_qrule_hf8, r_qrule_hf9
   public :: r_log_sum_exp, r_log_mean_exp
   public :: r_digamma, r_trigamma, r_log_beta, r_log_factorial, r_log_choose
   public :: r_regularized_beta, r_regularized_gamma_p, r_regularized_gamma_q
   public :: r_roll_mean_valid, r_roll_mean_right
   public :: r_roll_sum_valid, r_roll_sum_right
   public :: r_roll_min_valid, r_roll_min_right, r_roll_max_valid, r_roll_max_right
   public :: r_roll_covariance_valid, r_roll_covariance_right
   public :: r_roll_correlation_valid, r_roll_correlation_right
   public :: r_roll_variance_valid, r_roll_variance_right, r_roll_sd_valid, r_roll_sd_right
   public :: r_log1p, r_expm1, r_log1mexp, r_log1pexp
   public :: r_logistic, r_logit
   public :: r_autocovariance, r_autocorrelation, r_cross_covariance, r_cross_correlation
   public :: r_difference

end module r_mod
