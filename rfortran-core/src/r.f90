! SPDX-License-Identifier: MIT
! Compatibility facade. Prefer importing directly from the defining module.
module r_mod
   use r_descriptive, only : r_count_nonmissing, r_mean, r_sd, r_variance
   use r_distributions, only : r_dnorm, r_pnorm
   use r_kinds, only : dp, i64, r_pi
   use r_missing, only : r_is_finite, r_is_na
   use r_rolling, only : r_roll_mean_right, r_roll_mean_valid
   use r_special, only : r_digamma, r_log_beta, r_trigamma
   use r_status, only : r_invalid_input, r_no_data, r_ok, r_singular
   use r_sorting, only : r_average_ranks, r_quantile_type7, r_sort
   use r_stability, only : r_log_mean_exp, r_log_sum_exp
   use r_time_series, only : r_autocorrelation, r_autocovariance, r_cross_correlation, r_cross_covariance
   use r_vectors, only : r_difference
   implicit none
   private

   public :: dp, i64, r_pi
   public :: r_ok, r_invalid_input, r_no_data, r_singular
   public :: r_is_na, r_is_finite
   public :: r_count_nonmissing, r_mean, r_variance, r_sd
   public :: r_dnorm, r_pnorm
   public :: r_sort, r_quantile_type7, r_average_ranks
   public :: r_log_sum_exp, r_log_mean_exp
   public :: r_digamma, r_trigamma, r_log_beta
   public :: r_roll_mean_valid, r_roll_mean_right
   public :: r_autocovariance, r_autocorrelation, r_cross_covariance, r_cross_correlation
   public :: r_difference

end module r_mod
