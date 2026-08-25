! SPDX-License-Identifier: MIT
! Compatibility facade. Prefer importing directly from the defining module.
module r_mod
   use r_descriptive, only : r_count_nonmissing, r_mean, r_sd, r_variance
   use r_kinds, only : dp, i64, r_pi
   use r_missing, only : r_is_finite, r_is_na
   use r_status, only : r_invalid_input, r_no_data, r_ok, r_singular
   use r_time_series, only : r_autocorrelation, r_autocovariance, r_cross_correlation, r_cross_covariance
   use r_vectors, only : r_difference
   implicit none
   private

   public :: dp, i64, r_pi
   public :: r_ok, r_invalid_input, r_no_data, r_singular
   public :: r_is_na, r_is_finite
   public :: r_count_nonmissing, r_mean, r_variance, r_sd
   public :: r_autocovariance, r_autocorrelation, r_cross_covariance, r_cross_correlation
   public :: r_difference

end module r_mod
