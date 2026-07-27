! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
module markowitzr
   use markowitzr_kinds, only: dp
   use markowitzr_types, only: theta_result, markowitz_result
   use markowitzr_linalg, only: symmetric_vech, symmetric_ivech
   use markowitzr_linalg, only: duplication_matrix, kronecker_product
   use markowitzr_moments, only: theta_vcov, itheta_vcov
   use markowitzr_moments, only: covariance_empirical, covariance_normal
   use markowitzr_moments, only: covariance_hac, hac_covariance_of_mean
   use markowitzr_moments, only: moment_covariance_callback
   use markowitzr_portfolio, only: mp_vcov, markowitz_weights
   use markowitzr_portfolio, only: weights_upstream, weights_all_columns
   implicit none
   public
end module markowitzr
