! SPDX-License-Identifier: GPL-2.0-or-later
module skellam
   use skellam_kinds, only : dp, i8
   use skellam_special, only : seed_random_number
   use skellam_distribution, only : dskellam, pskellam, qskellam, rskellam, &
      dskellam_sp, pskellam_sp, skellam_log_pmf, skellam_mean, skellam_variance, &
      skellam_skewness, skellam_excess_kurtosis
   use skellam_estimation, only : skellam_mle_result, skellam_regression_result, &
      fit_skellam_mle, fit_skellam_regression, skellam_log_likelihood
   implicit none
   public
end module skellam
