! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail
   use fitheavytail_kinds, only: dp
   use fitheavytail_status
   use fitheavytail_types, only: heavy_tail_fit
   use fitheavytail_elliptical, only: fit_tyler, fit_cauchy
   use fitheavytail_mvt, only: fit_mvt, mvt_log_likelihood, nu_mle
   use fitheavytail_mvst, only: fit_mvst, mvst_log_likelihood, sample_skewness
   use fitheavytail_tail, only: default_nu_min, default_nu_max, &
      nu_opp_estimator, nu_pop_estimator, excess_kurtosis_unbiased, &
      nu_from_average_marginal_kurtosis, nu_from_cross_cumulants, &
      nu_from_all_cumulants, nu_hill_estimator, nu_pareto_tail_index
   use fitheavytail_special, only: digamma_dp, log_bessel_k, &
      bessel_k_ratio
   use fitheavytail_rng, only: random_mvt_identity
   implicit none
   public
end module fitheavytail
