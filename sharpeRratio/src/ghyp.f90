! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran numerical core derived from ghyp 1.6.5.
module ghyp
   use ghyp_kinds, only : dp, i8
   use ghyp_special, only : normal_pdf, normal_cdf, normal_quantile, student_cdf, &
      bessel_k, log_bessel_k
   use ghyp_gig, only : gig_valid, dgig, log_dgig, pgig, qgig, rgig, &
      gig_raw_moment, gig_mean, gig_variance, gig_mean_log, gig_mean_inverse, esgig
   use ghyp_model, only : ghyp_model_type, moments_result, make_ghyp, ghyp_uv, ghyp_mv, &
      hyp_uv, nig_uv, student_t_uv, vg_uv, gaussian_uv, gaussian_mv, ghyp_ad, &
      alpha_bar_to_chi_psi, ghyp_moments, transform_ghyp, ghyp_family_name, &
      model_ghyp, model_hyp, model_nig, model_student, model_vg, model_gaussian
   use ghyp_distribution, only : probability_result, log_dghyp, dghyp, pghyp, &
      pghyp_rectangle, qghyp, rghyp, rghyp_one
   use ghyp_risk, only : attribution_result, ghyp_moment, ghyp_skewness, &
      ghyp_kurtosis, esghyp, ghyp_omega, esghyp_attribution
   use ghyp_fitting, only : fit_result, likelihood_ratio_result, model_selection_result, &
      fit_ghyp_uv, fit_ghyp_mv, fit_gaussian_uv, fit_gaussian_mv, &
      likelihood_ratio_test, step_aic_ghyp
   use ghyp_portfolio, only : portfolio_result, portfolio_optimize
   use ghyp_utilities, only : alpha_delta_result, qq_result, subset_ghyp, &
      standardize_ghyp, ghyp_alpha_delta, qqghyp_data
   implicit none
   public
end module ghyp
