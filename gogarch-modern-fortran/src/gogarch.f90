! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of gogarch, copyright (C) 2008-2026 Bernhard Pfaff.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module gogarch
   use gogarch_kinds, only : dp
   use gogarch_types, only : univariate_spec, garch11_fit, gogarch_fit
   use gogarch_orthogonal, only : rd2, uprod_r, umatch, unvech, vech, angle_dimension
   use gogarch_core, only : initialize_gogarch, cora, build_covariance_path
   use gogarch_core, only : conditional_variances, conditional_correlations
   use gogarch_distributions, only : distribution_is_valid, innovation_logpdf, innovation_pdf
   use gogarch_distributions, only : random_innovation, innovation_asym_power_moment
   use gogarch_distributions, only : symmetric_absolute_moment, fs_location_scale
   use gogarch_univariate, only : filter_aparch, filter_garchpq, filter_garch11
   use gogarch_univariate, only : fit_univariate, fit_garchpq, fit_garch11
   use gogarch_univariate, only : forecast_univariate, forecast_garch11
   use gogarch_univariate, only : simulate_aparch, simulate_garchpq, simulate_garch11
   use gogarch_univariate, only : validate_specification
   use gogarch_ica, only : ica_result, fastica
   use gogarch_estimators, only : fit_gogarch, fit_gogarch_ica, fit_gogarch_mm
   use gogarch_estimators, only : fit_gogarch_nls, fit_gogarch_ml, gonls_objective
   use gogarch_estimators, only : gogarch_from_angles, gogarch_negloglik
   use gogarch_model, only : forecast_gogarch, standardized_residuals
   use gogarch_model, only : simulate_fitted_gogarch, reconstruction_error, factor_coefficients
   use gogarch_model, only : factor_coefficients_full
   use gogarch_rng, only : seed_rng, random_uniform, random_normal, random_gamma
   use gogarch_rng, only : random_student_t, fill_random_normal
   implicit none
   public
end module gogarch
