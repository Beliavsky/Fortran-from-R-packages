! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam
   use vgam_kinds
   use vgam_special
   use vgam_links
   use vgam_random
   use vgam_distributions
   use vgam_actuarial
   use vgam_extremes
   use vgam_linalg
   use vgam_optim
   use vgam_vglm
   use vgam_categorical
   use vgam_extended_models
   use vgam_constraints
   use vgam_information
   use vgam_count_models
   use vgam_inflated
   use vgam_zero_altered
   use vgam_zero_altered_models
   use vgam_zoa_beta_models
   use vgam_gaitd
   use vgam_gaitd_regression
   use vgam_gaitd_mlm
   use vgam_gaitd_mix
   use vgam_gaitd_mix_regression
   use vgam_gaitd_nb_dispersion
   use vgam_copulas
   use vgam_student_t
   use vgam_bivariate_extra
   use vgam_multivariate_extra
   use vgam_dirichlet
   use vgam_normal_special
   use vgam_positive_count
   use vgam_censored
   use vgam_qreg
   use vgam_reduced_rank
   use vgam_drr
   use vgam_quadratic_rr
   use vgam_cqo
   use vgam_cao
   use vgam_timeseries
   use vgam_garma
   use vgam_rrar
   use vgam_smoothing
   use splines, only : b_spline_t, poly_spline_t, spline_design, &
      spline_basis_nonzero, linear_interp, fit_interpolating_spline, &
      fit_periodic_spline, to_polynomial_spline, inverse_monotone_spline, &
      type7_quantile, bs_basis, natural_spline_basis
   implicit none
   public
end module vgam
