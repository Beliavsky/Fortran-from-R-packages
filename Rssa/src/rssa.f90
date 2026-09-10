! SPDX-License-Identifier: GPL-2.0-or-later
! Public umbrella module for the modern Fortran translation of Rssa 1.1.
module rssa
   use rssa_kinds, only : dp
   use rssa_types, only : cssa_result, gap_summary, grouping_result, hbhmat_type, hmat_type, mssa_result
   use rssa_types, only : period_estimate, rssa_invalid_input, rssa_not_supported, rssa_numerical_failure, rssa_success
   use rssa_types, only : ssa2d_result, ssa_result, tmat_type
   use rssa_matrices
   use rssa_decomposition
   use rssa_projection
   use rssa_reconstruction
   use rssa_metrics
   use rssa_oblique
   use rssa_iterative_oblique
   use rssa_forecast
   use rssa_parestimate
   use rssa_gapfill
   use rssa_cadzow
   use rssa_autogroup
   use rssa_hmatr
   implicit none
   private
   public :: dp
   public :: rssa_success, rssa_invalid_input, rssa_numerical_failure, rssa_not_supported
   public :: ssa_result, mssa_result, ssa2d_result, cssa_result
   public :: hmat_type, hbhmat_type, tmat_type, period_estimate, gap_summary, grouping_result
   public :: hankel, hankel_matrix, hankelize_matrix, complex_hankel_matrix, complex_hankelize_matrix
   public :: new_hmat, hmatmul, hcols, hrows, is_hmat
   public :: new_hbhmat, hbhmatmul, hbhcols, hbhrows, is_hbhmat
   public :: new_tmat, tmatmul, tcols, trows, is_tmat, lag_covariance
   public :: trajectory_2d, hankelize_2d, hankel_weights, hankel_weights_2d
   public :: ssa, ssa_1d, ssa_mssa, ssa_2d, ssa_toeplitz, ssa_complex
   public :: decompose_ssa, decompose_mssa, decompose_2d, decompose_toeplitz, decompose_complex
   public :: decompose_pssa, calc_v_pssa, nspecial_pssa
   public :: calc_v_ssa, calc_v_complex, clone_ssa, nu, nv, nlambda, nsigma, nspecial, contributions
   public :: reconstruct_ssa, reconstruct_mssa, reconstruct_2d, reconstruct_complex
   public :: elementary_series_ssa, elementary_series_mssa, elementary_field_2d, elementary_series_complex
   public :: residuals_ssa, residuals_mssa, residuals_2d, residuals_complex
   public :: wnorm, wnorm_default, wnorm_complex, wnorm_ssa, wnorm_mssa, wnorm_2d
   public :: wcor_default, wcor_ssa, frobenius_cor
   public :: decompose_ossa, decompose_wossa, fossa, fossa_ssa, owcor_ssa, wcor_ossa
   public :: iossa, iossa_ssa, eossa, eossa_ssa
   public :: lrr, lrr_default, lrr_ssa, lrr_mssa, lrr_complex
   public :: roots_lrr, roots_lrr_complex, apply_lrr, apply_lrr_complex
   public :: rforecast_ssa, rforecast_mssa, rforecast_complex, rforecast_pssa
   public :: vforecast_ssa, vforecast_mssa, vforecast_complex, vforecast_pssa, bforecast_ssa
   public :: roots_to_parameters, parestimate_ssa, parestimate_mssa, parestimate_complex, parestimate_2d
   public :: summarize_gaps, summarize_gaps_complex
   public :: gapfill_ssa, gapfill_mssa_channel, gapfill_complex
   public :: igapfill, igapfill_ssa, igapfill_2d, igapfill_mssa, igapfill_complex
   public :: cadzow, cadzow_ssa
   public :: grouping_auto, grouping_auto_wcor_ssa, grouping_auto_pgram_ssa
   public :: hmatr
end module rssa
