! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yieldcurves
  use yc_kinds, only : dp
  use yc_types, only : status_t, curve_t, series_t, duration_result_t
  use yc_types, only : bond_duration_result_t, zspread_result_t, carry_result_t
  use yc_types, only : slope_result_t, factor_result_t, pca_result_t
  use yc_models, only : yc_curve, yc_nelson_siegel, yc_svensson, yc_cubic_spline, yc_fit
  use yc_models, only : ns_rate_scalar, sv_rate_scalar, ns_forward_scalar, sv_forward_scalar
  use yc_models, only : ns_loadings_matrix, sv_loadings_matrix
  use yc_curve_ops, only : yc_predict, yc_interpolate, yc_discount, yc_forward
  use yc_analysis, only : yc_par_to_zero, yc_zero_to_par, yc_duration, yc_bond_duration
  use yc_analysis, only : yc_zspread, yc_key_rate_duration, yc_carry, yc_slope
  use yc_analysis, only : yc_level_slope_curvature
  use yc_pca_mod, only : yc_pca
  implicit none
  public
end module yieldcurves
