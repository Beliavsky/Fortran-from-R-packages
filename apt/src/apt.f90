! SPDX-License-Identifier: GPL-2.0-or-later
module apt
  use apt_kinds, only : dp
  use apt_special, only : normal_cdf, student_t_cdf, f_cdf, chi_square_cdf
  use apt_regression, only : regression_result, hypothesis_result, &
    residual_diagnostics, fit_ols, linear_f_test, zero_coefficient_f_test, &
    ljung_box_test, durbin_watson_test
  use apt_cointegration, only : apt_tar, apt_mtar, ci_tar_fit_result, &
    ci_tar_lag_result, ci_tar_threshold_result, ci_tar_fit, ci_tar_lag, &
    ci_tar_threshold, ciTarFit, ciTarLag, ciTarThd
  use apt_ecm, only : apt_linear, ecm_fit_result, ecm_equation_diagnostics, &
    ecm_diagnostics_result, paired_hypothesis_result, &
    ecm_asymmetry_test_result, ecm_symmetric_fit, ecm_asymmetric_fit, &
    ecm_asymmetry_tests, ecm_diagnostics, ecmSymFit, ecmAsyFit, &
    ecmAsyTest, ecmDiag
  implicit none
  public
end module apt
