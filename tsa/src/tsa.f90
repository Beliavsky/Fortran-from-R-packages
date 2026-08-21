! SPDX-License-Identifier: GPL-2.0-or-later
module tsa
  use tsa_kinds, only : dp
  use tsa_types, only : tsa_test_result, spectrum_result, ar_fit_result, &
    arimax_result, tar_result, runs_result, outlier_result, bootstrap_result, &
    transfer_spec, spectral_estimate, tar_multi_result
  use tsa_utils, only : lag_vector
  use tsa_statistics, only : skewness, kurtosis, autocorrelation, &
    autocovariance, partial_autocorrelation, cross_correlation, periodogram, &
    harmonic_matrix, season_index, runs
  use tsa_arma, only : ar_ols_fit, arma_spectrum, eacf, armasubsets_fit, &
    boxcox_ar, prewhiten_filter
  use tsa_spectral, only : spec_pgram, spec_ar, modified_daniell_weights, &
    spectral_kernel_df, spectral_kernel_bandwidth, tskernel_weights, &
    daniell_kernel, modified_daniell_kernel, fejer_kernel, dirichlet_kernel
  use tsa_tests, only : lb_test, mcleod_li_test, keenan_test, tsay_test, &
    detect_io, detect_ao
  use tsa_simulation, only : qar_sim, garch_sim, tar_sim
  use tsa_tar, only : tar_fit, tar_fit_multi, tar_skeleton, tar_predict, tlrt_test, &
    tlrt_p_value
  use tsa_arimax, only : arima_fit, arimax_fit, arima_sim, arima_bootstrap, arima_bootstrap_sample, &
    transfer_filter, io_regressor
  use tsa_garch_diagnostics, only : gbox_test
  implicit none
  public
end module tsa
