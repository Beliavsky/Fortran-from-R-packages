! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran computational translation of sharpeRratio 1.4.3.
module sharpe_rratio
   use ghyp_kinds, only : dp
   use sharpe_rratio_calibration, only : calibration_a, calibration_a_medium, calibration_f
   use sharpe_rratio_records, only : r0_result, num_records_up, num_records_down, &
      compute_r0bar, computeR0bar
   use sharpe_rratio_statistics, only : test_n, a_full, f_full, &
      correction_b, theta_snr, sample_variance, quantile_type7
   use sharpe_rratio_estimator, only : snr_result, gaussian_nu, estimate_snr, &
      estimateSNR, estimate_tail_exponent
   implicit none
   public
end module sharpe_rratio
