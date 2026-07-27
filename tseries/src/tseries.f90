! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries
   use tseries_kinds, only : dp
   use tseries_types, only : test_result, bds_result, arma_result, garch_result, drawdown_result, portfolio_result
   use tseries_models, only : arma_fit, arma_residuals, garch_fit, garch_variance
   use tseries_tests, only : runs_test, jarque_bera_test, adf_test, pp_test, kpss_test, po_test, bds_test, &
      terasvirta_test, white_test
   use tseries_resampling, only : quadratic_map, stationary_bootstrap, block_bootstrap, permutation_surrogate, &
      fft_surrogate, amplitude_surrogate
   use tseries_finance, only : maximum_drawdown, sharpe_ratio, sterling_ratio, portfolio_optimize
   use tseries_random, only : seed_random
   implicit none
   private

   public :: dp
   public :: test_result, bds_result, arma_result, garch_result, drawdown_result, portfolio_result
   public :: arma_fit, arma_residuals, garch_fit, garch_variance
   public :: runs_test, jarque_bera_test, adf_test, pp_test, kpss_test, po_test, bds_test
   public :: terasvirta_test, white_test
   public :: quadratic_map, stationary_bootstrap, block_bootstrap, permutation_surrogate
   public :: fft_surrogate, amplitude_surrogate
   public :: maximum_drawdown, sharpe_ratio, sterling_ratio, portfolio_optimize
   public :: seed_random

end module tseries
