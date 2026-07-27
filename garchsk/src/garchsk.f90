! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from GARCHSK 0.1.0, Copyright (C) 2021 Kei Nakagawa.
module garchsk
   use garchsk_kinds, only : dp
   use garchsk_types, only : moment_path, forecast_result, estimate_result
   use garchsk_stats, only : skewness => sample_skewness, kurtosis => sample_kurtosis
   use garchsk_models, only : garchsk_construct, gjrsk_construct, &
      garchsk_lik => garchsk_negative_log_likelihood, gjrsk_lik => gjrsk_negative_log_likelihood, &
      garchsk_ineqfun => garchsk_constraints, gjrsk_ineqfun => gjrsk_constraints, &
      garchsk_fcst => garchsk_forecast, gjrsk_fcst => gjrsk_forecast, &
      garchsk_parameters_valid, gjrsk_parameters_valid
   use garchsk_estimation, only : garchsk_est => garchsk_estimate, gjrsk_est => gjrsk_estimate, &
      garchsk_initial_parameters, gjrsk_initial_parameters
   implicit none
   public
end module garchsk
