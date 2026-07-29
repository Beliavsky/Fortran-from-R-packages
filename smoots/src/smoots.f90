! SPDX-License-Identifier: GPL-3.0-only
module smoots
   use smoots_kinds
   use smoots_status
   use smoots_types
   use smoots_stats, only : seed_rng
   use smoots_smoothing, only : gsmooth, knsmooth, local_polynomial_smooth, kernel_smooth, &
      lag_window_variance, rescale_derivative
   use smoots_arma, only : fit_arma, arma_residuals, arma_point_forecast, simulate_arma, &
      information_criterion_matrix, optimal_order, ma_infinity, estimate_cf0_ar, &
      estimate_cf0_ma, estimate_cf0_arma
   use smoots_estimation, only : tsmooth, msmooth, dsmooth, fixed_gsmooth, fixed_knsmooth, conf_bounds
   use smoots_forecast, only : trend_forecast, normal_forecast, bootstrap_forecast, model_forecast, rolling_backtest
   implicit none
   public
end module smoots
