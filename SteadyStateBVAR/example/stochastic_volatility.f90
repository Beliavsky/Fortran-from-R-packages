! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from SteadyStateBVAR 0.2.0, Copyright (c) 2026 Mark Becker.
program stochastic_volatility
   use steadystatebvar, only : dp, bvar_model, forecast_result, bvar_summary_result
   use steadystatebvar, only : bvar_create, setup_model, set_priors, fit_model
   use steadystatebvar, only : forecast_model, bvar_summary, ssbvar_success
   implicit none

   integer, parameter :: n = 50, k = 2
   real(dp) :: data(n, k), theta_a(1), omega_a(1, 1)
   real(dp) :: theta_log_lambda_1(k), omega_log_lambda_1(k, k)
   real(dp) :: alpha_phi(k), beta_phi(k)
   integer :: t, info
   type(bvar_model) :: model
   type(forecast_result) :: forecast
   type(bvar_summary_result) :: summary

   do t = 1, n
      data(t, 1) = 0.4_dp * sin(0.15_dp * real(t, dp)) + 0.01_dp * real(t, dp)
      data(t, 2) = 0.3_dp * data(t, 1) + 0.5_dp * cos(0.11_dp * real(t, dp))
   end do

   theta_a = 0.0_dp
   omega_a = 0.0_dp
   omega_a(1, 1) = 10.0_dp
   theta_log_lambda_1 = -1.0_dp
   omega_log_lambda_1 = 0.0_dp
   omega_log_lambda_1(1, 1) = 1.0_dp
   omega_log_lambda_1(2, 2) = 1.0_dp
   alpha_phi = 5.0_dp
   beta_phi = 0.4_dp

   call bvar_create(data, model, info)
   if (info /= ssbvar_success) error stop "bvar_create failed"
   call setup_model(model, 1, "constant", info=info)
   if (info /= ssbvar_success) error stop "setup_model failed"
   call set_priors(model, sv=.true., sv_type="RW", theta_a=theta_a, omega_a=omega_a, &
      theta_log_lambda_1=theta_log_lambda_1, omega_log_lambda_1=omega_log_lambda_1, &
      alpha_phi=alpha_phi, beta_phi=beta_phi, info=info)
   if (info /= ssbvar_success) error stop "set_priors failed"
   call fit_model(model, h=3, n_draws=30, burnin=20, seed=2468, sv_state_scale=0.15_dp, info=info)
   if (info /= ssbvar_success) error stop "fit_model failed"
   call forecast_model(model, pi=0.90_dp, result=forecast, info=info)
   if (info /= ssbvar_success) error stop "forecast_model failed"
   call bvar_summary(model, result=summary, info=info, t=n)
   if (info /= ssbvar_success) error stop "bvar_summary failed"

   print '(a,2(1x,f12.6))', "RW innovation variances:", summary%phi
   print '(a)', "Posterior predictive mean forecast:"
   do t = 1, size(forecast%forecast, 1)
      print '(i3,2(1x,f12.6))', t, forecast%forecast(t, :)
   end do
end program stochastic_volatility
