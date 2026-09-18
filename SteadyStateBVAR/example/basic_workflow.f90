! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from SteadyStateBVAR 0.2.0, Copyright (c) 2026 Mark Becker.
program basic_workflow
   use steadystatebvar, only : dp, bvar_model, forecast_result
   use steadystatebvar, only : bvar_create, setup_model, set_priors, fit_model, forecast_model
   use steadystatebvar, only : ssbvar_success
   implicit none

   integer, parameter :: n = 60, k = 2
   real(dp) :: data(n, k), own_mean(k)
   integer :: t, info
   type(bvar_model) :: model
   type(forecast_result) :: forecast

   do t = 1, n
      data(t, 1) = sin(0.12_dp * real(t, dp)) + 0.01_dp * real(t, dp)
      data(t, 2) = 0.35_dp * data(t, 1) + cos(0.08_dp * real(t, dp))
   end do
   own_mean = 0.8_dp

   call bvar_create(data, model, info)
   if (info /= ssbvar_success) error stop "bvar_create failed"
   call setup_model(model, 1, "constant", info=info)
   if (info /= ssbvar_success) error stop "setup_model failed"
   call set_priors(model, first_own_lag_prior_mean=own_mean, jeffreys=.true., info=info)
   if (info /= ssbvar_success) error stop "set_priors failed"
   call fit_model(model, h=4, n_draws=100, burnin=50, seed=12345, info=info)
   if (info /= ssbvar_success) error stop "fit_model failed"
   call forecast_model(model, pi=0.90_dp, result=forecast, info=info)
   if (info /= ssbvar_success) error stop "forecast_model failed"

   print '(a)', "Posterior predictive mean forecast:"
   do t = 1, size(forecast%forecast, 1)
      print '(i3,2(1x,f12.6))', t, forecast%forecast(t, :)
   end do
end program basic_workflow
