! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

program test_models
   use fgarch, only : dp, garch_spec, garch_fit_result, make_garch_spec, &
      simulate_garch, garch_log_likelihood, fit_garch11, forecast_volatility, &
      seed_rng, dist_norm, dist_sged, model_aparch, garch_kappa, true_persistence
   implicit none

   integer, parameter :: n = 800
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   real(dp) :: y(n), sigma(n), residuals(n), llh, forecast(3), kappa, persistence

   call seed_rng(1953)
   spec = make_garch_spec(1,1,cond_dist=dist_norm)
   spec%mean = 0.0_dp
   spec%omega = 2.0e-6_dp
   spec%alpha = [0.10_dp]
   spec%beta = [0.85_dp]
   call simulate_garch(spec,n,y,sigma,residuals,burn_in=300)

   if (any(sigma <= 0.0_dp)) error stop 'nonpositive simulated volatility'
   llh = garch_log_likelihood(y,spec)
   if (llh < -1.0e7_dp .or. llh > 1.0e7_dp) error stop 'invalid likelihood'

   fit = fit_garch11(y,cond_dist=dist_norm,max_iterations=1400)
   if (fit%log_likelihood < llh-100.0_dp) error stop 'fit has implausibly poor likelihood'
   if (fit%spec%alpha(1) < 0.0_dp .or. fit%spec%beta(1) < 0.0_dp) error stop 'negative coefficient'
   if (fit%spec%alpha(1)+fit%spec%beta(1) >= 0.999_dp) error stop 'nonstationary fit'
   call forecast_volatility(fit%spec,fit%residuals,fit%sigma,size(forecast),forecast)
   if (any(forecast <= 0.0_dp)) error stop 'invalid volatility forecast'


   spec = make_garch_spec(1,1,model=model_aparch,cond_dist=dist_sged)
   spec%mean = 0.0_dp
   spec%omega = 1.0e-5_dp
   spec%alpha = [0.08_dp]
   spec%gamma = [0.15_dp]
   spec%beta = [0.86_dp]
   spec%delta = 1.4_dp
   spec%shape = 1.7_dp
   spec%skew = 1.2_dp
   call simulate_garch(spec,n,y,sigma,residuals,burn_in=300)
   llh = garch_log_likelihood(y,spec)
   if (llh < -1.0e7_dp .or. llh > 1.0e7_dp) error stop 'invalid APARCH likelihood'

   spec = make_garch_spec(1,1,cond_dist=dist_norm)
   spec%delta = 2.0_dp
   spec%gamma = [0.0_dp]
   kappa = garch_kappa(spec,0.0_dp)
   if (abs(kappa-1.0_dp) > 2.0e-5_dp) error stop 'normal GARCH kappa is not one'
   persistence = true_persistence(spec)
   if (abs(persistence-(spec%alpha(1)+spec%beta(1))) > 2.0e-5_dp) error stop 'persistence mismatch'

   print '(a)', 'model tests passed'
end program test_models
