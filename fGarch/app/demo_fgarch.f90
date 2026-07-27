! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

program demo_fgarch
   use fgarch, only : dp, garch_spec, garch_fit_result, make_garch_spec, &
      simulate_garch, fit_garch11, forecast_volatility, value_at_risk, &
      expected_shortfall, seed_rng, dist_std, distribution_name
   implicit none

   integer, parameter :: n = 1200
   type(garch_spec) :: spec
   type(garch_fit_result) :: fit
   real(dp) :: y(n), sigma(n), residuals(n), forecast(5)
   real(dp) :: var01, es01

   call seed_rng(4711)
   spec = make_garch_spec(1,1,cond_dist=dist_std)
   spec%mean = 0.0002_dp
   spec%omega = 1.0e-6_dp
   spec%alpha = [0.08_dp]
   spec%beta = [0.90_dp]
   spec%shape = 7.0_dp

   call simulate_garch(spec,n,y,sigma,residuals)
   fit = fit_garch11(y,cond_dist=dist_std,fit_shape=.true.,max_iterations=1800)
   call forecast_volatility(fit%spec,fit%residuals,fit%sigma,size(forecast),forecast)
   var01 = value_at_risk(0.01_dp,fit%spec%mean,forecast(1),fit%spec%cond_dist, &
                         fit%spec%shape,fit%spec%skew)
   es01 = expected_shortfall(0.01_dp,fit%spec%mean,forecast(1),fit%spec%cond_dist, &
                             fit%spec%shape,fit%spec%skew)

   print '(a,i0)', 'observations: ', n
   print '(a,a)', 'conditional distribution: ', trim(distribution_name(fit%spec%cond_dist))
   print '(a,f12.8)', 'estimated mean:  ', fit%spec%mean
   print '(a,es14.6)', 'estimated omega: ', fit%spec%omega
   print '(a,f10.6)', 'estimated alpha: ', fit%spec%alpha(1)
   print '(a,f10.6)', 'estimated beta:  ', fit%spec%beta(1)
   print '(a,f10.4)', 'estimated shape: ', fit%spec%shape
   print '(a,f14.3)', 'log likelihood:  ', fit%log_likelihood
   print '(a,f12.8)', 'one-step sigma:  ', forecast(1)
   print '(a,f12.8)', 'one-percent VaR: ', var01
   print '(a,f12.8)', 'one-percent ES:  ', es01
   print '(a,a)', 'optimizer status: ', trim(fit%message)
end program demo_fgarch
