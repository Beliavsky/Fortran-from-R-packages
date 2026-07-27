! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

program complete_workflows
   use rugarch
   implicit none
   integer,parameter::n=250
   type(garch_spec)::spec
   type(extended_garch_fit_result)::fit
   type(bootstrap_forecast_result)::boot
   real(dp)::y(n),sigma(n),eps(n),mx(n,1),vx(n,1)
   integer::i

   call seed_rng(20260722)
   spec=make_garch_spec(1,1,model_sgarch,dist_norm)
   spec%omega=0.03_dp;spec%alpha(1)=0.08_dp;spec%beta(1)=0.88_dp
   call simulate_garch(spec,n,y,sigma,eps,burn_in=300)
   do i=1,n
      mx(i,1)=sin(0.02_dp*real(i,dp))
      vx(i,1)=0.01_dp+0.005_dp*cos(0.03_dp*real(i,dp))**2
   end do
   y=y+0.2_dp+0.3_dp*mx(:,1)

   fit=fit_garch_extended(y,model_sgarch,1,1,dist_norm, &
      mean_regressors=mx,variance_regressors=vx,variance_targeting=.true., &
      max_iterations=500)
   boot=garch_bootstrap_forecast(fit%fit,5,100,sampling_kernel,bootstrap_partial)

   print '(a,f10.5)','mean-regressor coefficient: ',fit%mean_beta(1)
   print '(a,f10.5)','target variance: ',fit%target_variance
   print '(a,i0)','covariance status: ',fit%covariance_status
   print '(a,5f10.5)','bootstrap mean: ',boot%mean
end program complete_workflows
