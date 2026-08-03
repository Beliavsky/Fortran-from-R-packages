program example_rolling_marginals
   use portvine
   use rugarch, only : garch_spec,simulate_garch,seed_rng,dist_norm
   implicit none
   real(dp)::r(1,90),sigma0(90),eps0(90),mu,sigma,u
   type(garch_spec)::spec
   type(marginal_settings_type)::settings
   type(asset_marginal_result),allocatable::fit(:)
   integer::status
   call seed_rng(31415)
   spec=make_portvine_spec(0,0,1,1,cond_dist=dist_norm)
   spec%omega=2.0e-5_dp;spec%alpha=[0.08_dp];spec%beta=[0.88_dp]
   call simulate_garch(spec,90,r(1,:),sigma0,eps0,burn_in=300)
   settings=make_marginal_settings(60,15,1,spec);settings%max_iterations=250
   call fit_rolling_marginals(r,settings,30,fit,status)
   call get_marginal_point(fit(1),1,61,mu,sigma,u,status=status)
   print '(a,f10.6)','forecast mean:  ',mu
   print '(a,f10.6)','forecast sigma: ',sigma
   print '(a,f10.6)','PIT residual:   ',u
end program example_rolling_marginals
