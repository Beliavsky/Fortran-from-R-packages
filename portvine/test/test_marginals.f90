program test_marginals
   use portvine, only : dp, marginal_settings_type, asset_marginal_result, &
      make_portvine_spec, fit_rolling_marginals, get_marginal_point, portvine_success
   use rugarch, only : garch_spec, simulate_garch, seed_rng, dist_norm
   implicit none
   real(dp) :: returns(1,90),sigma0(90),eps0(90),mu,sigma,u
   type(garch_spec) :: simspec
   type(marginal_settings_type) :: settings
   type(asset_marginal_result),allocatable :: fit(:)
   integer :: status
   call seed_rng(99173)
   simspec=make_portvine_spec(0,0,1,1,cond_dist=dist_norm)
   simspec%mean=0.001_dp;simspec%omega=2.0e-5_dp
   simspec%alpha=[0.08_dp];simspec%beta=[0.88_dp]
   call simulate_garch(simspec,90,returns(1,:),sigma0,eps0,burn_in=300)
   settings%train_size=60;settings%refit_size=15;settings%max_iterations=250
   allocate(settings%spec(1));settings%spec(1)=simspec
   call fit_rolling_marginals(returns,settings,30,fit,status)
   if(status/=portvine_success .or. size(fit)/=1 .or. size(fit(1)%window)/=2)error stop 1
   call get_marginal_point(fit(1),1,61,mu,sigma,u,status=status)
   if(status/=0 .or. sigma<=0.0_dp .or. u<=0.0_dp .or. u>=1.0_dp)error stop 2
   print '(a)', 'test_marginals: PASS'
end program test_marginals
