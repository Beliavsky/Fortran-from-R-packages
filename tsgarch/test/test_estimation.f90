program test_estimation
  use ghyp_kinds, only : dp, i8
  use tsgarch
  use test_support
  implicit none
  real(dp)::dummy(30),ll0
  type(garch_spec)::spec
  type(garch_parameters)::par,start
  type(garch_simulation)::sim
  type(garch_fit)::fit
  type(fit_options)::opt
  integer::i
  dummy=[(0.01_dp*sin(real(i,dp)),i=1,30)]
  spec=standard_spec('garch','norm')
  par=standard_parameters(dummy,spec)
  sim=simulate_garch(spec,par,260,burn=60,seed=2026_i8)
  call assert_true(sim%status==tsg_success,'simulation for estimation failed')
  start=initialize_parameters(sim%series(:,1),spec)
  ll0=garch_loglikelihood(sim%series(:,1),spec,start)
  opt%max_iterations=900
  opt%tolerance=1.0e-6_dp
  opt%compute_inference=.false.
  fit=estimate_garch(sim%series(:,1),spec,start=start,options=opt)
  call assert_true(fit%filtered%status==tsg_success,'estimation produced invalid filter')
  call assert_true(fit%log_likelihood>ll0,'estimation did not improve likelihood')
  call assert_true(fit%parameters%alpha(1)>=0.0_dp.and.fit%parameters%beta(1)>=0.0_dp,'invalid fitted coefficients')
  call assert_true(fit%filtered%persistence<1.0_dp,'fitted persistence not stationary')
  write(*,'(a,3f11.5)')'test_estimation: PASS parameters ',fit%parameters%omega,fit%parameters%alpha(1),fit%parameters%beta(1)
end program test_estimation
