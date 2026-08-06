program demo_tsgarch
  use ghyp_kinds, only : dp, i8
  use tsgarch
  implicit none
  real(dp)::dummy(30)
  type(garch_spec)::spec
  type(garch_parameters)::par
  type(garch_simulation)::sim
  type(garch_fit)::fit
  type(garch_forecast)::fc
  type(fit_options)::opt
  integer::i
  dummy=[(0.0_dp,i=1,30)]
  spec%model='gjrgarch'
  spec%distribution='std'
  spec%p=1
  spec%q=1
  spec%constant=.true.
  par=initialize_parameters(dummy,spec)
  par%omega=0.02_dp
  par%alpha=0.05_dp
  par%gamma=0.06_dp
  par%beta=0.86_dp
  par%dist%shape=8.0_dp
  sim=simulate_garch(spec,par,300,burn=80,seed=20260804_i8)
  opt%max_iterations=900
  opt%compute_inference=.false.
  fit=estimate_garch(sim%series(:,1),spec,options=opt)
  fc=forecast_garch(sim%series(:,1),fit,5,paths=300,seed=17_i8)
  write(*,'(a,f12.4)')'log likelihood: ',fit%log_likelihood
  write(*,'(a,4f10.5)')'omega alpha gamma beta: ',fit%parameters%omega, &
    fit%parameters%alpha(1),fit%parameters%gamma(1),fit%parameters%beta(1)
  write(*,'(a,5f10.5)')'forecast sigma: ',fc%sigma
end program demo_tsgarch
