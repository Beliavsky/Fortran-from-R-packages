program test_simulation_forecast
  use ghyp_kinds, only : dp, i8
  use tsgarch
  use test_support
  implicit none
  real(dp)::dummy(40),probs(3)
  type(garch_spec)::spec
  type(garch_parameters)::par
  type(garch_simulation)::sim,csim
  type(garch_fit)::fit
  type(garch_forecast)::fc
  integer::i
  dummy=[(0.01_dp*sin(real(i,dp)),i=1,40)]
  spec=standard_spec('garch','std')
  par=standard_parameters(dummy,spec)
  par%dist%shape=8.0_dp
  sim=simulate_garch(spec,par,120,paths=2,burn=30,seed=12345_i8)
  call assert_true(sim%status==tsg_success,'unconditional simulation failed')
  call assert_true(all(shape(sim%series)==[120,2]),'simulation dimensions')
  csim=simulate_conditional(sim%series(:,1),spec,par,5,paths=30,seed=77_i8)
  call assert_true(csim%status==tsg_success,'conditional simulation failed')
  fit%spec=spec
  fit%parameters=par
  fit%filtered=filter_garch(sim%series(:,1),spec,par)
  fit%status=tsg_success
  probs=[0.05_dp,0.5_dp,0.95_dp]
  fc=forecast_garch(sim%series(:,1),fit,5,probs,paths=100,seed=99_i8)
  call assert_true(fc%status==tsg_success,'forecast failed')
  call assert_true(all(fc%quantiles(:,1)<=fc%quantiles(:,2)).and. &
    all(fc%quantiles(:,2)<=fc%quantiles(:,3)),'forecast quantile order')
  call assert_true(all(fc%sigma>0.0_dp),'forecast standard deviations')
  write(*,'(a)')'test_simulation_forecast: PASS'
end program test_simulation_forecast
