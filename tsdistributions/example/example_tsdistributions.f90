program example_tsdistributions
  use tsdistributions
  implicit none
  type(rng_state) :: rng
  type(distribution_parameters) :: truth
  type(parameter_specification) :: spec
  type(distribution_fit) :: fit
  type(spd_fit) :: tail_fit
  real(dp), allocatable :: y(:)

  call seed_rng(rng,20260804_i8)
  truth=distribution_parameters(mu=0.02_dp,sigma=1.1_dp,skew=1.25_dp,shape=7.0_dp)
  y=rdist('sstd',1000,rng,truth)
  spec=distribution_modelspec(y,'sstd')
  fit=estimate_distribution(y,spec,max_iterations=1200)

  write(*,'(a,i0)') 'fit status: ',fit%status
  write(*,'(a,4f12.5)') 'mu sigma skew shape: ',fit%parameters%mu,fit%parameters%sigma, &
    fit%parameters%skew,fit%parameters%shape
  write(*,'(a,3f12.5)') 'logLik AIC BIC: ',fit%log_likelihood,fit%aic,fit%bic
  write(*,'(a,f12.5)') 'one-percent quantile: ',qdist('sstd',0.01_dp,fit%parameters)

  tail_fit=estimate_spd(spd_modelspec(y,0.1_dp,0.9_dp,'normal'))
  write(*,'(a,2f12.5)') 'SPD tail shapes: ',tail_fit%lower_shape,tail_fit%upper_shape
  write(*,'(a,f12.5)') 'SPD one-percent quantile: ',qspd(0.01_dp,tail_fit)
end program example_tsdistributions
