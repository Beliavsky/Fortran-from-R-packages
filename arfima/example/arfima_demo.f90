program arfima_demo
  use arfima
  implicit none
  type(arfima_spec)::spec
  type(arfima_parameters)::truth
  type(arfima_fit_result)::fit
  type(arfima_forecast_result)::forecast
  type(arfima_error)::error
  real(dp),allocatable::series(:)
  integer::i

  call set_random_seed(20260804)
  spec%p=1
  spec%lmodel=long_memory_fd
  spec%estimate_mean=.true.
  allocate(truth%phi(1),truth%theta(0),truth%phiseas(0),truth%thetaseas(0))
  allocate(truth%beta(0),truth%delta(0),truth%omega(0))
  truth%phi=[0.35_dp]
  truth%dfrac=0.22_dp
  truth%mean=0.5_dp

  call arfima_simulate(spec,truth,120,0.8_dp,series,error,center_sample=.false.)
  if(error%code/=arfima_ok) error stop trim(error%message)
  call fit_arfima(spec,series,fit,max_iterations=1500,tolerance=2.0e-6_dp)
  if(fit%error%code/=arfima_ok .and. fit%error%code/=arfima_no_convergence) error stop trim(fit%error%message)
  call predict_from_fitted(fit,series,5,forecast,error)
  if(error%code/=arfima_ok) error stop trim(error%message)

  write(*,'(a,f10.5)') 'true phi:       ',truth%phi(1)
  write(*,'(a,f10.5)') 'estimated phi:  ',fit%parameters%phi(1)
  write(*,'(a,f10.5)') 'true d:         ',truth%dfrac
  write(*,'(a,f10.5)') 'estimated d:    ',fit%parameters%dfrac
  write(*,'(a,f10.5)') 'estimated mean: ',fit%parameters%mean
  write(*,'(a,f12.5)') 'log likelihood: ',fit%loglik
  write(*,'(a)') 'five-step forecasts:'
  do i=1,5
    write(*,'(i3,2f14.6)') i,forecast%mean(i),forecast%standard_error(i)
  end do
end program arfima_demo
