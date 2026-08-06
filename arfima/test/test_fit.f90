program test_fit
  use arfima
  use test_support
  implicit none
  type(arfima_spec)::spec
  type(arfima_parameters)::par
  type(arfima_fit_result)::fit
  type(arfima_error)::err
  real(dp),allocatable::z(:)

  call set_random_seed(24680)
  spec%p=1; spec%lmodel=long_memory_none; spec%estimate_mean=.true.
  allocate(par%phi(1),par%theta(0),par%phiseas(0),par%thetaseas(0),par%beta(0),par%delta(0),par%omega(0))
  par%phi=[0.55_dp]; par%mean=1.2_dp
  call arfima_simulate(spec,par,140,1.0_dp,z,err,center_sample=.false.)
  call assert_true(err%code==arfima_ok,'AR fit simulation')
  call fit_arfima(spec,z,fit,max_iterations=1200,tolerance=2.0e-6_dp)
  call assert_true(fit%error%code==arfima_ok .or. fit%error%code==arfima_no_convergence,'AR fit status')
  call assert_close(fit%parameters%phi(1),0.55_dp,0.16_dp,'AR coefficient recovery')
  call assert_close(fit%parameters%mean,1.2_dp,0.35_dp,'mean recovery')
  call assert_true(fit%sigma2>0.0_dp,'positive fitted variance')

  call set_random_seed(13579)
  spec%p=0; spec%lmodel=long_memory_fd; spec%estimate_mean=.true.
  deallocate(par%phi); allocate(par%phi(0)); par%dfrac=0.25_dp; par%mean=0.3_dp
  call arfima_simulate(spec,par,90,0.8_dp,z,err,center_sample=.false.)
  call fit_arfima(spec,z,fit,max_iterations=900,tolerance=4.0e-6_dp)
  call assert_true(fit%error%code==arfima_ok .or. fit%error%code==arfima_no_convergence,'FD fit status')
  call assert_close(fit%parameters%dfrac,0.25_dp,0.20_dp,'fractional d recovery')
  write(*,'(a)') 'test_fit: PASS'
end program test_fit
