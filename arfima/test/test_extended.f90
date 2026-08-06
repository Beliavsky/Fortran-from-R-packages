program test_extended
  use arfima
  use test_support
  implicit none
  type(arfima_spec)::spec
  type(arfima_fit_result)::fit
  real(dp),allocatable::z(:)
  integer::i

  spec%lmodel=long_memory_none
  spec%estimate_mean=.false.
  spec%use_transfer=.true.
  allocate(spec%transfer%r(1),spec%transfer%s(1),spec%transfer%b(1))
  allocate(spec%transfer%delta(0),spec%transfer%omega(1),spec%transfer%x(80,1))
  spec%transfer%r=0
  spec%transfer%s=1
  spec%transfer%b=0
  spec%transfer%omega=0.0_dp
  do i=1,80
    spec%transfer%x(i,1)=sin(0.1_dp*real(i,dp))
  end do
  allocate(z(80))
  z=1.8_dp*spec%transfer%x(:,1)
  do i=1,80
    z(i)=z(i)+0.05_dp*cos(0.7_dp*real(i,dp))
  end do
  call fit_arfima(spec,z,fit,max_iterations=1000,tolerance=2.0e-7_dp)
  call assert_true(fit%error%code==arfima_ok .or. fit%error%code==arfima_no_convergence,'transfer fit status')
  call assert_close(fit%parameters%omega(1),1.8_dp,0.12_dp,'transfer coefficient recovery')

  spec%use_transfer=.false.
  spec%use_regression=.true.
  allocate(spec%xreg(80,2))
  do i=1,80
    spec%xreg(i,1)=1.0_dp+0.02_dp*real(i,dp)
    spec%xreg(i,2)=cos(0.13_dp*real(i,dp))
  end do
  z=0.7_dp*spec%xreg(:,1)-1.2_dp*spec%xreg(:,2)+0.03_dp*[(sin(0.8_dp*real(i,dp)),i=1,80)]
  call fit_arfima(spec,z,fit,max_iterations=1000,tolerance=2.0e-7_dp)
  call assert_close(fit%parameters%beta(1),0.7_dp,0.08_dp,'regression beta 1')
  call assert_close(fit%parameters%beta(2),-1.2_dp,0.08_dp,'regression beta 2')
  write(*,'(a)') 'test_extended: PASS'
end program test_extended
