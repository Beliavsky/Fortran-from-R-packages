program test_information_transfer
  use arfima
  use test_support
  implicit none
  type(arfima_spec)::spec
  type(arfima_parameters)::par
  type(arfima_error)::err
  type(transfer_spec)::tr
  real(dp),allocatable::info(:,:),res(:),eff(:),w(:)
  real(dp)::y(6),x(6,1)
  integer::i

  spec%p=1; spec%lmodel=long_memory_none
  allocate(par%phi(1),par%theta(0),par%phiseas(0),par%thetaseas(0)); par%phi=0.0_dp
  call arfima_information(spec,par,info,err,nfreq=512,maxlag=256)
  call assert_true(err%code==arfima_ok,'information matrix status')
  call assert_close(info(1,1),1.0_dp,2.0e-2_dp,'AR(1) information at zero')
  call assert_true(identifiable_invertible(spec,par,.true.,err),'identifiability check')

  do i=1,6; x(i,1)=real(i,dp); end do
  y=2.0_dp*x(:,1)+[1.0_dp,-1.0_dp,0.5_dp,0.0_dp,1.0_dp,-0.5_dp]
  allocate(tr%r(1),tr%s(1),tr%b(1),tr%delta(0),tr%omega(1),tr%x(6,1))
  tr%r=0; tr%s=1; tr%b=0; tr%omega=2.0_dp; tr%x=x
  call apply_transfer_function(y,tr,0.0_dp,res,eff,err)
  call assert_vector_close(res,[1.0_dp,-1.0_dp,0.5_dp,0.0_dp,1.0_dp,-0.5_dp],1.0e-12_dp,'static transfer residual')
  w=psi_weights([0.5_dp],[real(dp)::],[real(dp)::],[real(dp)::],0.0_dp,0.0_dp,0,0,0,5)
  call assert_vector_close(w,[1.0_dp,0.5_dp,0.25_dp,0.125_dp,0.0625_dp],1.0e-12_dp,'AR psi weights')
  write(*,'(a)') 'test_information_transfer: PASS'
end program test_information_transfer
