program test_autocov
  use arfima
  use test_support
  implicit none
  real(dp),allocatable::p(:),a(:),g(:)
  type(arfima_error)::err
  type(arfima_spec)::spec
  type(arfima_parameters)::par
  integer::k

  a=[0.30_dp,-0.20_dp,0.10_dp]
  p=ar_to_pacf(a)
  call assert_vector_close(pacf_to_ar(p),a,1.0e-12_dp,'AR/PACF round trip')
  call assert_true(is_stationary_polynomial(a),'stable AR polynomial')

  call tacvf_arma([0.5_dp],[real(dp)::],5,1.0_dp,g,err)
  call assert_true(err%code==arfima_ok,'AR(1) autocovariance status')
  do k=0,5
    call assert_close(g(k+1),(0.5_dp**k)/(1.0_dp-0.25_dp),2.0e-12_dp,'AR(1) autocovariance')
  end do
  call tacvf_arma([real(dp)::],[0.4_dp],3,1.0_dp,g,err)
  call assert_vector_close(g,[1.16_dp,-0.4_dp,0.0_dp,0.0_dp],2.0e-12_dp,'MA(1) autocovariance')

  call tacvf_fdwn(0.2_dp,4,g,err)
  call assert_close(g(2),0.2_dp/0.8_dp*g(1),1.0e-13_dp,'FDWN recurrence')
  call tacvf_fgn(0.5_dp,4,g,err)
  call assert_vector_close(g,[1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],2.0e-13_dp,'FGN H=0.5')
  call tacvf_pla(1.5_dp,3,g,err)
  call assert_close(g(1),1.0_dp,0.0_dp,'PLA variance')

  spec%period=4; spec%pseas=1; spec%slmodel=long_memory_none; spec%lmodel=long_memory_none
  allocate(par%phi(0),par%theta(0),par%phiseas(1),par%thetaseas(0))
  par%phiseas=[0.4_dp]
  call tacvf_arfima(spec,par,8,1.0_dp,g,err)
  call assert_true(err%code==arfima_ok,'seasonal autocovariance status')
  call assert_close(g(1),1.0_dp/(1.0_dp-0.16_dp),2.0e-10_dp,'seasonal AR variance')
  call assert_close(g(5),0.4_dp*g(1),2.0e-10_dp,'seasonal AR lag period')
  call assert_close(g(2),0.0_dp,2.0e-10_dp,'seasonal AR off-season lag')
  write(*,'(a)') 'test_autocov: PASS'
end program test_autocov
