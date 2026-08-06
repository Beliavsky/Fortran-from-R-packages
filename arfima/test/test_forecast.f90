program test_forecast
  use arfima
  use test_support
  implicit none
  real(dp),allocatable::z(:),y(:),roundtrip(:),w(:)
  type(arfima_error)::err
  type(arfima_spec)::spec
  type(arfima_parameters)::par
  type(arfima_forecast_result)::fc
  integer::h,i
  real(dp)::phi,varh

  z=[1.0_dp,1.3_dp,1.1_dp,1.8_dp,2.2_dp,2.0_dp,2.7_dp,3.1_dp,2.9_dp,3.4_dp]
  call difference_series(z,1,1,3,y,err)
  call integrate_series(y,z(1:4),1,1,3,roundtrip,err)
  call assert_vector_close(roundtrip,z,1.0e-12_dp,'difference/integration round trip')
  w=exact_integration_weights(1,0,0,5)
  call assert_vector_close(w,[1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp],1.0e-14_dp,'I(1) integration weights')

  phi=0.6_dp
  spec%p=1; spec%lmodel=long_memory_none
  allocate(par%phi(1),par%theta(0),par%phiseas(0),par%thetaseas(0)); par%phi=[phi]; par%mean=0.0_dp
  y=[0.1_dp,-0.2_dp,0.4_dp,0.3_dp]
  call arfima_forecast(spec,par,y,3,fc,err)
  call assert_true(err%code==arfima_ok,'stationary forecast status')
  do h=1,3
    call assert_close(fc%mean(h),phi**h*y(size(y)),2.0e-10_dp,'AR(1) forecast mean')
    varh=sum([(phi**(2*(i-1)),i=1,h)])
    call assert_close(fc%covariance(h,h),varh,2.0e-9_dp,'AR(1) forecast variance')
  end do

  spec%dint=1
  z=[0.0_dp,0.2_dp,0.1_dp,0.5_dp,0.8_dp]
  call arfima_forecast(spec,par,z,1,fc,err)
  call assert_close(fc%mean(1),z(5)+phi*(z(5)-z(4)),2.0e-9_dp,'ARIMA(1,1,0) first forecast')
  write(*,'(a)') 'test_forecast: PASS'
end program test_forecast
