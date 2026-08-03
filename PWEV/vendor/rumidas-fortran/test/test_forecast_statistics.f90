program test_forecast_statistics
  use rumidas
  implicit none
  real(dp) :: p(5),aic,bic,mse,qlike,v(3),proxy(3),expected
  real(dp),allocatable::forecast(:)
  type(garch_midas_spec)::spec
  integer::status

  call information_criteria(-100.0_dp,5,200,aic,bic)
  call check_close(aic,210.0_dp,1.0e-13_dp,'AIC')
  call check_close(bic,5.0_dp*log(200.0_dp)+200.0_dp,1.0e-13_dp,'BIC')
  v=[1.0_dp,2.0_dp,4.0_dp];proxy=[1.1_dp,1.8_dp,4.2_dp]
  call volatility_loss_functions(v,proxy,mse,qlike)
  call check(mse>0.0_dp.and.qlike>0.0_dp,'losses')

  p=[0.1_dp,0.8_dp,-2.0_dp,0.1_dp,2.0_dp]
  spec=garch_midas_spec(RUMIDAS_GM,RUMIDAS_NORMAL,RUMIDAS_BETA_LAG,2,0,.false.)
  call multi_step_ahead_pred(p,spec,1.4_dp,0.04_dp,4,forecast,status)
  call check(status==0,'forecast status')
  expected=(1.0_dp+(p(1)+p(2))**2*(1.4_dp-1.0_dp))*0.04_dp
  call check_close(forecast(3),expected,1.0e-13_dp,'forecast recursion')
  call check(all(forecast>0.0_dp),'forecast positive')

  print '(a)', 'test_forecast_statistics: PASS'
contains
  subroutine check(condition,message)
    logical,intent(in)::condition
    character(len=*),intent(in)::message
    if(.not.condition) error stop message
  end subroutine check
  subroutine check_close(actual,expected,tol,message)
    real(dp),intent(in)::actual,expected,tol
    character(len=*),intent(in)::message
    call check(abs(actual-expected)<=tol,message)
  end subroutine check_close
end program test_forecast_statistics
