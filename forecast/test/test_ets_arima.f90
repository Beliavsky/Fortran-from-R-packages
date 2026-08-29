program test_ets_arima
   use forecast, only : dp, ets_model, arima_model, forecast_result, ets_fit, ets_forecast, ets_auto, &
      arima_fit, arima_forecast, auto_arima, ETS_ADD, ETS_NONE
   implicit none
   real(dp) :: y(60)
   integer :: i
   type(ets_model) :: em
   type(arima_model) :: am
   type(forecast_result) :: fc
   do i=1,60
   y(i)=10.0_dp+0.2_dp*real(i,dp)+0.5_dp*sin(0.4_dp*real(i,dp))
   end do
   em=ets_fit(y,1,ETS_ADD,ETS_ADD,ETS_NONE,.false.,.true.)
   fc=ets_forecast(em,5)
   call check(size(fc%mean)==5 .and. all(fc%mean==fc%mean),'ETS forecast finite')
   call check(em%sigma2>=0.0_dp,'ETS variance')
   am=arima_fit(y,1,1,0,m=1,include_mean=.false.,optimize=.true.)
   fc=arima_forecast(am,y,5)
   call check(size(fc%mean)==5 .and. all(fc%mean==fc%mean),'ARIMA forecast finite')
   am=auto_arima(y,m=1,max_p=2,max_q=1,seasonal=.false.,stepwise=.true.)
   call check(am%aicc<huge(1.0_dp),'auto_arima finite criterion')
   print '(a)','test_ets_arima: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      write(*,'(a)')'FAIL: '//trim(msg)
      error stop 1
      end if
   end subroutine
end program
