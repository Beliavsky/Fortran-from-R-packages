program test_regression_spline
   use forecast, only : dp, regression_model, spline_model_t, linear_ar_model, forecast_result, &
      tslm_fit, regression_forecast, trend_season_matrix, modelar_fit, modelar_forecast, spline_model_fit, spline_forecast
   implicit none
   real(dp)::y(50)
   real(dp),allocatable::x(:,:),xn(:,:)
   integer::i
   type(regression_model)::lm
   type(linear_ar_model)::ar
   type(spline_model_t)::sp
   type(forecast_result)::fc
   do i=1,50
   y(i)=3.0_dp+0.5_dp*i+sin(0.2_dp*i)
   end do
   x=trend_season_matrix(50,1,include_season=.false.)
   lm=tslm_fit(y,x,.true.)
   call check(abs(lm%coefficients(2)-0.5_dp)<0.03_dp,'linear trend coefficient')
   xn=trend_season_matrix(3,1,start_index=51,include_season=.false.)
   fc=regression_forecast(lm,xn)
   call check(size(fc%mean)==3 .and. all(fc%mean==fc%mean),'regression forecast')
   ar=modelar_fit(y,p=2)
   fc=modelar_forecast(ar,y,3)
   call check(all(fc%mean==fc%mean),'modelAR')
   sp=spline_model_fit(y)
   fc=spline_forecast(sp,3,[80.0_dp,95.0_dp])
   call check(size(fc%mean)==3 .and. all(fc%se>=0.0_dp),'spline forecast')
   print '(a)','test_regression_spline: PASS'
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
