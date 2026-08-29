program demo_forecast
   use forecast, only : dp, ets_model, forecast_result, ets_auto, ets_forecast
   implicit none
   real(dp) :: y(48)
   integer :: i
   type(ets_model) :: model
   type(forecast_result) :: fc
   do i=1,size(y)
      y(i)=100.0_dp+0.4_dp*i+8.0_dp*sin(2.0_dp*acos(-1.0_dp)*i/12.0_dp)
   end do
   model=ets_auto(y,12)
   fc=ets_forecast(model,6,[80.0_dp,95.0_dp])
   print '(a,6f11.3)','Forecasts:',fc%mean
end program demo_forecast
