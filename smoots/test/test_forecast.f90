! SPDX-License-Identifier: GPL-3.0-only
program test_forecast
   use smoots
   implicit none
   integer,parameter::n=180
   real(dp)::y(n),x
   integer::i
   type(smooth_result)::fit
   type(forecast_result)::nf,bf,mf
   type(rolling_result)::roll
   do i=1,n
      x=real(i,dp)/real(n,dp)
      y(i)=3.0_dp+x+0.25_dp*sin(0.4_dp*real(i,dp))+0.08_dp*cos(1.7_dp*real(i,dp))
   end do
   call msmooth(y,fit)
   call normal_forecast(fit%residuals,1,1,.false.,5,0.95_dp,nf)
   if(nf%status/=sm_ok.or.any(nf%lower>nf%upper))error stop 'normal forecast'
   call bootstrap_forecast(fit%residuals,1,0,.false.,3,80,100,0.90_dp,bf,987654321_8)
   if(bf%status/=sm_ok.or.size(bf%errors,1)/=80)error stop 'bootstrap forecast'
   call model_forecast(fit,1,0,3,.false.,0.95_dp,mf)
   if(mf%status/=sm_ok)error stop 'model forecast'
   call rolling_backtest(y,5,1,0,.false.,0.95_dp,roll)
   if(roll%status/=sm_ok.or.roll%mase<0.0_dp.or.roll%rmsse<0.0_dp)error stop 'rolling backtest'
   print '(a)','test_forecast: PASS'
end program test_forecast
