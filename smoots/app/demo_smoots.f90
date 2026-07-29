! SPDX-License-Identifier: GPL-3.0-only
program demo_smoots
   use smoots
   implicit none
   integer, parameter :: n = 240
   real(dp) :: y(n), t, pi
   integer :: i
   type(smooth_result) :: fit
   type(forecast_result) :: fc
   pi = acos(-1.0_dp)
   call seed_rng(1234567_8)
   do i=1,n
      t=real(i,dp)/real(n,dp)
      y(i)=2.0_dp+1.5_dp*t+0.35_dp*sin(4.0_dp*pi*t)+0.15_dp*sin(37.0_dp*t)
   end do
   call msmooth(y,fit,p=1,mu=1,algorithm='A')
   call model_forecast(fit,1,0,5,.false.,0.95_dp,fc)
   print '(a,f10.6)', 'selected bandwidth: ',fit%b0
   print '(a,5f12.6)', 'forecasts: ',fc%point
end program demo_smoots
