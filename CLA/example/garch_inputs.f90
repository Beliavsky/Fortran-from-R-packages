! SPDX-License-Identifier: GPL-3.0-or-later
program garch_inputs
   use kind_mod, only: dp
   use cla, only: cla_garch_result_t, mu_sigma_garch, cla_distribution_normal
   implicit none
   integer,parameter::n=80
   real(dp)::prices(n,2)
   type(cla_garch_result_t)::inputs
   integer::t
   prices(1,:)=[100.0_dp,90.0_dp]
   do t=2,n
      prices(t,1)=prices(t-1,1)*exp(0.0003_dp+0.01_dp*sin(0.2_dp*real(t,dp)))
      prices(t,2)=prices(t-1,2)*exp(0.0001_dp+0.008_dp*cos(0.17_dp*real(t,dp)))
   end do
   inputs=mu_sigma_garch(prices,distribution=cla_distribution_normal,max_iterations=100)
   if(inputs%info/=0)error stop 'GARCH estimation failed'
   write(*,'(a,2f12.7)')'Forecast means: ',inputs%mu
   write(*,'(a,2f12.7)')'Forecast sigmas:',inputs%forecast_sigma
   write(*,'(a)')'Forecast covariance:'
   write(*,'(2f14.9)')inputs%covariance(1,:)
   write(*,'(2f14.9)')inputs%covariance(2,:)
end program garch_inputs
