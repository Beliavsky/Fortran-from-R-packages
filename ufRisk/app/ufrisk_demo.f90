! SPDX-License-Identifier: GPL-3.0-only
program ufrisk_demo
   use ufrisk
   implicit none
   integer, parameter :: n = 180
   real(dp) :: prices(n),return_value
   integer :: i
   type(ufrisk_options) :: options
   type(ufrisk_result) :: forecast
   type(ufrisk_traffic_result) :: traffic

   prices(1) = 100.0_dp
   do i = 2,n
      return_value = 0.0003_dp+0.011_dp*sin(0.21_dp*real(i,dp)) + &
         0.004_dp*cos(0.073_dp*real(i*i,dp))
      prices(i) = prices(i-1)*exp(return_value)
   end do
   options%n_out = 20
   options%model = ufrisk_model_sgarch
   options%distribution = ufrisk_distribution_student
   options%max_fit_iterations = 80
   forecast = varcast(prices,options)
   if (forecast%status /= ufrisk_ok) then
      write(*,'(a)') trim(forecast%message)
      error stop 1
   end if
   traffic = trafftest(forecast)
   write(*,'(a,a)') 'model: ',trim(forecast%model_name)
   write(*,'(a,f10.4)') 'estimated Student-t df: ',forecast%degrees_freedom
   write(*,'(a,3f12.6)') 'first sigma, VaR, ES: ',forecast%sigma_forecast(1), &
      forecast%var_var_level(1),forecast%expected_shortfall(1)
   write(*,'(a,i0)') '99% VaR violations: ',traffic%violations_var
end program ufrisk_demo
