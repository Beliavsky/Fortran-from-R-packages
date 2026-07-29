! SPDX-License-Identifier: GPL-3.0-only
program semiparametric_figarch
   use ufrisk
   implicit none
   integer, parameter :: n = 160
   real(dp) :: prices(n),return_value,scale_value
   integer :: i
   type(ufrisk_options) :: options
   type(ufrisk_result) :: result

   prices(1) = 100.0_dp
   do i = 2,n
      scale_value = 0.007_dp+0.00003_dp*real(i,dp)
      return_value = scale_value*(sin(0.34_dp*real(i,dp))+0.35_dp*cos(0.17_dp*real(i,dp)))
      prices(i) = prices(i-1)*exp(return_value)
   end do
   options%model = ufrisk_model_figarch
   options%smooth = ufrisk_smooth_lpr
   options%distribution = ufrisk_distribution_normal
   options%n_out = 15
   options%smoothing_iterations = 6
   options%max_fit_iterations = 60
   result = varcast(prices,options)
   if (result%status /= ufrisk_ok) error stop trim(result%message)
   write(*,'(a,f10.6)') 'selected scale bandwidth: ',result%long_memory_smooth%bandwidth
   write(*,'(a,f10.6)') 'estimated scale-memory d: ',result%long_memory_smooth%d
   write(*,'(a,3f12.6)') 'first forecast sigma, VaR, ES: ',result%sigma_forecast(1), &
      result%var_var_level(1),result%expected_shortfall(1)
end program semiparametric_figarch
