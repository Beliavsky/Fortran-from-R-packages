! SPDX-License-Identifier: GPL-2.0-or-later
program estimate_and_forecast
   use garchsk, only : dp, estimate_result, forecast_result, garchsk_est, garchsk_fcst
   implicit none
   real(dp) :: x(80)
   type(estimate_result) :: fit
   type(forecast_result) :: fcst
   integer :: i
   do i = 1, size(x)
      x(i) = 0.01_dp*sin(0.23_dp*real(i, dp)) + 0.003_dp*cos(0.61_dp*real(i, dp))
   end do
   fit = garchsk_est(x, max_iterations=1200, tolerance=1.0e-5_dp)
   fcst = garchsk_fcst(fit%params, x, 3)
   print '(a,l1)', 'converged: ', fit%converged
   print '(a,es14.6)', 'negative log likelihood: ', fit%negative_log_likelihood
   print '(a,3(1x,es12.4))', 'variance forecasts:', fcst%h
end program estimate_and_forecast
