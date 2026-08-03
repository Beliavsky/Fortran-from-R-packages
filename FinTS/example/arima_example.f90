! SPDX-License-Identifier: GPL-2.0-or-later
program arima_example
   use fints
   implicit none

   integer, parameter :: n = 240
   real(dp) :: x(n), innovation
   type(arima_result) :: fit
   integer :: i

   x(1) = 1.5_dp
   do i = 2, n
      innovation = 0.15_dp * sin(0.17_dp * real(i * i, dp))
      x(i) = 1.5_dp + 0.6_dp * (x(i - 1) - 1.5_dp) + innovation
   end do

   call ARIMA(x, [1, 0, 0], fit, include_mean=.true., &
      box_test_type='Ljung-Box')

   print '(a,f10.5)', 'estimated mean = ', fit%intercept
   print '(a,f10.5)', 'estimated AR(1) = ', fit%ar(1)
   print '(a,f10.6)', 'innovation variance = ', fit%sigma2
   print '(a,f10.6)', 'residual-test p-value = ', fit%box_test%p_value
end program arima_example
