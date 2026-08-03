! SPDX-License-Identifier: GPL-2.0-or-later
program time_series_diagnostics
   use fints
   implicit none

   integer, parameter :: n = 250
   real(dp) :: x(n)
   type(acf_result) :: correlations, partial
   type(test_result) :: box, arch
   integer :: i

   x(1) = 0.0_dp
   do i = 2, n
      x(i) = 0.7_dp * x(i - 1) + 0.2_dp * sin(0.43_dp * real(i * i, dp))
   end do

   call acf(x, correlations, lag_max=10)
   call acf(x, partial, lag_max=10, acf_type='partial')
   call AutocorTest(x, box, lag=10, test_type='Ljung-Box')
   call ArchTest(x, arch, lags=8, demean=.true.)

   print '(a,10f9.4)', 'ACF:  ', correlations%value(2:11)
   print '(a,10f9.4)', 'PACF: ', partial%value
   print '(a,f12.6)', 'Ljung-Box p-value: ', box%p_value
   print '(a,f12.6)', 'ARCH LM p-value:   ', arch%p_value
end program time_series_diagnostics
