! SPDX-License-Identifier: GPL-2.0-or-later
program demo_fints
   use fints
   implicit none

   integer, parameter :: n = 180
   real(dp) :: returns(n)
   type(summary_result) :: stats
   type(acf_result) :: correlations
   type(test_result) :: arch
   type(arima_result) :: fit
   integer :: i

   returns(1) = 0.01_dp
   do i = 2, n
      returns(i) = 0.45_dp * returns(i - 1) + &
         0.012_dp * sin(0.37_dp * real(i * i, dp))
   end do

   call FinTS_stats(returns, stats, start=2000.01_dp)
   call acf(returns, correlations, lag_max=5)
   call ArchTest(returns, arch, lags=5, demean=.true.)
   call ARIMA(returns, [1, 0, 0], fit, include_mean=.true.)

   print '(a,i0)', 'observations: ', stats%size
   print '(a,f12.8)', 'mean:         ', stats%mean
   print '(a,f12.8)', 'std. dev.:    ', stats%standard_deviation
   print '(a,5f10.5)', 'ACF 1:5:     ', correlations%value(2:6)
   print '(a,f10.5)', 'ARCH p-value: ', arch%p_value
   print '(a,f10.5)', 'AR(1):        ', fit%ar(1)
   print '(a,f10.5)', 'Box p-value:  ', fit%box_test%p_value
end program demo_fints
