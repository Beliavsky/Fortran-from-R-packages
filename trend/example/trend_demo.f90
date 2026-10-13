! SPDX-License-Identifier: GPL-3.0-only
program trend_demo
   use trend_module, only : dp, trend_result, trend_ok, mk_test, sens_slope, pettitt_test, br_test
   implicit none
   real(dp), parameter :: observations(12) = [8.0_dp, 3.0_dp, 5.0_dp, 1.0_dp, 9.0_dp, 2.0_dp, &
      7.0_dp, 4.0_dp, 12.0_dp, 10.0_dp, 6.0_dp, 11.0_dp]
   type(trend_result) :: test
   call mk_test(observations, test)
   if (test%status /= trend_ok) error stop 'MK test failed'
   print '(a,f9.5,a,f9.5)', 'Mann-Kendall z = ', test%statistic, ', p = ', test%p_value
   call sens_slope(observations, test)
   if (test%status /= trend_ok) error stop 'Sen slope failed'
   print '(a,f9.5,a,2f9.5)', 'Sen slope = ', test%estimate, ', 95% interval = ', test%conf_int
   call pettitt_test(observations, test)
   if (test%status /= trend_ok) error stop 'Pettitt test failed'
   print '(a,*(i0,1x))', 'Pettitt change-point candidates: ', test%change_points
   call br_test(observations, test, m = 20000, seed = 8675309)
   if (test%status /= trend_ok) error stop 'Buishand range test failed'
   print '(a,f9.5,a,f9.5)', 'Buishand R/sqrt(n) = ', test%statistic, ', Monte Carlo p = ', test%p_value
end program trend_demo
