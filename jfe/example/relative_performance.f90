! SPDX-License-Identifier: GPL-2.0-or-later
program relative_performance
   use jfe
   implicit none

   real(dp), parameter :: fund(8) = [0.014_dp, -0.006_dp, 0.018_dp, 0.003_dp, &
      -0.009_dp, 0.011_dp, 0.007_dp, -0.001_dp]
   real(dp), parameter :: benchmark(8) = [0.010_dp, -0.005_dp, 0.012_dp, 0.002_dp, &
      -0.008_dp, 0.008_dp, 0.005_dp, -0.002_dp]

   print '(a,f12.6)', 'Active premium:    ', active_premium(fund, benchmark, 12.0_dp)
   print '(a,f12.6)', 'Tracking error:    ', tracking_error(fund, benchmark, 12.0_dp)
   print '(a,f12.6)', 'Information ratio: ', information_ratio(fund, benchmark, 12.0_dp)
   print '(a,f12.6)', 'Jensen alpha:      ', capm_jensen_alpha(fund, benchmark, scale=12.0_dp)
   print '(a,f12.6)', 'Treynor ratio:     ', treynor_ratio(fund, benchmark, scale=12.0_dp)
end program relative_performance
