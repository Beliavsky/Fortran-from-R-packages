! SPDX-License-Identifier: GPL-2.0-or-later
program lps_demo
   use lindley_power_series
   implicit none

   real(dp), parameter :: lambda = 1.0_dp, theta = 0.5_dp
   real(dp) :: x, p

   x = 1.25_dp
   p = 0.75_dp
   print '(a,f12.8)', 'Lindley-geometric CDF:      ', &
      plindleygeometric(x,lambda,theta)
   print '(a,f12.8)', 'Lindley-logarithmic PDF:   ', &
      dlindleylogarithmic(x,lambda,theta)
   print '(a,f12.8)', 'Lindley-Poisson 75% q:     ', &
      qlindleypoisson(p,lambda,theta)
   print '(a,f12.8)', 'Lindley-binomial hazard:   ', &
      hlindleybinomial(x,lambda,theta,4)
end program lps_demo
