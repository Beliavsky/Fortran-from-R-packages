! SPDX-License-Identifier: MIT
program distribution_functions
   use vasicekfit, only : dp, vasicek_density, vasicek_cdf, vasicek_quantile
   implicit none
   real(dp) :: x

   x = 0.05_dp
   print '(a,f12.8)', 'density at 5% loss = ', vasicek_density(x, 0.03_dp, 0.10_dp)
   print '(a,f12.8)', 'CDF at 5% loss     = ', vasicek_cdf(x, 0.03_dp, 0.10_dp)
   print '(a,f12.8)', '99% loss quantile  = ', vasicek_quantile(0.99_dp, 0.03_dp, 0.10_dp)
end program distribution_functions
