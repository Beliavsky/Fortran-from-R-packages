! SPDX-License-Identifier: GPL-2.0-or-later
program random_example
   use gnorm
   implicit none
   real(dp), allocatable :: x(:)
   real(dp) :: average, variance

   x = rgnorm(10000, mu=0.0_dp, alpha=1.0_dp, beta=1.5_dp, seed=42_i8)
   average = sum(x) / real(size(x), dp)
   variance = sum((x - average)**2) / real(size(x) - 1, dp)
   print '(a,f12.6)', 'sample mean:     ', average
   print '(a,f12.6)', 'sample variance: ', variance
   print '(a,f12.6)', 'theory variance: ', gnorm_variance(beta=1.5_dp)
end program random_example
