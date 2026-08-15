! SPDX-License-Identifier: LGPL-3.0-only
program basic
   use distr
   implicit none

   type(distribution_t) :: x, y, z
   real(dp), allocatable :: sample(:)

   x = normal_dist(mean=1.0_dp, sd=2.0_dp)
   y = exponential_dist(rate=0.5_dp)
   z = x + y

   print '(a,f12.6)', 'normal density at zero: ', x%density(0.0_dp)
   print '(a,f12.6)', 'normal CDF at zero:     ', x%cdf(0.0_dp)
   print '(a,f12.6)', 'normal 95% quantile:    ', x%quantile(0.95_dp)
   print '(a,f12.6)', 'sum mean:               ', z%mean()
   print '(a,f12.6)', 'sum standard deviation: ', z%sd()

   call seed_rng(12345)
   sample = z%random(10000)
   print '(a,f12.6)', 'sample mean:             ', &
        sum(sample)/real(size(sample), dp)
end program basic
