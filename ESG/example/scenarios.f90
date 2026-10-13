! SPDX-License-Identifier: GPL-2.0-or-later
program scenarios
   use esg
   implicit none
   type(esg_parameters) :: parameters
   type(esg_scenario) :: paths
   type(randompack_rng_type) :: stream
   real(dp) :: zero(73), normals(1000,4,5), deviation
   real(dp), allocatable :: prices(:), means(:)
   integer :: status
   zero = 0.03_dp
   stream = randompack_rng('pcg64', seed=20261012)
   call esg_innovations(stream, normals, status)
   if (status /= 0) error stop 'innovation generation failed'
   call esg_generate(parameters, zero, normals, paths, status)
   if (status /= 0) error stop 'scenario generation failed'
   call esg_zc_bond(paths, 1, 3.5_dp, prices, status)
   if (status /= 0) error stop 'zero-coupon valuation failed'
   print '(a,f12.8)', 'Mean year-one zero-coupon price: ', sum(prices) / real(size(prices), dp)
   call esg_stock_call(paths, 1, 3.5_dp, 100.0_dp, prices, status)
   if (status /= 0) error stop 'stock option valuation failed'
   print '(a,f12.8)', 'Mean source-compatible stock call: ', sum(prices) / real(size(prices), dp)
   call esg_martingale(paths, means, deviation, status)
   if (status /= 0) error stop 'martingale diagnostic failed'
   print '(a,f12.8)', 'Maximum discounted-stock mean deviation: ', deviation
end program
