! SPDX-License-Identifier: GPL-2.0-only
program test_rng
   use poibin, only : dp, rpoibin_sample, poibin_seed, poibin_mean
   implicit none
   integer, parameter :: m = 50000
   real(dp), parameter :: pp(5) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
   integer, parameter :: wts(5) = [2, 2, 2, 2, 2]
   integer :: x(m)
   real(dp) :: empirical, truth

   call poibin_seed(12345)
   call rpoibin_sample(m, pp, x, wts)
   empirical = sum(real(x, dp))/real(m, dp)
   truth = poibin_mean(pp, wts)
   if (abs(empirical - truth) > 0.03_dp) error stop 1
   if (minval(x) < 0 .or. maxval(x) > 10) error stop 2
   print '(a)', 'test_rng: PASS'
end program test_rng
