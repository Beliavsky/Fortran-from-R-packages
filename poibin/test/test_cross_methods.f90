! SPDX-License-Identifier: GPL-2.0-only
program test_cross_methods
   use poibin, only : dp, poibin_pmf_dft_cf, poibin_pmf_rf
   implicit none
   real(dp), parameter :: pp(8) = [ &
      0.01_dp, 0.07_dp, 0.19_dp, 0.33_dp, 0.51_dp, 0.72_dp, 0.91_dp, 0.99_dp]
   integer, parameter :: wts(8) = [3, 1, 4, 2, 5, 1, 2, 3]
   real(dp), allocatable :: a(:), b(:)
   integer :: status

   call poibin_pmf_dft_cf(pp, a, wts, status)
   if (status /= 0) error stop 1
   call poibin_pmf_rf(pp, b, wts, status)
   if (status /= 0) error stop 2
   if (abs(sum(a) - 1.0_dp) > 2.0e-14_dp) error stop 3
   if (abs(sum(b) - 1.0_dp) > 2.0e-14_dp) error stop 4
   if (maxval(abs(a - b)) > 3.0e-13_dp) error stop 5
   if (minval(a) < 0.0_dp .or. minval(b) < 0.0_dp) error stop 6
   print '(a)', 'test_cross_methods: PASS'
end program test_cross_methods
