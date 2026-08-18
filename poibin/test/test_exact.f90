! SPDX-License-Identifier: GPL-2.0-only
program test_exact
   use poibin, only : dp, poibin_pmf_dft_cf, poibin_pmf_rf, ppoibin
   implicit none
   real(dp), parameter :: pp(5) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
   integer, parameter :: wts(5) = [2, 2, 2, 2, 2]
   real(dp), parameter :: ref(0:10) = [ &
      0.02286144_dp, 0.11231136_dp, 0.23694372_dp, 0.28171848_dp, &
      0.20829484_dp, 0.09968032_dp, 0.03113484_dp, 0.00623848_dp, &
      0.00076372_dp, 0.00005136_dp, 0.00000144_dp]
   real(dp), allocatable :: a(:), b(:)
   integer :: status

   call poibin_pmf_dft_cf(pp, a, wts, status)
   if (status /= 0) error stop 1
   call poibin_pmf_rf(pp, b, wts, status)
   if (maxval(abs(a - ref)) > 2.0e-13_dp) error stop 2
   if (maxval(abs(b - ref)) > 2.0e-15_dp) error stop 3
   if (maxval(abs(a - b)) > 2.0e-13_dp) error stop 4
   if (abs(ppoibin(3, pp, 'DFT-CF', wts) - 0.653835_dp) > 2.0e-13_dp) error stop 5
   if (abs(ppoibin(3, pp, 'RF', wts) - 0.653835_dp) > 2.0e-15_dp) error stop 6
   print '(a)', 'test_exact: PASS'
end program test_exact
