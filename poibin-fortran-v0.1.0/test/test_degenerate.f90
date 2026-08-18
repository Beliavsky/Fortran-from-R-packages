! SPDX-License-Identifier: GPL-2.0-only
program test_degenerate
   use poibin, only : dp, dpoibin, ppoibin, qpoibin
   implicit none
   real(dp), parameter :: pp(4) = [1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]

   if (abs(dpoibin(2, pp) - 1.0_dp) > 1.0e-13_dp) error stop 1
   if (dpoibin(1, pp) > 1.0e-13_dp) error stop 2
   if (abs(ppoibin(1, pp, 'NA')) > 1.0e-15_dp) error stop 3
   if (abs(ppoibin(2, pp, 'RNA') - 1.0_dp) > 1.0e-15_dp) error stop 4
   if (qpoibin(0.5_dp, pp) /= 2) error stop 5
   print '(a)', 'test_degenerate: PASS'
end program test_degenerate
