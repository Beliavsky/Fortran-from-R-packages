! SPDX-License-Identifier: GPL-2.0-only
program test_approximations
   use poibin, only : dp, ppoibin
   implicit none
   real(dp), parameter :: pp(5) = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp]
   integer, parameter :: wts(5) = [2, 2, 2, 2, 2]
   real(dp) :: pna, prna, ppa

   pna = ppoibin(3, pp, 'NA', wts)
   prna = ppoibin(3, pp, 'RNA', wts)
   ppa = ppoibin(3, pp, 'PA', wts)
   if (abs(pna - 0.6415997415027527_dp) > 2.0e-13_dp) error stop 1
   if (abs(prna - 0.6539859390494451_dp) > 2.0e-13_dp) error stop 2
   if (abs(ppa - 0.6472318887822313_dp) > 2.0e-13_dp) error stop 3
   print '(a)', 'test_approximations: PASS'
end program test_approximations
