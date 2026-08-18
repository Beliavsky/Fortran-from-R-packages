! SPDX-License-Identifier: GPL-2.0-only
program test_quantile
   use poibin, only : dp, qpoibin, dpoibin, ppoibin
   implicit none
   real(dp), parameter :: pp(3) = [0.2_dp, 0.5_dp, 0.8_dp]
   integer :: k
   real(dp) :: q

   if (abs(dpoibin(0, pp) - 0.08_dp) > 1.0e-13_dp) error stop 1
   if (abs(dpoibin(1, pp) - 0.42_dp) > 1.0e-13_dp) error stop 2
   if (abs(dpoibin(2, pp) - 0.42_dp) > 1.0e-13_dp) error stop 3
   if (abs(dpoibin(3, pp) - 0.08_dp) > 1.0e-13_dp) error stop 4

   do k = 0, 3
      q = ppoibin(k, pp)
      if (qpoibin(q, pp) > k) error stop 5
      if (k > 0) then
         if (qpoibin(max(0.0_dp, ppoibin(k - 1, pp) + 1.0e-10_dp), pp) /= k) error stop 6
      end if
   end do
   if (qpoibin(0.0_dp, pp) /= 0) error stop 7
   if (qpoibin(1.0_dp, pp) /= 3) error stop 8
   print '(a)', 'test_quantile: PASS'
end program test_quantile
