program test_tails
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use bridgedist, only : dp, dbridge, pbridge, qbridge
   implicit none
   real(dp) :: p1, p2, q1
   integer :: fails

   fails = 0
   p1 = pbridge(1000.0_dp, 0.9_dp, .false.)
   p2 = pbridge(-1000.0_dp, 0.9_dp)
   if (.not. ieee_is_finite(p1) .or. p1 < 0.0_dp .or. p1 > 1.0_dp) fails = fails + 1
   if (abs(p1 - p2) > tiny(1.0_dp)) fails = fails + 1
   if (.not. ieee_is_finite(dbridge(1000.0_dp, 0.9_dp, .true.))) fails = fails + 1
   q1 = qbridge(1.0e-300_dp, 0.6_dp)
   if (.not. ieee_is_finite(q1) .or. q1 >= 0.0_dp) fails = fails + 1
   if (abs(pbridge(q1, 0.6_dp) / 1.0e-300_dp - 1.0_dp) > 5.0e-13_dp) fails = fails + 1

   if (fails /= 0) then
      print '(a,i0)', 'test_tails: FAIL ', fails
      error stop 1
   end if
   print '(a)', 'test_tails: PASS'
end program test_tails
