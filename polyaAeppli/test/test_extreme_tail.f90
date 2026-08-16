program test_extreme_tail
   use polya_aeppli
   implicit none
   integer :: fails,q
   real(dp) :: llo,lhi,lsf,p

   fails = 0

   ! Mirrors the extreme-tail regime in the upstream package examples.
   llo = p_polya_aeppli(3000.0_dp,4000.0_dp,0.005_dp,.true.,.true.)
   lhi = p_polya_aeppli(5000.0_dp,4000.0_dp,0.005_dp,.false.,.true.)
   if (.not. (llo < -100.0_dp)) fails=fails+1
   if (.not. (lhi < -50.0_dp)) fails=fails+1

   p = 1.0e-12_dp
   q = q_polya_aeppli(p,4000.0_dp,0.005_dp,.false.,.false.)
   lsf = p_polya_aeppli(real(q,dp),4000.0_dp,0.005_dp,.false.,.true.)
   if (lsf > log(p)+2.0e-12_dp) fails=fails+1
   if (q > 0) then
      lsf = p_polya_aeppli(real(q-1,dp),4000.0_dp,0.005_dp,.false.,.true.)
      if (lsf <= log(p)-2.0e-12_dp) fails=fails+1
   end if

   if (fails /= 0) then
      print *, "test_extreme_tail: FAIL", fails, llo, lhi, q
      error stop 1
   end if
   print *, "test_extreme_tail: PASS"
end program test_extreme_tail
