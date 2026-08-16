program test_reference
   use polya_aeppli
   implicit none
   integer :: fails
   real(dp) :: p

   fails = 0

   ! Independent high-precision references from the compound-Poisson sum.
   if (abs(d_polya_aeppli(0.0_dp,2.5_dp,0.3_dp)- &
           8.208499862389880e-2_dp) > 2.0e-14_dp) fails=fails+1
   if (abs(d_polya_aeppli(1.0_dp,2.5_dp,0.3_dp)- &
           1.436487475918229e-1_dp) > 2.0e-14_dp) fails=fails+1
   if (abs(d_polya_aeppli(5.0_dp,2.5_dp,0.3_dp)- &
           1.040521573607770e-1_dp) > 2.0e-14_dp) fails=fails+1
   if (abs(d_polya_aeppli(12.0_dp,2.5_dp,0.3_dp)- &
           3.934601727939026e-3_dp) > 2.0e-14_dp) fails=fails+1

   p = p_polya_aeppli(5.0_dp,2.5_dp,0.3_dp)
   if (abs(p-7.961198656871009e-1_dp) > 3.0e-14_dp) fails=fails+1

   p = p_polya_aeppli(12.0_dp,2.5_dp,0.3_dp,.false.)
   if (abs(p-4.668840127275420e-3_dp) > 3.0e-14_dp) fails=fails+1


   ! prob=0 reduces exactly to Poisson(lambda).
   if (abs(d_polya_aeppli(0.0_dp,2.5_dp,0.0_dp)-exp(-2.5_dp)) > 2.0e-14_dp) fails=fails+1
   if (abs(d_polya_aeppli(1.0_dp,2.5_dp,0.0_dp)-2.5_dp*exp(-2.5_dp)) > 2.0e-14_dp) fails=fails+1

   if (q_polya_aeppli(0.1_dp,2.5_dp,0.3_dp) /= 1) fails=fails+1
   if (q_polya_aeppli(0.5_dp,2.5_dp,0.3_dp) /= 3) fails=fails+1
   if (q_polya_aeppli(0.9_dp,2.5_dp,0.3_dp) /= 7) fails=fails+1
   if (q_polya_aeppli(0.99_dp,2.5_dp,0.3_dp) /= 11) fails=fails+1

   if (fails /= 0) then
      print *, "test_reference: FAIL", fails
      error stop 1
   end if
   print *, "test_reference: PASS"
end program test_reference
