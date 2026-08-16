program test_countermeasure
   use genbinomapps
   implicit none
   integer(i64) :: sizes(2),nreq
   real(dp) :: effect(2)
   type(confidence_interval) :: ci
   integer :: fails

   fails = 0
   sizes = [1_i64,1_i64]
   effect = [0.5_dp,0.8_dp]

   ci = cm_clopper_pearson_ci(110000_i64,sizes,effect,0.1_dp,"upper")
   if (ci%status /= 0) fails=fails+1
   if (abs(ci%upper-3.320868588998336e-5_dp) > 2.0e-12_dp) fails=fails+1

   ci = cm_clopper_pearson_ci(110000_i64,sizes,effect,0.1_dp,"two.sided")
   if (ci%status /= 0) fails=fails+1
   if (abs(ci%upper-4.144537100217226e-5_dp) > 2.0e-12_dp) fails=fails+1

   nreq = cm_n_clopper_pearson(1.0e-5_dp,sizes,effect,0.1_dp)
   if (nreq /= 365299_i64) fails=fails+1

   if (fails /= 0) then
      print *, "test_countermeasure: FAIL", fails, ci%lower, ci%upper, nreq
      error stop 1
   end if
   print *, "test_countermeasure: PASS"
end program test_countermeasure
