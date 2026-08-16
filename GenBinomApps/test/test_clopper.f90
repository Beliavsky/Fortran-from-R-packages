program test_clopper
   use genbinomapps
   implicit none
   type(confidence_interval) :: ci
   integer(i64) :: nreq
   integer :: fails

   fails = 0

   ci = clopper_pearson_ci(5_i64,100000_i64,0.05_dp,"upper")
   if (abs(ci%upper-0.0001051274511724034_dp) > 2.0e-12_dp) fails=fails+1

   ci = clopper_pearson_ci(5_i64,100000_i64,0.1_dp,"two.sided")
   if (abs(ci%lower-1.9701695640052194e-5_dp) > 2.0e-12_dp) fails=fails+1
   if (abs(ci%upper-0.0001051274511724034_dp) > 2.0e-12_dp) fails=fails+1

   ci = clopper_pearson_ci(0_i64,100_i64,0.1_dp,"two.sided")
   if (abs(ci%lower) > epsilon(1.0_dp)) fails=fails+1

   ci = clopper_pearson_ci(100_i64,100_i64,0.1_dp,"two.sided")
   if (abs(ci%upper-1.0_dp) > epsilon(1.0_dp)) fails=fails+1

   nreq = n_clopper_pearson(8_i64,0.0002_dp,0.1_dp)
   if (nreq /= 64972_i64) fails=fails+1

   if (fails /= 0) then
      print *, "test_clopper: FAIL", fails, ci%lower, ci%upper, nreq
      error stop 1
   end if
   print *, "test_clopper: PASS"
end program test_clopper
