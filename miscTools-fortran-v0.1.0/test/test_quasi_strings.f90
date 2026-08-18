program test_quasi_strings
   use misc_tools
   implicit none
   real(dp) :: m1(3,3),m2(3,3),z(3,3)
   character(len=8) :: alln(3),testn(2)
   integer :: fails,missing

   fails = 0
   m1 = reshape([0.0_dp,-1.0_dp,-1.0_dp, &
                 -1.0_dp,-2.0_dp,3.0_dp, &
                 -1.0_dp,3.0_dp,5.0_dp],[3,3])
   m2 = reshape([0.0_dp,1.0_dp,-1.0_dp, &
                 1.0_dp,-2.0_dp,3.0_dp, &
                 -1.0_dp,3.0_dp,5.0_dp],[3,3])
   z = 0.0_dp

   if (.not. quasiconcavity(m1)) fails=fails+1
   if (quasiconvexity(m1)) fails=fails+1
   if (quasiconcavity(m2)) fails=fails+1
   if (.not. quasiconvexity(m2)) fails=fails+1
   if (.not. quasiconcavity(z)) fails=fails+1
   if (.not. quasiconvexity(z)) fails=fails+1

   alln = ["alpha   ","beta    ","gamma   "]
   testn = ["alpha   ","gamma   "]
   if (.not. check_names(testn,alln,missing)) fails=fails+1
   testn = ["alpha   ","delta   "]
   if (check_names(testn,alln,missing)) fails=fails+1
   if (missing /= 2) fails=fails+1

   if (fails /= 0) then
      print *, "test_quasi_strings: FAIL", fails
      error stop 1
   end if
   print *, "test_quasi_strings: PASS"
end program test_quasi_strings
