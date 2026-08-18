program test_reference
   use frbinom
   implicit none
   real(dp), allocatable :: pmf(:)
   integer :: fails,status

   fails = 0

   ! The only upstream testthat regression value.
   if (abs(dfrbinom(4.0_dp,10,0.6_dp,0.6_dp,0.2_dp)- &
           0.07919005624468024_dp) > 2.0e-14_dp) fails=fails+1

   ! README distribution for size=50, p=.6, h=.7, c=.2.
   call frbinom_pmf_table(50,0.6_dp,0.7_dp,0.2_dp,.false.,pmf,status)
   if (status /= 0) fails=fails+1
   if (abs(pmf(0)-0.03309256052426757_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(1)-0.0014318634702621782_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(pmf(30)-0.0468892146069507_dp) > 5.0e-14_dp) fails=fails+1
   if (abs(pmf(50)-1.0704357695294731e-5_dp) > 2.0e-16_dp) fails=fails+1
   if (abs(sum(pmf)-1.0_dp) > 3.0e-14_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_reference: FAIL", fails
      error stop 1
   end if
   print *, "test_reference: PASS"
end program test_reference
