program test_properties
   use ccd, only : dp, i8, dcc, pcc, qcc
   implicit none
   integer :: fail, k
   real(dp) :: total, cprev, c
   integer(i8) :: qi
   fail = 0
   total = 0.0_dp
   do k = -200000, 200000
      total = total + dcc(real(k,dp), 0.0_dp, 0.7_dp)
   end do
   if (abs(total-1.0_dp) > 3.0e-6_dp) fail = fail + 1

   cprev = 0.0_dp
   do k = -50, 50
      c = pcc(int(k,i8), 0.0_dp, 0.7_dp)
      if (c < cprev) fail = fail + 1
      cprev = c
   end do

   do k = 1, 99
      qi = qcc(real(k,dp)/100.0_dp, 0.0_dp, 0.7_dp)
      if (pcc(qi, 0.0_dp, 0.7_dp) < real(k,dp)/100.0_dp) fail = fail + 1
      if (qi > -huge(0_i8)) then
         if (pcc(qi-1_i8, 0.0_dp, 0.7_dp) >= real(k,dp)/100.0_dp) fail = fail + 1
      end if
   end do

   if (fail == 0) then
      print *, 'test_properties: PASS'
   else
      print *, 'test_properties: FAIL', fail
      error stop 1
   end if
end program test_properties
