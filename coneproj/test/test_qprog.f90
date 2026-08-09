program test_qprog
   use coneproj
   implicit none
   real(dp) :: q(2,2), c(2), a(2,2), b(2)
   type(qprog_result) :: ans
   q = 0.0_dp; q(1,1)=1.0_dp; q(2,2)=1.0_dp
   c = [-1.0_dp,2.0_dp]
   a = 0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp
   b = 0.0_dp
   call qprog(q,c,a,b,ans)
   if (ans%status /= coneproj_success) error stop 'qprog status'
   if (maxval(abs(ans%theta-[0.0_dp,2.0_dp])) > 1.0e-8_dp) then
      print *, ans%theta
      error stop 'qprog nonnegative'
   end if
   c = 0.0_dp
   b = [1.0_dp,0.0_dp]
   call qprog(q,c,a,b,ans)
   if (maxval(abs(ans%theta-[1.0_dp,0.0_dp])) > 1.0e-8_dp) then
      print *, ans%theta
      error stop 'qprog shifted constraints'
   end if
   print *, 'test_qprog: PASS'
end program test_qprog
