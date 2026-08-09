program test_cone_b
   use coneproj
   implicit none
   real(dp) :: y(3), delta(3,2), v(3,1)
   type(cone_result) :: ans
   y = [3.0_dp,1.0_dp,2.0_dp]
   delta(:,1) = [-2.0_dp/3.0_dp, 1.0_dp/3.0_dp, 1.0_dp/3.0_dp]
   delta(:,2) = [-1.0_dp/3.0_dp,-1.0_dp/3.0_dp, 2.0_dp/3.0_dp]
   v(:,1) = 1.0_dp
   call cone_b(y,delta,ans,vmat=v)
   if (ans%status /= coneproj_success) error stop 'cone_b status'
   if (maxval(abs(ans%fit-[2.0_dp,2.0_dp,2.0_dp])) > 2.0e-5_dp) then
      print *, ans%fit
      error stop 'cone_b isotonic projection'
   end if
   print *, 'test_cone_b: PASS'
end program test_cone_b
