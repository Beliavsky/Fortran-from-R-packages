program test_cone_a
   use coneproj
   implicit none
   real(dp) :: y(3), a(2,3)
   type(cone_result) :: ans
   y = [3.0_dp,1.0_dp,2.0_dp]
   a = reshape([-1.0_dp,0.0_dp, 1.0_dp,-1.0_dp, 0.0_dp,1.0_dp],[2,3])
   call cone_a(y,a,ans)
   if (ans%status /= coneproj_success) error stop 'cone_a status'
   if (maxval(abs(ans%fit-[2.0_dp,2.0_dp,2.0_dp])) > 1.0e-8_dp) then
      print *, ans%fit
      error stop 'cone_a isotonic projection'
   end if
   print *, 'test_cone_a: PASS'
end program test_cone_a
