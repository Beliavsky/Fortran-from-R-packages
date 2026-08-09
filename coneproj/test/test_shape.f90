program test_shape
   use coneproj
   implicit none
   real(dp) :: x(5), y(5)
   real(dp), allocatable :: delta(:,:)
   type(cone_result) :: ans
   integer :: status
   x = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   y = [3.0_dp,1.0_dp,2.0_dp,4.0_dp,5.0_dp]
   call make_delta(x,shape_increasing,delta,status)
   if (status /= coneproj_success) error stop 'make_delta status'
   call cone_b(y,transpose(delta),ans,vmat=reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp],[5,1]))
   if (ans%status /= coneproj_success) error stop 'shape cone status'
   if (any(ans%fit(2:5) < ans%fit(1:4)-1.0e-7_dp)) then
      print *, ans%fit
      error stop 'increasing shape violation'
   end if
   print *, 'test_shape: PASS'
end program test_shape
