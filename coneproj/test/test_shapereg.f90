program test_shapereg
   use coneproj
   implicit none
   integer, parameter :: n=8
   real(dp) :: x(n), y(n)
   type(shapereg_result) :: ans
   integer :: i
   do i=1,n
      x(i)=real(i,dp)
      y(i)=0.5_dp+0.25_dp*x(i)
   end do
   call shapereg_fit(y,x,shape_increasing,ans)
   if (ans%status /= coneproj_success) error stop 'shapereg status'
   if (maxval(abs(ans%constrained_fit-y)) > 2.0e-5_dp) then
      print *, ans%constrained_fit
      error stop 'shapereg exact increasing fit'
   end if
   print *, 'test_shapereg: PASS'
end program test_shapereg
