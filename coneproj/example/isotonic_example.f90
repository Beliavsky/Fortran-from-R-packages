program isotonic_example
   use coneproj
   implicit none
   real(dp) :: y(5), a(4,5)
   type(cone_result) :: ans
   integer :: i
   y=[3.0_dp,1.0_dp,2.0_dp,4.0_dp,3.5_dp]
   a=0.0_dp
   do i=1,4
      a(i,i)=-1.0_dp
      a(i,i+1)=1.0_dp
   end do
   call cone_a(y,a,ans)
   print '(a,*(f10.5,1x))','projection: ',ans%fit
   print '(a,i0)','active-set steps: ',ans%steps
end program isotonic_example
