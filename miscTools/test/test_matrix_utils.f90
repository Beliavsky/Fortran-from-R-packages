program test_matrix_utils
   use misc_tools
   implicit none
   real(dp), allocatable :: m(:,:),v(:),m2(:,:),out(:,:)
   real(dp) :: base(3,3)
   integer :: fails,status,i

   fails = 0

   call sym_matrix([(real(i,dp),i=1,10)],m,4)
   if (maxval(abs(m-reshape([ &
      1.0_dp,2.0_dp,3.0_dp,4.0_dp, &
      2.0_dp,5.0_dp,6.0_dp,7.0_dp, &
      3.0_dp,6.0_dp,8.0_dp,9.0_dp, &
      4.0_dp,7.0_dp,9.0_dp,10.0_dp],[4,4]))) > 1.0e-14_dp) fails=fails+1

   base = reshape([11.0_dp,12.0_dp,13.0_dp, &
                   12.0_dp,22.0_dp,23.0_dp, &
                   13.0_dp,23.0_dp,33.0_dp],[3,3])
   call vecli(base,v)
   if (maxval(abs(v-[11.0_dp,12.0_dp,13.0_dp,22.0_dp,23.0_dp,33.0_dp])) > &
       1.0e-14_dp) fails=fails+1

   call vecli2m(v,m2,status)
   if (status /= 0 .or. maxval(abs(m2-base)) > 1.0e-14_dp) fails=fails+1
   if (veclipos(1,2,3) /= 2) fails=fails+1

   call triang([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],3,m2)
   if (maxval(abs(m2-reshape([ &
      1.0_dp,0.0_dp,0.0_dp, &
      2.0_dp,4.0_dp,0.0_dp, &
      3.0_dp,5.0_dp,0.0_dp],[3,3]))) > 1.0e-14_dp) fails=fails+1

   base = reshape([(real(i,dp),i=1,9)],[3,3])
   call insert_row(base,2,[10.0_dp,11.0_dp,12.0_dp],out)
   if (size(out,1) /= 4 .or. size(out,2) /= 3) fails=fails+1
   if (maxval(abs(out(2,:)-[10.0_dp,11.0_dp,12.0_dp])) > 1.0e-14_dp) fails=fails+1
   if (maxval(abs(out(3,:)-base(2,:))) > 1.0e-14_dp) fails=fails+1

   call insert_col(base,2,[10.0_dp,11.0_dp,12.0_dp],out)
   if (size(out,1) /= 3 .or. size(out,2) /= 4) fails=fails+1
   if (maxval(abs(out(:,2)-[10.0_dp,11.0_dp,12.0_dp])) > 1.0e-14_dp) fails=fails+1
   if (maxval(abs(out(:,3)-base(:,2))) > 1.0e-14_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_matrix_utils: FAIL", fails
      error stop 1
   end if
   print *, "test_matrix_utils: PASS"

end program test_matrix_utils
