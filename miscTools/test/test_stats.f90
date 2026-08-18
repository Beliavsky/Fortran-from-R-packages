program test_stats
   use misc_tools
   implicit none
   real(dp) :: x,eps,num,r2
   real(dp), allocatable :: cm(:),rm(:),se(:)
   real(dp) :: m(6,4),a3(4,3,2),cov(2,2)
   type(coef_table_result) :: ct
   integer :: i,fails

   fails = 0

   eps = 1.0e-7_dp
   x = 1.3_dp
   num = (normal_density(x+eps)-normal_density(x-eps))/(2.0_dp*eps)
   if (abs(ddnorm(x)-num) > 2.0e-9_dp) fails=fails+1

   num = (normal_density2(4.0_dp+eps)-normal_density2(4.0_dp-eps))/(2.0_dp*eps)
   if (abs(ddnorm(4.0_dp,5.0_dp,2.0_dp)-num) > 2.0e-9_dp) fails=fails+1

   m = reshape([(real(i,dp),i=1,24)],shape(m))
   call col_medians(m,cm)
   if (maxval(abs(cm-[3.5_dp,9.5_dp,15.5_dp,21.5_dp])) > 1.0e-14_dp) fails=fails+1

   call row_medians(m,rm)
   if (maxval(abs(rm-[10.0_dp,11.0_dp,12.0_dp,13.0_dp,14.0_dp,15.0_dp])) > &
       1.0e-14_dp) fails=fails+1

   a3 = reshape([(real(i,dp),i=1,24)],shape(a3))
   block
      real(dp), allocatable :: c3(:,:)
      call col_medians_3d(a3,c3)
      if (abs(c3(2,2)-18.5_dp) > 1.0e-14_dp) fails=fails+1
   end block

   r2 = r_squared([1.0_dp,2.0_dp,3.0_dp],[0.1_dp,-0.2_dp,0.1_dp])
   if (abs(r2-0.97_dp) > 1.0e-14_dp) fails=fails+1

   cov = reshape([4.0_dp,1.0_dp,1.0_dp,9.0_dp],[2,2])
   call std_er(cov,se)
   if (maxval(abs(se-[2.0_dp,3.0_dp])) > 1.0e-14_dp) fails=fails+1

   ct = coef_table([1.0_dp,2.0_dp],[0.5_dp,0.25_dp],20.0_dp)
   if (abs(ct%table(1,3)-2.0_dp) > 1.0e-14_dp) fails=fails+1
   if (abs(ct%table(2,3)-8.0_dp) > 1.0e-14_dp) fails=fails+1
   if (abs(ct%table(1,4)-0.05926553544657044_dp) > 2.0e-13_dp) fails=fails+1
   if (abs(ct%table(2,4)-1.165662827148851e-7_dp) > 2.0e-16_dp) fails=fails+1

   if (n_obs_matrix(m) /= 6) fails=fails+1
   if (n_param_vector([1.0_dp,2.0_dp,3.0_dp]) /= 3) fails=fails+1

   if (fails /= 0) then
      print *, "test_stats: FAIL", fails
      error stop 1
   end if
   print *, "test_stats: PASS"

contains

   pure real(dp) function normal_density(z) result(v)
      real(dp), intent(in) :: z
      v = exp(-0.5_dp*z*z)/sqrt(2.0_dp*acos(-1.0_dp))
   end function normal_density

   pure real(dp) function normal_density2(z) result(v)
      real(dp), intent(in) :: z
      real(dp) :: q
      q = (z-5.0_dp)/2.0_dp
      v = exp(-0.5_dp*q*q)/(2.0_dp*sqrt(2.0_dp*acos(-1.0_dp)))
   end function normal_density2

end program test_stats
