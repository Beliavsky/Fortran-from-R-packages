program test_native_reference
   use lsei
   implicit none
   real(dp) :: a(5,3), b(5), bb(5,2)
   type(ls_result) :: r
   type(hfti_result) :: h
   a=reshape([1d0,2d0,0d0,1d0,3d0, 0d0,1d0,2d0,1d0,0d0, &
              1d0,0d0,1d0,2d0,1d0],[5,3])
   b=[1d0,2d0,-1d0,.5d0,3d0]
   call nnls_solve(a,b,r)
   if (r%mode/=LSEI_SUCCESS) error stop 1
   if (maxval(abs(r%x-[.96666666666666645_dp,0.0_dp,0.0_dp]))>2e-14_dp) error stop 2
   if (abs(r%rnorm-1.1105554165971787_dp)>2e-14_dp) error stop 3
   call pnnls_solve(a,b,1,r)
   if (maxval(abs(r%x-[.96666666666666645_dp,0.0_dp,0.0_dp]))>2e-14_dp) error stop 4
   bb(:,1)=b; bb(:,2)=[2d0,-1d0,1d0,0d0,4d0]
   call hfti_solve(a,bb,h,1d-10)
   if (h%krank/=3) error stop 5
   if (maxval(abs(h%rnorm-[.26196841599779169_dp,3.0524821144150107_dp]))>3e-14_dp) error stop 6
   if (maxval(abs(h%x(:,1)-[1.1078431372549018_dp,-.35294117647058809_dp, &
       -.17647058823529427_dp]))>3e-14_dp) error stop 7
   if (maxval(abs(h%x(:,2)-[.57647058823529396_dp,-.74117647058823521_dp, &
       .92941176470588238_dp]))>3e-14_dp) error stop 8
   print *, 'PASS test_native_reference'
end program
