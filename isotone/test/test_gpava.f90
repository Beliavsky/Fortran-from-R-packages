program test_gpava
   use isotone
   implicit none
   type(gpava_result) :: r, rt
   real(dp) :: y(9), w(9), z(6), yt(6), expected(9), ex2(6)
   integer :: i
   y=[-0.3_dp,1.2_dp,-0.5_dp,0.7_dp,0.6_dp,-1.0_dp,2.0_dp,1.5_dp,2.2_dp]
   w=1.0_dp
   expected=[-0.3_dp,0.2_dp,0.2_dp,0.2_dp,0.2_dp,0.2_dp,1.75_dp,1.75_dp,2.2_dp]
   call gpava_fit(y,r,weights=w)
   call assert_close(maxval(abs(r%x-expected)),1.0e-13_dp,'mean PAVA')
   do i=1,8
      if(r%x(i)>r%x(i+1)+1.0e-14_dp) error stop 'PAVA not monotone'
   end do
   call gpava_fit(y,r,weights=w,solver=GPAVA_MEDIAN)
   do i=1,8
      if(r%x(i)>r%x(i+1)+1.0e-14_dp) error stop 'median PAVA not monotone'
   end do
   call gpava_fit(y,r,weights=w,solver=GPAVA_FRACTILE,p=0.25_dp)
   do i=1,8
      if(r%x(i)>r%x(i+1)+1.0e-14_dp) error stop 'fractile PAVA not monotone'
   end do

   z=[1.0_dp,1.0_dp,2.0_dp,3.0_dp,3.0_dp,4.0_dp]
   yt=[3.0_dp,1.0_dp,2.0_dp,5.0_dp,1.0_dp,4.0_dp]
   call gpava_fit(yt,rt,z=z,ties=GPAVA_SECONDARY)
   ex2=[2.0_dp,2.0_dp,2.0_dp,3.0_dp,3.0_dp,4.0_dp]
   call assert_close(maxval(abs(rt%x-ex2)),1.0e-13_dp,'secondary ties')
   call gpava_fit(yt,rt,z=z,ties=GPAVA_TERTIARY)
   ex2=[3.0_dp,1.0_dp,2.0_dp,5.0_dp,1.0_dp,4.0_dp]
   ! Group fits are [2,2,3,4], so tertiary preserves within-tie deviations.
   call assert_close(maxval(abs(rt%x-ex2)),1.0e-13_dp,'tertiary ties')
   print *, 'test_gpava: PASS'
contains
   subroutine assert_close(err,tol,msg)
      real(dp),intent(in)::err,tol
      character(*),intent(in)::msg
      if(err>tol) then
         print *,trim(msg),err
         error stop 1
      end if
   end subroutine
end program
