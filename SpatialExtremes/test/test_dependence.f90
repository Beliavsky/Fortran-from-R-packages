program test_dependence
   use SpatialExtremes
   implicit none
   real(dp)::x(4,3),m(3),v(3),fm(3),theta(3)
   x(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   x(:,2)=[1.0_dp,3.0_dp,2.0_dp,5.0_dp]
   x(:,3)=[4.0_dp,3.0_dp,2.0_dp,1.0_dp]
   m=madogram(x)
   v=variogram(x)
   call check(abs(m(1)-0.375_dp)<1e-14_dp,'madogram')
   call check(abs(v(1)-0.375_dp)<1e-14_dp,'variogram')
   fm=[0.0_dp,0.1_dp,1.0_dp/6.0_dp]
   theta=fmadogram_extcoeff(fm)
   call check(maxval(abs(theta-[1.0_dp,1.5_dp,2.0_dp]))<1e-14_dp,'F-madogram theta')
   call check(abs(extremal_coefficient_schlather(1.0_dp)-1.0_dp)<1e-14_dp,'complete dependence')
   call check(abs(extremal_coefficient_schlather_ind(0.0_dp,1.0_dp)-2.0_dp)<1e-14_dp,'independence')
   print *,'test_dependence: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      print *,'FAIL: ',msg
      error stop 1
      end if
   end subroutine
end program
