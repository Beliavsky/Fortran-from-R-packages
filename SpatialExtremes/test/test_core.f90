program test_core
   use SpatialExtremes
   implicit none
   real(dp)::p,x,c(3,3),coord(3,2),d(3)
   p=pgev(0.7_dp,0.2_dp,1.3_dp,0.15_dp)
   x=qgev(p,0.2_dp,1.3_dp,0.15_dp)
   call check(abs(x-0.7_dp)<1e-11_dp,'GEV p/q')
   p=pgpd(1.5_dp,0.1_dp,0.8_dp,0.2_dp)
   x=qgpd(p,0.1_dp,0.8_dp,0.2_dp)
   call check(abs(x-1.5_dp)<1e-11_dp,'GPD p/q')
   coord(1,:)=[0.0_dp,0.0_dp]
   coord(2,:)=[1.0_dp,0.0_dp]
   coord(3,:)=[0.0_dp,1.0_dp]
   d=euclidean_distances(coord)
   call check(maxval(abs(d-[1.0_dp,1.0_dp,sqrt(2.0_dp)]))<1e-14_dp,'distances')
   c=covariance_matrix(coord,COV_POWEREXP,0.1_dp,0.9_dp,2.0_dp,1.5_dp)
   call check(abs(c(1,1)-1.0_dp)<1e-14_dp,'cov diagonal')
   call check(abs(c(1,2)-0.9_dp*exp(-(0.5_dp)**1.5_dp))<1e-14_dp,'power exponential')
   call check(abs(extremal_coefficient_smith(0.8_dp)-1.3108434832_dp)<1e-9_dp,'Smith theta')
   print *,'test_core: PASS'
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
