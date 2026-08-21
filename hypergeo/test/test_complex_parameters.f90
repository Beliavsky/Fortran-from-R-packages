! SPDX-License-Identifier: GPL-2.0-only
program test_complex_parameters
 use hypergeo_fortran, only: dp, hypergeo
 implicit none
 complex(dp) :: a(5),b(5),c(5),z(5),ref(5),v
 integer :: i
 a=[cmplx(1.125_dp,-.105_dp,dp),cmplx(-.045_dp,1.005_dp,dp),cmplx(-.195_dp,.675_dp,dp), &
    cmplx(.015_dp,-.585_dp,dp),cmplx(-.795_dp,.015_dp,dp)]
 b=[cmplx(-.705_dp,-.825_dp,dp),cmplx(.225_dp,.165_dp,dp),cmplx(.765_dp,1.245_dp,dp), &
    cmplx(-.375_dp,-1.335_dp,dp),cmplx(1.125_dp,-.855_dp,dp)]
 c=[cmplx(-.825_dp,-.675_dp,dp),cmplx(-.345_dp,.915_dp,dp),cmplx(.825_dp,.735_dp,dp), &
    cmplx(-.765_dp,-1.125_dp,dp),cmplx(-1.455_dp,.105_dp,dp)]
 z=[cmplx(.135_dp,.285_dp,dp),cmplx(1.305_dp,-1.305_dp,dp),cmplx(-.465_dp,.165_dp,dp), &
    cmplx(1.215_dp,-1.035_dp,dp),cmplx(-.195_dp,.915_dp,dp)]
 ref=[cmplx(.964064283226315_dp,.406653008068512_dp,dp), &
      cmplx(.82707630709444_dp,-.370937724712521_dp,dp), &
      cmplx(1.03665840072975_dp,-.3745212608222_dp,dp), &
      cmplx(.311643703830213_dp,-.059274615332188_dp,dp), &
      cmplx(1.31889143955355_dp,.98358042791712_dp,dp)]
 do i=1,5
   v=hypergeo(a(i),b(i),c(i),z(i))
   if (abs(v-ref(i)) > 1.0e-12_dp) then
      write(*,'(a,i0,1x,es14.6)') 'FAIL complex SAGE case ',i,abs(v-ref(i))
      error stop 1
   end if
 end do
print *, 'test_complex_parameters: PASS'
end program test_complex_parameters
