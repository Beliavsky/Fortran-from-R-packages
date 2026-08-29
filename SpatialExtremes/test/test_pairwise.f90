program test_pairwise
   use SpatialExtremes
   implicit none
   real(dp)::z(1,2),jac(1,2),a(1),rho(1),ll
   z(1,:)=[1.2_dp,2.3_dp]
   jac=0.0_dp
   a=0.8_dp
   rho=0.4_dp
   ll=lplik_smith(z,a,jac)
   call check(abs(ll-(-3.079153187092178_dp))<3e-11_dp,'Smith pair likelihood')
   ll=lplik_schlather(z,rho,jac)
   call check(abs(ll-(-3.203171220785689_dp))<3e-11_dp,'Schlather pair likelihood')
   ll=lplik_extremalt(z,rho,3.0_dp,jac)
   call check(abs(ll-(-3.252530499964498_dp))<2e-9_dp,'extremal-t pair likelihood')
   ll=lplik_schlather_ind(z,1.0_dp,rho,jac)
   call check(abs(ll-(-1/1.2_dp-1/2.3_dp-2*log(1.2_dp*2.3_dp)))<1e-13_dp,'independence mixture')
   print *,'test_pairwise: PASS'
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
