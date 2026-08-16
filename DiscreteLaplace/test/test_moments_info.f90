program test_moments_info
   use discrete_laplace
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   type(edlaplace_result) :: a
   type(edlaplace2_result) :: b
   real(dp) :: m(2,2),mi(2,2),eye(2,2),e1,e2
   integer :: k,fails,st
   real(dp), parameter :: p=0.4_dp,q=0.3_dp
   fails=0
   a=edlaplace(p,q)
   e1=0.0_dp; e2=0.0_dp
   do k=-200,200
      e1=e1+real(k,dp)*ddlaplace(k,p,q)
      e2=e2+real(k*k,dp)*ddlaplace(k,p,q)
   end do
   if(abs(a%e1-e1)>1.0e-11_dp) fails=fails+1
   if(abs(a%v-(e2-e1*e1))>1.0e-10_dp) fails=fails+1
   b=edlaplace2(p,q)
   e1=0.0_dp; e2=0.0_dp
   do k=-250,250
      e1=e1+real(k,dp)*ddlaplace2(k,p,q)
      e2=e2+real(k*k,dp)*ddlaplace2(k,p,q)
   end do
   if(abs(b%e1-e1)>1.0e-10_dp) fails=fails+1
   if(abs(b%e2-e2)>1.0e-9_dp) fails=fails+1
   m=ifi2(p,q,st)
   if(st/=0 .or. (.not. all(ieee_is_finite(m)))) fails=fails+1
   ! iFI is already the inverse Fisher information in the upstream package.
   mi=ifi(p,q)
   if((.not. all(ieee_is_finite(mi))) .or. mi(1,1)<=0.0_dp .or. mi(2,2)<=0.0_dp) fails=fails+1
   eye=0.0_dp; eye(1,1)=1.0_dp; eye(2,2)=1.0_dp
   if(fails/=0) then
      print *, 'test_moments_info: FAIL',fails
      error stop 1
   end if
   print *, 'test_moments_info: PASS'
end program
