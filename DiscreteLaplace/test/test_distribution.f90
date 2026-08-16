program test_distribution
   use discrete_laplace
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   integer :: fails,k
   real(dp) :: p,q,s,c
   integer(int64) :: z
   fails=0; p=0.4_dp; q=0.3_dp
   s=0.0_dp
   do k=-100,100
      s=s+ddlaplace(k,p,q)
   end do
   if(abs(s-1.0_dp)>1.0e-12_dp) fails=fails+1
   if(abs(pdlaplace(0.0_dp,p,q)-sum([(ddlaplace(k,p,q),k=-100,0)]))>1.0e-12_dp) fails=fails+1
   do k=-6,6
      c=pdlaplace(real(k,dp),p,q)
      z=qdlaplace(min(c-1.0e-10_dp,1.0_dp-1.0e-12_dp),p,q)
      if(z/=k) fails=fails+1
   end do
   s=0.0_dp
   do k=-100,100
      s=s+ddlaplace2(k,p,q)
   end do
   if(abs(s-1.0_dp)>1.0e-12_dp) fails=fails+1
   do k=-6,6
      c=pdlaplace2(real(k,dp),p,q)
      z=qdlaplace2(min(c-1.0e-10_dp,1.0_dp-1.0e-12_dp),p,q)
      if(z/=k) fails=fails+1
   end do
   if(fails/=0) then
      print *, 'test_distribution: FAIL',fails
      error stop 1
   end if
   print *, 'test_distribution: PASS'
end program
