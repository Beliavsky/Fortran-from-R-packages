program test_rng
   use discrete_laplace
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   integer, parameter :: n=100000
   integer(int64) :: x(n)
   real(dp) :: m,v,tm
   type(edlaplace_result) :: e
   type(edlaplace2_result) :: e2
   integer :: fails
   fails=0
   call set_rng_seed(111)
   call rdlaplace(x,0.4_dp,0.3_dp)
   m=sum(real(x,dp))/real(n,dp)
   v=sum((real(x,dp)-m)**2)/real(n-1,dp)
   e=edlaplace(0.4_dp,0.3_dp)
   if(abs(m-e%e1)>0.03_dp .or. abs(v-e%v)>0.08_dp) fails=fails+1
   call set_rng_seed(222)
   call rdlaplace2(x,0.4_dp,0.3_dp)
   m=sum(real(x,dp))/real(n,dp)
   tm=sum(real(x,dp)**2)/real(n,dp)
   e2=edlaplace2(0.4_dp,0.3_dp)
   if(abs(m-e2%e1)>0.03_dp .or. abs(tm-e2%e2)>0.10_dp) fails=fails+1
   if(fails/=0) then
      print *, 'test_rng: FAIL',fails,m,tm,e2%e1,e2%e2
      error stop 1
   end if
   print *, 'test_rng: PASS'
end program
