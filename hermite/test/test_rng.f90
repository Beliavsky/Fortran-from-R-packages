program test_rng
   use hermite
   implicit none
   integer, parameter :: n=150000
   integer(i64), allocatable :: x(:)
   real(dp) :: em,ev,mu,v
   integer :: fails

   fails=0
   allocate(x(n))
   call set_hermite_seed(12345)
   call rhermite(x,1.7_dp,0.8_dp,3)
   em=sum(real(x,dp))/real(n,dp)
   ev=sum((real(x,dp)-em)**2)/real(n-1,dp)
   mu=hermite_mean(1.7_dp,0.8_dp,3)
   v=hermite_variance(1.7_dp,0.8_dp,3)

   if (abs(em-mu) > 0.025_dp) fails=fails+1
   if (abs(ev-v) > 0.10_dp) fails=fails+1
   if (any(x<0_i64)) fails=fails+1

   if (fails/=0) then
      print *,'test_rng: FAIL',fails,em,ev,mu,v
      error stop 1
   end if
   print *,'test_rng: PASS'
end program test_rng
