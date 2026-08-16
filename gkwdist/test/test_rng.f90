program test_rng
   use gkwdist
   implicit none
   real(dp),allocatable :: x(:)
   real(dp) :: m,theory
   integer :: fails
   fails=0
   call seed_rng(12345)
   x=rbeta_(50000,2.0_dp,2.0_dp)
   m=sum(x)/real(size(x),dp)
   theory=2.0_dp/(2.0_dp+3.0_dp)
   if(abs(m-theory)>0.008_dp) then
      print '(a,2f12.6)','beta RNG mean mismatch: ',m,theory; fails=fails+1
   end if
   if(any(x<0.0_dp) .or. any(x>1.0_dp)) then
      print '(a)','beta RNG support failure'; fails=fails+1
   end if
   x=rgkw(50000,1.7_dp,2.4_dp,1.3_dp,0.8_dp,1.2_dp)
   if(any(x<0.0_dp) .or. any(x>1.0_dp)) then
      print '(a)','gkw RNG support failure'; fails=fails+1
   end if
   if(fails==0) then
      print '(a)','test_rng: PASS'
   else
      print '(a,i0)','test_rng: FAIL ',fails
      error stop 1
   end if
end program test_rng
