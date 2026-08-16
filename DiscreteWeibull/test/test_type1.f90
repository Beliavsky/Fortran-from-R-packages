program test_type1
   use discrete_weibull
   implicit none
   real(dp) :: q,beta,s,p
   integer(i64) :: x,k
   integer :: fails

   fails = 0
   q = 0.6_dp
   beta = 0.8_dp

   s = 0.0_dp
   do k = 1_i64, 10000_i64
      s = s+ddweibull(k,q,beta,.false.)
   end do
   if (abs(s-1.0_dp) > 2.0e-11_dp) fails=fails+1

   s = 0.0_dp
   do k = 0_i64, 10000_i64
      s = s+ddweibull(k,q,beta,.true.)
   end do
   if (abs(s-1.0_dp) > 2.0e-11_dp) fails=fails+1

   do k = 1_i64, 30_i64
      p = pdweibull(real(k,dp),q,beta,.false.)
      x = qdweibull(p-1.0e-12_dp,q,beta,.false.)
      if (x > k) fails=fails+1
      if (pdweibull(real(k,dp),q,beta,.false.) < &
          pdweibull(real(k-1_i64,dp),q,beta,.false.)) fails=fails+1
   end do

   if (qdweibull(0.0_dp,q,beta,.false.) /= 1_i64) fails=fails+1
   if (qdweibull(0.0_dp,q,beta,.true.) /= 0_i64) fails=fails+1

   if (fails /= 0) then
      print *, "test_type1: FAIL", fails
      error stop 1
   end if
   print *, "test_type1: PASS"
end program test_type1
