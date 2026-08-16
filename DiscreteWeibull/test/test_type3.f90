program test_type3
   use discrete_weibull
   implicit none
   real(dp) :: c,beta,s,p,h
   integer(i64) :: k,x
   integer :: fails

   fails = 0
   c = 0.3_dp
   beta = 0.75_dp

   s = 0.0_dp
   do k = 0_i64, 200_i64
      s = s+ddweibull3(k,c,beta)
   end do
   if (abs(s-1.0_dp) > 1.0e-12_dp) fails=fails+1

   do k = 0_i64, 20_i64
      p = pdweibull3(real(k,dp),c,beta)
      x = qdweibull3(max(0.0_dp,p-1.0e-12_dp),c,beta)
      if (x > k) fails=fails+1
   end do

   h = hdweibull3(3_i64,c,beta)
   if (abs(h-(1.0_dp-exp(-c*4.0_dp**beta))) > 1.0e-14_dp) fails=fails+1

   ! beta=0 reduces to a geometric distribution on {0,1,...}.
   beta = 0.0_dp
   if (abs(ddweibull3(0_i64,c,beta)-(1.0_dp-exp(-c))) > 1.0e-14_dp) &
      fails=fails+1
   if (abs(pdweibull3(2.0_dp,c,beta)-(1.0_dp-exp(-3.0_dp*c))) > 1.0e-14_dp) &
      fails=fails+1

   ! beta=-1 quantile uses harmonic bracketing and must invert the CDF.
   beta = -1.0_dp
   x = qdweibull3(0.8_dp,2.0_dp,beta)
   if (pdweibull3(real(x,dp),2.0_dp,beta) < 0.8_dp) fails=fails+1
   if (x > 0_i64) then
      if (pdweibull3(real(x-1_i64,dp),2.0_dp,beta) >= 0.8_dp) fails=fails+1
   end if

   if (fails /= 0) then
      print *, "test_type3: FAIL", fails
      error stop 1
   end if
   print *, "test_type3: PASS"
end program test_type3
