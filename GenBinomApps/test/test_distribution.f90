program test_distribution
   use genbinomapps
   implicit none
   integer(i64) :: sizes(3),x(3),q
   real(dp) :: probs(3),d(3),tail
   integer :: fails

   fails = 0
   sizes = [100_i64,100_i64,200_i64]
   probs = [0.001_dp,0.005_dp,0.01_dp]
   x = [0_i64,1_i64,2_i64]
   call dgbinom_vec(x,sizes,probs,d)

   if (abs(d(1)-0.07343376860548392_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(d(2)-0.19260316653501980_dp) > 3.0e-14_dp) fails=fails+1
   if (abs(d(3)-0.25173556277946479_dp) > 3.0e-14_dp) fails=fails+1

   tail = pgbinom(2.0_dp,sizes,probs,.false.)
   if (abs(tail-0.4822275020800305_dp) > 5.0e-14_dp) fails=fails+1

   q = qgbinom(1.0_dp-tail,sizes,probs)
   if (q /= 2_i64) fails=fails+1

   if (fails /= 0) then
      print *, "test_distribution: FAIL", fails, d, tail
      error stop 1
   end if
   print *, "test_distribution: PASS"
end program test_distribution
