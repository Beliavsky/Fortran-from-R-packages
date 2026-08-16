program test_rng
   use genbinomapps
   implicit none
   integer, parameter :: n = 150000
   integer(i64), allocatable :: x(:)
   integer(i64) :: sizes(3)
   real(dp) :: probs(3),mu,v,em,ev
   integer :: fails

   fails = 0
   sizes = [2_i64,5_i64,3_i64]
   probs = [0.8_dp,0.7_dp,0.3_dp]
   allocate(x(n))
   call set_genbinom_seed(12345)
   call rgbinom(x,sizes,probs)

   mu = generalized_binomial_mean(sizes,probs)
   v = generalized_binomial_variance(sizes,probs)
   em = sum(real(x,dp))/real(n,dp)
   ev = sum((real(x,dp)-em)**2)/real(n-1,dp)

   if (abs(em-mu) > 0.02_dp) fails=fails+1
   if (abs(ev-v) > 0.03_dp) fails=fails+1
   if (any(x < 0_i64) .or. any(x > sum(sizes))) fails=fails+1

   if (fails /= 0) then
      print *, "test_rng: FAIL", fails, em, ev, mu, v
      error stop 1
   end if
   print *, "test_rng: PASS"
end program test_rng
