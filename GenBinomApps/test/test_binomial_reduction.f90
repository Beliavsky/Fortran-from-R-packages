program test_binomial_reduction
   use genbinomapps
   implicit none
   integer(i64) :: sizes(1),k
   real(dp) :: probs(1),p,ref,s
   integer :: fails,n

   fails = 0
   n = 12
   sizes = [int(n,i64)]
   probs = [0.37_dp]
   s = 0.0_dp
   do k = 0_i64, int(n,i64)
      ref = exp(log_gamma(real(n+1,dp))-log_gamma(real(k+1_i64,dp)) - &
            log_gamma(real(n-int(k)+1,dp)) + real(k,dp)*log(probs(1)) + &
            real(n-int(k),dp)*log(1.0_dp-probs(1)))
      p = dgbinom(k,sizes,probs)
      if (abs(p-ref) > 2.0e-14_dp) fails=fails+1
      s = s+p
   end do
   if (abs(s-1.0_dp) > 2.0e-14_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_binomial_reduction: FAIL", fails
      error stop 1
   end if
   print *, "test_binomial_reduction: PASS"
end program test_binomial_reduction
