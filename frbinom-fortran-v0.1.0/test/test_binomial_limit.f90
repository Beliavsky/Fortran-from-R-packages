program test_binomial_limit
   use frbinom
   implicit none
   integer, parameter :: n=20
   real(dp), allocatable :: pmf(:)
   real(dp) :: p,ref
   integer :: k,fails,status

   fails = 0
   p = 0.37_dp
   call frbinom_pmf_table(n,p,0.73_dp,0.0_dp,.false.,pmf,status)
   if (status /= 0) fails=fails+1

   do k = 0, n
      ref = exp(log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp)) - &
            log_gamma(real(n-k+1,dp)) + real(k,dp)*log(p) + &
            real(n-k,dp)*log(1.0_dp-p))
      if (abs(pmf(k)-ref) > 3.0e-13_dp) fails=fails+1
   end do

   if (abs(sum(pmf)-1.0_dp) > 5.0e-14_dp) fails=fails+1

   ! size=1 is Bernoulli(prob), including start=.true.
   call frbinom_pmf_table(1,p,0.73_dp,0.2_dp,.true.,pmf,status)
   if (abs(pmf(0)-(1.0_dp-p)) > 1.0e-15_dp) fails=fails+1
   if (abs(pmf(1)-p) > 1.0e-15_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_binomial_limit: FAIL", fails
      error stop 1
   end if
   print *, "test_binomial_limit: PASS"
end program test_binomial_limit
