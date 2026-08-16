program test_beta
   use genbinomapps
   implicit none
   real(dp) :: p,x
   integer :: fails

   fails = 0
   p = beta_cdf(0.3_dp,2.5_dp,4.0_dp)
   if (abs(p-0.3521975859067671_dp) > 2.0e-13_dp) fails=fails+1

   x = beta_quantile(0.95_dp,6.0_dp,99995.0_dp)
   if (abs(x-0.0001051274511724034_dp) > 2.0e-12_dp) fails=fails+1

   if (fails /= 0) then
      print *, "test_beta: FAIL", fails, p, x
      error stop 1
   end if
   print *, "test_beta: PASS"
end program test_beta
