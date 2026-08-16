program test_rng
   use polya_aeppli
   implicit none
   integer, parameter :: n = 150000
   integer, allocatable :: x(:)
   real(dp) :: lambda,prob,mean_emp,var_emp,mean_ref,var_ref
   integer :: fails

   fails = 0
   lambda = 4.0_dp
   prob = 0.25_dp
   allocate(x(n))
   call set_polya_aeppli_seed(12345)
   call r_polya_aeppli(x,lambda,prob)

   mean_emp = sum(real(x,dp))/real(n,dp)
   var_emp = sum((real(x,dp)-mean_emp)**2)/real(n-1,dp)
   mean_ref = polya_aeppli_mean(lambda,prob)
   var_ref = polya_aeppli_variance(lambda,prob)

   if (abs(mean_emp-mean_ref) > 0.035_dp) fails=fails+1
   if (abs(var_emp-var_ref) > 0.11_dp) fails=fails+1
   if (any(x < 0)) fails=fails+1

   if (fails /= 0) then
      print *, "test_rng: FAIL", fails, mean_emp, var_emp, mean_ref, var_ref
      error stop 1
   end if
   print *, "test_rng: PASS"
end program test_rng
