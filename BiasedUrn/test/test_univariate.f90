program test_univariate
   use biasedurn
   implicit none
   real(dp), parameter :: tol = 5.0e-11_dp
   real(dp) :: p, s
   integer :: x, failures

   failures = 0

   call check_close(dfnchypergeo(4, 10, 15, 8, 2.5_dp), &
      0.30973837935631643_dp, tol, 'Fisher pmf', failures)
   call check_close(pfnchypergeo(4, 10, 15, 8, 2.5_dp), &
      0.509059198458439_dp, tol, 'Fisher cdf', failures)
   call check_close(meanfnchypergeo(10, 15, 8, 2.5_dp), &
      4.458638750003264_dp, tol, 'Fisher mean', failures)
   call check_close(varfnchypergeo(10, 15, 8, 2.5_dp), &
      1.3533504970824788_dp, tol, 'Fisher variance', failures)
   if (qfnchypergeo(0.5_dp, 10, 15, 8, 2.5_dp) /= 4) failures = failures + 1

   call check_close(dwnchypergeo(4, 10, 15, 8, 2.5_dp), &
      0.2799878275368307_dp, 2.0e-10_dp, 'Wallenius pmf', failures)
   call check_close(pwnchypergeo(4, 10, 15, 8, 2.5_dp), &
      0.42668352985243596_dp, 2.0e-10_dp, 'Wallenius cdf', failures)
   call check_close(meanwnchypergeo(10, 15, 8, 2.5_dp), &
      4.690532085568021_dp, 2.0e-10_dp, 'Wallenius mean', failures)
   call check_close(varwnchypergeo(10, 15, 8, 2.5_dp), &
      1.336614996492543_dp, 2.0e-10_dp, 'Wallenius variance', failures)
   if (qwnchypergeo(0.5_dp, 10, 15, 8, 2.5_dp) /= 5) failures = failures + 1

   s = 0.0_dp
   do x = 0, 8
      s = s + dwnchypergeo(x, 10, 15, 8, 2.5_dp)
   end do
   call check_close(s, 1.0_dp, 2.0e-10_dp, 'Wallenius normalization', failures)

   ! Equal odds reduce to the ordinary hypergeometric distribution.
   p = exp(log_gamma(11.0_dp) - log_gamma(5.0_dp) - log_gamma(7.0_dp) &
      + log_gamma(16.0_dp) - log_gamma(5.0_dp) - log_gamma(12.0_dp) &
      - log_gamma(26.0_dp) + log_gamma(9.0_dp) + log_gamma(18.0_dp))
   call check_close(dfnchypergeo(4, 10, 15, 8, 1.0_dp), p, tol, &
      'Fisher central case', failures)
   call check_close(dwnchypergeo(4, 10, 15, 8, 1.0_dp), p, tol, &
      'Wallenius central case', failures)

   if (failures == 0) then
      print *, 'test_univariate: PASS'
   else
      print *, 'test_univariate: FAIL', failures
      error stop 1
   end if

contains

   subroutine check_close(got, expected, atol, label, failures)
      real(dp), intent(in) :: got, expected, atol
      character(*), intent(in) :: label
      integer, intent(inout) :: failures
      if (abs(got - expected) > atol) then
         print *, trim(label), ' got=', got, ' expected=', expected
         failures = failures + 1
      end if
   end subroutine check_close

end program test_univariate
