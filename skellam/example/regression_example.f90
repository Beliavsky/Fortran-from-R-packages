program regression_example
   use skellam, only : dp, i8, seed_random_number, skellam_regression_result, &
      fit_skellam_regression
   use skellam_special, only : random_poisson
   implicit none
   integer, parameter :: n = 300
   real(dp) :: x(n,1), eta1, eta2
   integer(i8) :: y(n), count1, count2
   integer :: i, status
   type(skellam_regression_result) :: fit

   call seed_random_number(77)
   do i = 1, n
      x(i,1) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
      eta1 = 1.0_dp + 0.5_dp*x(i,1)
      eta2 = -0.5_dp + 0.5_dp*x(i,1)
      count1 = random_poisson(exp(eta1), status)
      count2 = random_poisson(exp(eta2), status)
      y(i) = count2 - count1
   end do
   call fit_skellam_regression(y, x, fit)
   print '(a,l1)', 'converged: ', fit%converged
   print '(a,*(f11.5,1x))', 'beta1: ', fit%beta1
   print '(a,*(f11.5,1x))', 'beta2: ', fit%beta2
   print '(a,*(f11.5,1x))', 'se1:   ', fit%standard_error1
   print '(a,*(f11.5,1x))', 'se2:   ', fit%standard_error2
end program regression_example
