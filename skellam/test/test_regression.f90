program test_regression
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use skellam, only : dp, i8, seed_random_number, skellam_regression_result, &
      fit_skellam_regression
   use skellam_special, only : random_poisson
   implicit none

   integer, parameter :: n = 500
   integer(i8) :: response(n), count1, count2
   real(dp) :: predictors(n,1), eta1, eta2
   type(skellam_regression_result) :: fit
   integer :: i, status

   call seed_random_number(271828)
   do i = 1, n
      predictors(i,1) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
      eta1 = 0.7_dp + 0.4_dp*predictors(i,1)
      eta2 = -0.2_dp + 0.4_dp*predictors(i,1)
      count1 = random_poisson(exp(eta1), status)
      count2 = random_poisson(exp(eta2), status)
      response(i) = count2 - count1
   end do

   call fit_skellam_regression(response, predictors, fit)
   if (.not. fit%converged) then
      print *, 'regression did not converge:', fit%status
      error stop 1
   end if
   if (maxval(abs(fit%beta1 - [0.7_dp, 0.4_dp])) > 0.22_dp) then
      print *, 'beta1 failed:', fit%beta1
      error stop 1
   end if
   if (maxval(abs(fit%beta2 - [-0.2_dp, 0.4_dp])) > 0.28_dp) then
      print *, 'beta2 failed:', fit%beta2
      error stop 1
   end if
   if (.not. all(ieee_is_finite(fit%standard_error1))) error stop 'regression standard errors failed'
   if (any(fit%standard_error1 <= 0.0_dp) .or. any(fit%standard_error2 <= 0.0_dp)) then
      error stop 'regression covariance failed'
   end if
   print '(a)', 'test_regression: PASS'
end program test_regression
