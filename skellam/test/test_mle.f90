program test_mle
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use skellam, only : dp, i8, rskellam, seed_random_number, skellam_mle_result, fit_skellam_mle
   implicit none

   integer(i8), allocatable :: sample(:)
   type(skellam_mle_result) :: fit

   call seed_random_number(314159)
   sample = rskellam(3000, 10.0_dp, 6.0_dp)
   call fit_skellam_mle(sample, fit)
   if (.not. fit%converged) then
      print *, 'MLE did not converge:', fit%status
      error stop 1
   end if
   if (abs(fit%lambda1 - 10.0_dp) > 0.55_dp .or. abs(fit%lambda2 - 6.0_dp) > 0.55_dp) then
      print *, 'MLE estimates failed:', fit%lambda1, fit%lambda2
      error stop 1
   end if
   if (.not. ieee_is_finite(fit%standard_error1 + fit%standard_error2)) error stop 'MLE standard errors failed'
   if (fit%standard_error1 <= 0.0_dp .or. fit%standard_error2 <= 0.0_dp) error stop 'MLE covariance failed'
   print '(a)', 'test_mle: PASS'
end program test_mle
