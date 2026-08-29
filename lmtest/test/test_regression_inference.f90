program test_regression_inference
   use lmtest, only : dp, lm_result, coefficient_test_result, &
      confidence_interval_result, lm_fit, coefficient_tests, &
      coefficient_confint, likelihood_ratio_test, test_result
   implicit none
   integer, parameter :: n = 80
   real(dp) :: x(n,3), y(n)
   type(lm_result) :: fit
   type(coefficient_test_result) :: ct
   type(confidence_interval_result) :: ci
   type(test_result) :: lr

   call make_data(x, y)
   fit = lm_fit(x, y)
   if (fit%info /= 0) error stop 1
   call check(fit%beta(1), 1.0093494567065597_dp, 2.0e-13_dp)
   call check(fit%beta(2), 2.0089355734605654_dp, 2.0e-13_dp)
   call check(fit%beta(3), -0.6930185737303571_dp, 2.0e-13_dp)
   call check(fit%rss, 3.084257180628489_dp, 2.0e-12_dp)
   call check(fit%loglik, 16.713548965499285_dp, 2.0e-12_dp)

   ct = coefficient_tests(fit%beta, fit%vcov, real(fit%df_resid,dp))
   call check(ct%std_error(1), 0.02239098200353376_dp, 2.0e-13_dp)
   call check(ct%statistic(3), -21.8716260581817_dp, 2.0e-11_dp)

   ci = coefficient_confint(fit%beta, fit%vcov, 0.95_dp, &
      real(fit%df_resid,dp))
   call check(ci%lower(2), 1.9698131073497256_dp, 2.0e-10_dp)
   call check(ci%upper(2), 2.0480580395714054_dp, 2.0e-10_dp)

   lr = likelihood_ratio_test(10.0_dp, 2, 14.0_dp, 4)
   call check(lr%statistic, 8.0_dp, 1.0e-14_dp)
contains
   subroutine make_data(a, b)
      real(dp), intent(out) :: a(:,:), b(:)
      integer :: j
      real(dp) :: v
      do j = 1, size(b)
         v = real(j,dp)
         a(j,1) = 1.0_dp
         a(j,2) = (v - 40.0_dp) / 20.0_dp
         a(j,3) = sin(0.3_dp * v)
         b(j) = 1.0_dp + 2.0_dp*a(j,2) - 0.7_dp*a(j,3) + &
            (0.18_dp + 0.002_dp*v)*sin(0.91_dp*v) + &
            0.08_dp*cos(0.17_dp*v)
      end do
   end subroutine make_data
   subroutine check(actual, expected, tol)
      real(dp), intent(in) :: actual, expected, tol
      if (abs(actual - expected) > tol) then
         print *, actual, expected
         error stop 1
      end if
   end subroutine check
end program test_regression_inference
