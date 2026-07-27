! SPDX-License-Identifier: GPL-2.0-or-later
program test_estimation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchsk, only : dp, estimate_result, garchsk_est, gjrsk_est, &
      garchsk_initial_parameters, gjrsk_initial_parameters, garchsk_lik, gjrsk_lik, &
      garchsk_parameters_valid, gjrsk_parameters_valid
   use test_support, only : assert_true
   implicit none
   real(dp) :: data(90), initial_g(10), initial_j(13), initial_value
   type(estimate_result) :: fit
   integer :: i

   do i = 1, size(data)
      data(i) = 0.008_dp*sin(0.31_dp*real(i, dp)) + 0.004_dp*cos(0.11_dp*real(i, dp)) + &
         0.0015_dp*sin(1.7_dp*real(i, dp))
   end do

   initial_g = garchsk_initial_parameters(data)
   initial_value = garchsk_lik(initial_g, data)
   fit = garchsk_est(data, max_iterations=700, tolerance=1.0e-7_dp)
   call assert_true(all(ieee_is_finite(fit%params)), 'non-finite GARCHSK estimate')
   call assert_true(garchsk_parameters_valid(fit%params), 'infeasible GARCHSK estimate')
   call assert_true(fit%negative_log_likelihood <= initial_value + 1.0e-8_dp, 'GARCHSK objective did not improve')
   call assert_true(abs(fit%aic - (2.0_dp*fit%negative_log_likelihood + 20.0_dp)) < 1.0e-8_dp, &
      'corrected AIC formula not exposed')

   initial_j = gjrsk_initial_parameters(data)
   initial_value = gjrsk_lik(initial_j, data)
   fit = gjrsk_est(data, max_iterations=900, tolerance=1.0e-7_dp)
   call assert_true(all(ieee_is_finite(fit%params)), 'non-finite GJRSK estimate')
   call assert_true(gjrsk_parameters_valid(fit%params), 'infeasible GJRSK estimate')
   call assert_true(fit%negative_log_likelihood <= initial_value + 1.0e-8_dp, 'GJRSK objective did not improve')
   print '(a)', 'test_estimation: PASS'
end program test_estimation
