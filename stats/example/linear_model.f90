program linear_model
   use r_stats, only: dp, lm_fit, lm_fit_t, predict_lm
   implicit none

   real(dp) :: x(5, 1), y(5), new_x(2, 1)
   type(lm_fit_t) :: fit

   x(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   y = 1.5_dp + 2.0_dp*x(:, 1)
   fit = lm_fit(y, x)

   print "(a, *(f10.4, 1x))", "coefficients: ", fit%coef
   new_x(:, 1) = [5.0_dp, 6.0_dp]
   print "(a, *(f10.4, 1x))", "predictions:  ", predict_lm(fit, new_x)
end program linear_model
