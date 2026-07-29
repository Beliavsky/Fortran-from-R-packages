! SPDX-License-Identifier: MIT
program fit_and_predict
   use vasicekfit, only : dp, vasicek_fit_result, prediction_result, fit_vasicek, predict_quantiles
   implicit none
   real(dp) :: y(8), macro(8,1), stressed(1,1)
   type(vasicek_fit_result) :: fit
   type(prediction_result) :: forecast

   y = [0.012_dp, 0.018_dp, 0.021_dp, 0.028_dp, 0.035_dp, 0.047_dp, 0.061_dp, 0.072_dp]
   macro(:,1) = [-1.4_dp, -1.0_dp, -0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp, 1.0_dp, 1.4_dp]
   stressed(1,1) = 2.0_dp

   fit = fit_vasicek(y, macro)
   if (.not. fit%ok) error stop fit%message
   forecast = predict_quantiles(fit, [0.50_dp, 0.99_dp], stressed)
   if (.not. forecast%ok) error stop forecast%message

   print '(a,2f12.6)', 'stressed median and 99% loss: ', forecast%values(1,:)
end program fit_and_predict
