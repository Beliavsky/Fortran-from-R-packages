! SPDX-License-Identifier: MIT
program vasicekfit_demo
   use vasicekfit, only : dp, vasicek_fit_result, confidence_interval_result, &
      fit_vasicek, vasicek_confidence_intervals
   implicit none
   real(dp) :: y(8), macro(8,1)
   type(vasicek_fit_result) :: fit
   type(confidence_interval_result) :: intervals

   y = [0.012_dp, 0.018_dp, 0.021_dp, 0.028_dp, 0.035_dp, 0.047_dp, 0.061_dp, 0.072_dp]
   macro(:,1) = [-1.4_dp, -1.0_dp, -0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp, 1.0_dp, 1.4_dp]
   fit = fit_vasicek(y, macro)
   if (.not. fit%ok) error stop fit%message
   intervals = vasicek_confidence_intervals(fit)

   print '(a,f10.6)', 'p     = ', fit%p
   print '(a,f10.6)', 'rho   = ', fit%rho
   print '(a,f10.6)', 'kappa = ', fit%kappa(1)
   print '(a,2f12.6)', '95% CI for kappa: ', intervals%lower(3), intervals%upper(3)
end program vasicekfit_demo
