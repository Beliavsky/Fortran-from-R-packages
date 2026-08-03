! SPDX-License-Identifier: MIT
program nelson_siegel_example
   use fixedincome
   implicit none
   type(spot_rate_curve_t) :: curve
   type(interpolation_t) :: initial
   type(fit_result_t) :: fit
   real(dp) :: years(12), rates(12)
   integer :: i, status

   years = [(0.5_dp*real(i,dp), i=1,12)]
   rates = nelson_siegel(years, 0.06_dp, -0.025_dp, 0.035_dp, 0.8_dp)
   curve = spotratecurve(rates, term(years, 'years'), 'continuous', &
                         'actual/365', 'actual', status=status)
   initial = interp_nelsonsiegel(0.05_dp, -0.01_dp, 0.01_dp, 0.4_dp)
   fit = fit_interpolation(initial, curve, max_iterations=3000)

   print '(a,es12.4)', 'sum of squared errors: ', fit%objective
   print '(a,4(f10.6,1x))', 'parameters: ', parameters(fit%model)
end program nelson_siegel_example
