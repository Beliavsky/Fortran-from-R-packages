! SPDX-License-Identifier: GPL-2.0-or-later
program durbin_h_example
   use jfe
   implicit none

   real(dp), parameter :: residuals(6) = [0.8_dp, 0.4_dp, -0.1_dp, -0.3_dp, 0.2_dp, 0.1_dp]
   type(durbin_h_result) :: result

   result = durbin_h(residuals, n_fitted=20, n_coefficients=4, lag_variance=0.01_dp)
   print '(a,f12.6)', 'Durbin-Watson: ', result%durbin_watson
   print '(a,f12.6)', 'Durbin h:      ', result%statistic
   print '(a,f12.6)', 'p-value:       ', result%p_value
end program durbin_h_example
