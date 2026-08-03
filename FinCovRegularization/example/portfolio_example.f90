! SPDX-License-Identifier: GPL-2.0-only
program portfolio_example
   use fincovregularization
   implicit none
   real(dp) :: covariance(3,3), minimum_variance(3), parity(3), objective
   integer :: status

   covariance = reshape([ &
      0.040_dp, 0.018_dp, 0.012_dp, &
      0.018_dp, 0.090_dp, 0.015_dp, &
      0.012_dp, 0.015_dp, 0.160_dp], [3,3])

   minimum_variance = gmvp(covariance, allow_short=.false., status=status)
   if (status /= fincov_ok) error stop fincov_status_message(status)
   parity = risk_parity(covariance, status=status, objective_value=objective)
   if (status /= fincov_ok) error stop fincov_status_message(status)

   print '(a,3f12.6)', 'Long-only GMVP: ', minimum_variance
   print '(a,3f12.6)', 'Risk parity:     ', parity
   print '(a,es14.5)', 'Risk-parity objective: ', objective
end program portfolio_example
