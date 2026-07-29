program active_risk_parity
   use risk_parity_portfolio_mod
   implicit none
   integer, parameter :: n = 4
   real(dp) :: sigma(n, n), budgets(n), expected_returns(n), weights(n)
   integer :: info, iterations

   sigma = reshape([ &
      0.0400_dp, 0.0060_dp, 0.0040_dp, 0.0020_dp, &
      0.0060_dp, 0.0625_dp, 0.0075_dp, 0.0030_dp, &
      0.0040_dp, 0.0075_dp, 0.0900_dp, 0.0090_dp, &
      0.0020_dp, 0.0030_dp, 0.0090_dp, 0.1225_dp], [n, n])
   budgets = 1.0_dp / real(n, dp)
   expected_returns = [0.050_dp, 0.060_dp, 0.075_dp, 0.085_dp]

   call active_risk_parity_ccd(sigma, budgets, expected_returns, 20.0_dp, &
                               0.02_dp, weights, info, iterations, 1.0e-8_dp, 5000)
   if (info /= RPP_OK .and. info /= RPP_MAX_ITER) error stop 'active solver failed'

   write(*, '(a,4f12.7)') 'weights: ', weights
   write(*, '(a,4f12.7)') 'relative risk contributions: ', &
                          relative_risk_contributions(sigma, weights)
   write(*, '(a,i0)') 'iterations: ', iterations
end program active_risk_parity
