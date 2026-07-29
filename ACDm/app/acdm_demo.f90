! SPDX-License-Identifier: GPL-3.0-or-later
program acdm_demo
  use acdm
  implicit none
  integer, parameter :: n = 600
  type(acd_order) :: order
  type(acd_fit_options) :: options
  type(acd_fit_result) :: fit
  type(rng_state) :: rng
  real(dp) :: durations(n), true_parameters(3), forecast(5)
  integer :: status

  order = acd_order(p=1, r=0, q=1)
  true_parameters = [0.20_dp, 0.15_dp, 0.70_dp]
  call seed_rng(rng, 20250716)
  call simulate_acd(n, MODEL_ACD, order, true_parameters, &
                    DIST_EXPONENTIAL, [real(dp) ::], durations, status, rng, &
                    burn=300)
  if (status /= ACDM_SUCCESS) error stop 'simulation failed'

  options%model = MODEL_ACD
  options%dist = DIST_EXPONENTIAL
  options%order = order
  options%seed = 1948
  options%restarts = 2
  call acd_fit_model(durations, options, fit)
  if (fit%status /= ACDM_SUCCESS) error stop 'fit failed'

  call forecast_acd(fit, options, size(forecast), forecast, status)
  if (status /= ACDM_SUCCESS) error stop 'forecast failed'

  print '(a,3f12.6)', 'True parameters:      ', true_parameters
  print '(a,3f12.6)', 'Estimated parameters: ', fit%parameters
  print '(a,f14.5)', 'Log likelihood:       ', fit%loglik
  print '(a,f14.5)', 'AIC:                  ', fit%aic
  print '(a,5f12.6)', 'Conditional forecasts:', forecast
end program acdm_demo
