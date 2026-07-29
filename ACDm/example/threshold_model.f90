! SPDX-License-Identifier: GPL-3.0-or-later
program threshold_model
  use acdm
  implicit none
  integer, parameter :: n = 300
  type(acd_order) :: order
  type(rng_state) :: rng
  real(dp), parameter :: breakpoints(2) = [0.8_dp, 1.2_dp]
  real(dp), allocatable :: parameters(:)
  real(dp) :: durations(n), mu(n), residuals(n), loglik
  integer :: status

  order = acd_order(p=1, r=0, q=1)
  parameters = default_model_parameters(MODEL_TACD, order, 1.0_dp, &
                                        nbreak=size(breakpoints))
  call seed_rng(rng, 31051998)
  call simulate_acd(n, MODEL_TACD, order, parameters, DIST_WEIBULL, &
                    [1.4_dp], durations, status, rng, burn=250, &
                    breakpoints=breakpoints)
  if (status /= ACDM_SUCCESS) error stop 'TACD simulation failed'
  loglik = acd_loglik(durations, MODEL_TACD, order, parameters, &
                      DIST_WEIBULL, [1.4_dp], .true., mu, residuals, status, &
                      breakpoints=breakpoints)
  if (status /= ACDM_SUCCESS) error stop 'TACD likelihood failed'
  print '(a,f14.5)', 'TACD log likelihood: ', loglik
  print '(a,3f10.5)', 'First conditional means: ', mu(1:3)
end program threshold_model
