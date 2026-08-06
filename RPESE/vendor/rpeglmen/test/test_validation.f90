! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

program test_validation
  use rpeglmen, only : dp, enet_options, fit_result, glmnet_exp_fixed, &
    fit_glm_gamma_mle, rpe_invalid_input
  implicit none

  real(dp) :: a(4, 2), b(4)
  type(enet_options) :: options
  type(fit_result) :: fit

  a(:, 1) = 1.0_dp
  a(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
  b = [1.0_dp, 2.0_dp, -1.0_dp, 3.0_dp]

  call glmnet_exp_fixed(a, b, 0.1_dp, fit, options)
  if (fit%status /= rpe_invalid_input) error stop 'negative response was not rejected'

  b = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp]
  call fit_glm_gamma_mle(a, b, fit, options)
  if (fit%status /= rpe_invalid_input) error stop 'zero Gamma response was not rejected'

  options%alpha = 1.5_dp
  b = [1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]
  call glmnet_exp_fixed(a, b, 0.1_dp, fit, options)
  if (fit%status /= rpe_invalid_input) error stop 'invalid alpha was not rejected'

  print '(a)', 'test_validation: PASS'
end program test_validation
