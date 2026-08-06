! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

program test_gradients
  use rpeglmen, only : dp, exp_negative_log_likelihood, &
    grad_exp_negative_log_likelihood, gamma_negative_log_likelihood, &
    grad_gamma_negative_log_likelihood, prox_en, regularizer_en
  implicit none

  real(dp) :: a(6, 3), b(6), x(3), gradient(3), numerical(3)
  real(dp) :: xp(3), xm(3), h, shape, value
  integer :: j

  a(:, 1) = 1.0_dp
  a(:, 2) = [-1.0_dp, -0.6_dp, -0.2_dp, 0.2_dp, 0.6_dp, 1.0_dp]
  a(:, 3) = [0.5_dp, -0.3_dp, 0.7_dp, -0.8_dp, 0.2_dp, 0.9_dp]
  b = [0.8_dp, 1.1_dp, 0.6_dp, 1.7_dp, 1.2_dp, 2.0_dp]
  x = [0.15_dp, -0.2_dp, 0.1_dp]
  h = 1.0e-6_dp

  gradient = grad_exp_negative_log_likelihood(x, a, b)
  do j = 1, 3
    xp = x
    xm = x
    xp(j) = xp(j) + h
    xm(j) = xm(j) - h
    numerical(j) = (exp_negative_log_likelihood(xp, a, b) &
      - exp_negative_log_likelihood(xm, a, b)) / (2.0_dp * h)
  end do
  if (maxval(abs(gradient - numerical)) > 2.0e-6_dp) error stop 'exponential gradient mismatch'

  shape = 2.5_dp
  gradient = grad_gamma_negative_log_likelihood(x, a, b, shape)
  do j = 1, 3
    xp = x
    xm = x
    xp(j) = xp(j) + h
    xm(j) = xm(j) - h
    numerical(j) = (gamma_negative_log_likelihood(xp, a, b, shape) &
      - gamma_negative_log_likelihood(xm, a, b, shape)) / (2.0_dp * h)
  end do
  if (maxval(abs(gradient - numerical)) > 2.0e-6_dp) error stop 'Gamma gradient mismatch'

  value = regularizer_en([2.0_dp, -3.0_dp], 0.5_dp, .true.)
  if (abs(value - 5.75_dp) > 1.0e-12_dp) error stop 'elastic-net penalty mismatch'
  if (maxval(abs(prox_en([2.0_dp, -0.4_dp], 0.5_dp, 1.0_dp, 0.5_dp, .true., .false.) &
    - [1.4_dp, -0.12_dp])) > 1.0e-12_dp) error stop 'corrected proximal mismatch'

  print '(a)', 'test_gradients: PASS'
end program test_gradients
