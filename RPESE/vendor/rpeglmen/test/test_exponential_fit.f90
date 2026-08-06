! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

program test_exponential_fit
  use rpeglmen, only : dp, enet_options, fit_result, glmnet_exp_fixed, &
    exp_negative_log_likelihood, predict_mean
  implicit none

  integer, parameter :: n = 500, p = 3
  real(dp) :: a(n, p), b(n), beta_true(p), u, mu(n), error
  type(enet_options) :: options
  type(fit_result) :: fit
  integer :: i
  integer(kind=8) :: state

  beta_true = [0.25_dp, 0.45_dp, -0.30_dp]
  do i = 1, n
    a(i, 1) = 1.0_dp
    a(i, 2) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
    a(i, 3) = sin(0.17_dp * real(i, dp))
  end do
  mu = predict_mean(a, beta_true)

  state = 24681357_8
  do i = 1, n
    state = modulo(16807_8 * state, 2147483647_8)
    u = (real(state, dp) + 0.5_dp) / 2147483647.0_dp
    b(i) = mu(i) * (-log(u))
  end do

  options%max_iter = 3000
  options%abs_tol = 1.0e-10_dp
  options%rel_tol = 1.0e-9_dp
  options%use_fista = .true.
  options%penalize_intercept = .false.
  call glmnet_exp_fixed(a, b, 0.0_dp, fit, options)

  if (.not. allocated(fit%coefficients)) error stop 'missing exponential coefficients'
  error = maxval(abs(fit%coefficients - beta_true))
  if (error > 0.18_dp) error stop 'exponential coefficient recovery failed'
  if (exp_negative_log_likelihood(fit%coefficients, a, b) &
    >= exp_negative_log_likelihood([0.0_dp, 0.0_dp, 0.0_dp], a, b)) then
    error stop 'exponential objective did not improve'
  end if

  print '(a)', 'test_exponential_fit: PASS'
end program test_exponential_fit
