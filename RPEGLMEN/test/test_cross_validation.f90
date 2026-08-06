! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

program test_cross_validation
  use rpeglmen, only : dp, enet_options, fit_result, path_result, &
    glmnet_exp, fit_regularization_path, model_exponential, predict_mean
  implicit none

  integer, parameter :: n = 180, p = 4
  real(dp) :: a(n, p), b(n), beta_true(p), mu(n), u
  type(enet_options) :: options
  type(fit_result) :: fit
  type(path_result) :: path
  integer :: i
  integer(kind=8) :: state

  beta_true = [0.2_dp, 0.35_dp, 0.0_dp, -0.25_dp]
  do i = 1, n
    a(i, 1) = 1.0_dp
    a(i, 2) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
    a(i, 3) = sin(0.07_dp * real(i, dp))
    a(i, 4) = cos(0.13_dp * real(i, dp))
  end do
  mu = predict_mean(a, beta_true)
  state = 13579_8
  do i = 1, n
    state = modulo(16807_8 * state, 2147483647_8)
    u = (real(state, dp) + 0.5_dp) / 2147483647.0_dp
    b(i) = mu(i) * (-log(u))
  end do

  options%num_lambda = 8
  options%k_fold = 4
  options%k_fold_iter = 2
  options%max_iter = 800
  options%abs_tol = 1.0e-8_dp
  options%rel_tol = 1.0e-7_dp
  options%cv_metric = 'nll'
  call glmnet_exp(a, b, fit, options)

  if (.not. allocated(fit%coefficients)) error stop 'CV coefficients missing'
  if (size(fit%lambda_grid) /= options%num_lambda) error stop 'wrong lambda-grid size'
  if (size(fit%cv_mean) /= options%num_lambda) error stop 'wrong CV result size'
  if (fit%selected_lambda <= 0.0_dp) error stop 'invalid selected lambda'
  if (any(fit%cv_sd < 0.0_dp)) error stop 'negative CV standard deviation'

  call fit_regularization_path(a, b, model_exponential, 1.0_dp, path, options, fit%lambda_grid)
  if (size(path%coefficients, 2) /= options%num_lambda) error stop 'path size mismatch'
  if (maxval(abs(path%coefficients(2:, 1))) > maxval(abs(path%coefficients(2:, options%num_lambda))) + 1.0e-6_dp) then
    error stop 'regularization path ordering failed'
  end if

  print '(a)', 'test_cross_validation: PASS'
end program test_cross_validation
