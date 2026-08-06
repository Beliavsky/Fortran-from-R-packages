! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the RPEGLMEN computational code.

program test_gamma_fit
  use rpeglmen, only : dp, enet_options, fit_result, fit_glm_gamma_mle, &
    glm_gamma_net_fixed, fit_glm_gamma_net, predict_mean
  implicit none

  integer, parameter :: n = 700, p = 3, shape_integer = 3
  real(dp) :: a(n, p), b(n), beta_true(p), mu(n), u, gamma_noise, error
  type(enet_options) :: options
  type(fit_result) :: mle, penalized, selected
  integer :: i, j
  integer(kind=8) :: state

  beta_true = [0.10_dp, -0.35_dp, 0.28_dp]
  do i = 1, n
    a(i, 1) = 1.0_dp
    a(i, 2) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
    a(i, 3) = cos(0.11_dp * real(i, dp))
  end do
  mu = predict_mean(a, beta_true)

  state = 97531_8
  do i = 1, n
    gamma_noise = 0.0_dp
    do j = 1, shape_integer
      state = modulo(16807_8 * state, 2147483647_8)
      u = (real(state, dp) + 0.5_dp) / 2147483647.0_dp
      gamma_noise = gamma_noise - log(u)
    end do
    b(i) = mu(i) * gamma_noise / real(shape_integer, dp)
  end do

  options%max_iter = 2000
  options%abs_tol = 1.0e-9_dp
  options%rel_tol = 1.0e-8_dp
  options%penalize_intercept = .false.
  call fit_glm_gamma_mle(a, b, mle, options)

  if (.not. allocated(mle%coefficients)) error stop 'missing Gamma MLE coefficients'
  error = maxval(abs(mle%coefficients - beta_true))
  if (error > 0.18_dp) error stop 'Gamma coefficient recovery failed'
  if (mle%shape < 1.5_dp .or. mle%shape > 5.0_dp) error stop 'Gamma shape recovery failed'

  call glm_gamma_net_fixed(a, b, mle%shape, 0.01_dp, penalized, options, mle%coefficients)
  if (.not. allocated(penalized%coefficients)) error stop 'missing penalized Gamma coefficients'
  if (maxval(abs(penalized%coefficients)) > 10.0_dp) error stop 'unstable penalized Gamma fit'

  options%num_lambda = 6
  options%k_fold = 3
  options%k_fold_iter = 1
  options%max_iter = 800
  call fit_glm_gamma_net(a, b, selected, options)
  if (.not. allocated(selected%coefficients)) error stop 'Gamma CV coefficients missing'
  if (size(selected%lambda_grid) /= 6) error stop 'Gamma CV lambda grid missing'

  print '(a)', 'test_gamma_fit: PASS'
end program test_gamma_fit
