! SPDX-License-Identifier: Apache-2.0
program test_kalman
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use intraday_model
  use test_support, only : check, make_parameters
  implicit none
  type(volume_parameters) :: par
  type(kalman_output) :: kf
  real(dp), allocatable :: volume(:, :), states(:, :), y(:, :)
  real(dp) :: mse_pred, mse_filt, mse_smooth, symmetry

  par = make_parameters(8)
  call simulate_intraday_volume(par, 40, volume, states, seed=2026)
  allocate(y(size(volume, 1), size(volume, 2)))
  y = log(volume)
  call uniss_kalman(y, par, kf, smooth=.true.)

  call check(kf%status == intraday_ok, 'kalman status')
  call check(all(ieee_is_finite(kf%x_smooth)), 'finite smoothed states')
  mse_pred = sum((kf%x_pred - states)**2) / real(size(states), dp)
  mse_filt = sum((kf%x_filt - states)**2) / real(size(states), dp)
  mse_smooth = sum((kf%x_smooth - states)**2) / real(size(states), dp)
  call check(mse_filt < mse_pred, 'filter improves predictions')
  call check(mse_smooth <= 1.05_dp * mse_filt, 'smoother does not degrade filtered states')
  symmetry = maxval(abs(kf%v_smooth(1, 2, :) - kf%v_smooth(2, 1, :)))
  call check(symmetry < 1.0e-12_dp, 'smoothed covariance symmetry')
  call check(maxval(abs(kf%x_pred(:, 1) - par%x0)) < 1.0e-14_dp, 'initial prediction')

  write(*, '(a)') 'test_kalman: PASS'
end program test_kalman
