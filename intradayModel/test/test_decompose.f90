! SPDX-License-Identifier: Apache-2.0
program test_decompose
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use intraday_model
  use test_support, only : check, make_parameters
  implicit none
  type(volume_parameters) :: truth
  type(volume_model) :: model
  type(volume_decomposition) :: analysis, forecast, wrapper
  real(dp), allocatable :: volume(:, :), states(:, :), dirty(:, :), clean(:, :)
  integer, allocatable :: kept(:)
  integer :: status
  character(len=100) :: message

  truth = make_parameters(8)
  call simulate_intraday_volume(truth, 30, volume, states, seed=91)
  model%par = truth
  model%converged = .true.
  model%status = intraday_ok

  call decompose_volume('analysis', model, volume, analysis)
  call check(analysis%status == intraday_ok, 'analysis status')
  call check(size(analysis%fitted_signal) == size(volume), 'analysis length')
  call check(maxval(abs(analysis%fitted_signal - analysis%daily * analysis%dynamic * &
                        analysis%seasonal)) < 1.0e-12_dp, 'component product')
  call check(maxval(abs(analysis%residual * analysis%fitted_signal - &
                        analysis%original_signal)) < 1.0e-10_dp, 'residual identity')
  call check(analysis%error%mae >= 0.0_dp .and. ieee_is_finite(analysis%error%rmse), &
             'analysis errors')

  call decompose_volume('forecast', model, volume, forecast, burn_in_days=5)
  call forecast_volume(model, volume, wrapper, burn_in_days=5)
  call check(forecast%status == intraday_ok, 'forecast status')
  call check(size(forecast%fitted_signal) == 8 * 25, 'forecast burn-in length')
  call check(maxval(abs(forecast%fitted_signal - wrapper%fitted_signal)) < 1.0e-14_dp, &
             'forecast wrapper')
  call check(forecast%is_forecast, 'forecast result flag')

  allocate(dirty(size(volume, 1), size(volume, 2)))
  dirty = volume
  dirty(1, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
  call clean_volume_data(dirty, clean, kept, status, message)
  call check(status == intraday_ok, 'clean data status')
  call check(size(clean, 2) == size(volume, 2) - 1, 'invalid day removed')
  call check(all(kept /= 2), 'kept-column indices')

  write(*, '(a)') 'test_decompose: PASS'
end program test_decompose
