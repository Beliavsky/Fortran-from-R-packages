! SPDX-License-Identifier: Apache-2.0
program forecast_example
  use intraday_model
  implicit none
  type(volume_parameters) :: truth
  type(volume_model) :: model
  type(volume_fit_control) :: control
  type(volume_decomposition) :: forecast
  real(dp), allocatable :: volume(:, :)
  integer :: i

  allocate(truth%phi(8))
  do i = 1, 8
    truth%phi(i) = -0.15_dp * cos(2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / 8.0_dp)
  end do
  truth%phi = truth%phi - sum(truth%phi) / 8.0_dp
  truth%x0 = [3.0_dp, 0.0_dp]
  truth%a_eta = 0.90_dp
  truth%a_mu = 0.50_dp
  truth%var_eta = 0.010_dp
  truth%var_mu = 0.010_dp
  truth%r = 0.004_dp

  call simulate_intraday_volume(truth, 65, volume, seed=2030)
  control%acceleration = .true.
  control%maxit = 250
  control%abstol = 3.0e-5_dp
  control%save_history = .false.
  call fit_volume(volume(:, 1:50), model, control=control)
  call forecast_volume(model, volume, forecast, burn_in_days=50)

  write(*, '(a,i0)') 'out-of-sample bins: ', size(forecast%fitted_signal)
  write(*, '(a,f12.6)') 'MAE:  ', forecast%error%mae
  write(*, '(a,f12.6)') 'MAPE: ', forecast%error%mape
  write(*, '(a,f12.6)') 'RMSE: ', forecast%error%rmse
end program forecast_example
