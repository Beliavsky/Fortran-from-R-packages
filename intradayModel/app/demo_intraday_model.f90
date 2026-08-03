! SPDX-License-Identifier: Apache-2.0
program demo_intraday_model
  use intraday_model
  implicit none
  type(volume_parameters) :: truth
  type(volume_model) :: model
  type(volume_fit_control) :: control
  type(volume_decomposition) :: analysis, forecast
  real(dp), allocatable :: volume(:, :)
  integer :: i

  allocate(truth%phi(12))
  do i = 1, 12
    truth%phi(i) = -0.22_dp * cos(2.0_dp * acos(-1.0_dp) * real(i - 1, dp) / 12.0_dp)
  end do
  truth%phi = truth%phi - sum(truth%phi) / 12.0_dp
  truth%x0 = [3.2_dp, 0.0_dp]
  truth%a_eta = 0.86_dp
  truth%a_mu = 0.42_dp
  truth%var_eta = 0.014_dp
  truth%var_mu = 0.009_dp
  truth%r = 0.005_dp

  call simulate_intraday_volume(truth, 75, volume, seed=31415)
  control%acceleration = .true.
  control%maxit = 250
  control%abstol = 3.0e-5_dp
  control%save_history = .false.
  call fit_volume(volume(:, 1:60), model, control=control)
  call decompose_volume('analysis', model, volume(:, 1:60), analysis)
  call forecast_volume(model, volume, forecast, burn_in_days=60)

  write(*, '(a)') 'intradayModel Fortran demonstration'
  write(*, '(a,l1,a,i0)') 'fit converged: ', model%converged, ', iterations: ', model%iterations
  write(*, '(a,f12.6)') 'training smooth RMSE: ', analysis%error%rmse
  write(*, '(a,f12.6)') '15-day forecast RMSE: ', forecast%error%rmse
  write(*, '(a,3f11.6)') 'estimated variances: ', &
    model%par%var_eta, model%par%var_mu, model%par%r
end program demo_intraday_model
