program bootstrap_forecast_example
  use multiatsm, only : dp, bootstrap_result, forecast_result, BOOTSTRAP_IID, set_random_seed, &
    random_normal, bootstrap_var, fit_var, var_model, forecast_yields
  implicit none
  real(dp) :: factors(2, 300), a(3), b(3, 2)
  type(var_model) :: model
  type(bootstrap_result) :: bootstrap
  type(forecast_result) :: forecast
  integer :: t, info

  call set_random_seed(77)
  factors(:, 1) = 0.0_dp
  do t = 2, size(factors, 2)
    factors(:, t) = [0.01_dp, -0.005_dp] + &
      matmul(reshape([0.8_dp, 0.05_dp, 0.1_dp, 0.6_dp], [2, 2]), factors(:, t - 1)) + &
      [0.02_dp * random_normal(), 0.015_dp * random_normal()]
  end do
  call fit_var(factors, model, info)
  if (info /= 0) error stop 'fit_var failed'
  call bootstrap_var(factors, 50, BOOTSTRAP_IID, bootstrap, info, seed=991)
  if (info /= 0) error stop 'bootstrap_var failed'
  a = [0.01_dp, 0.015_dp, 0.02_dp]
  b = reshape([1.0_dp, 0.5_dp, 0.2_dp, 0.0_dp, 0.5_dp, 1.0_dp], [3, 2])
  call forecast_yields(model%intercept, model%phi, a, b, factors(:, size(factors, 2)), 6, forecast, info)
  if (info /= 0) error stop 'forecast_yields failed'
  write(*, '(a,3f11.6)') 'Six-step yield forecast: ', forecast%yields(:, 6)
  write(*, '(a,4f11.6)') 'Bootstrap median VAR feedback: ', bootstrap%median
end program bootstrap_forecast_example
