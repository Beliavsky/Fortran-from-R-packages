program pso_weights
  use pwev
  implicit none
  real(dp) :: forecasts(30, 2), actual(30)
  type(pwev_control) :: control
  type(pwev_pso_result) :: result
  integer :: i
  do i = 1, 30
    forecasts(i, 1) = real(i, dp) / 30.0_dp
    forecasts(i, 2) = 1.0_dp + sin(0.3_dp * real(i, dp))
  end do
  actual = 0.25_dp * forecasts(:, 1) + 0.8_dp * forecasts(:, 2)
  control%pso_iterations = 250
  control%pso_population = 30
  call pso_ensemble_weights(actual, forecasts, control, result)
  if (result%status /= PWEV_SUCCESS) error stop 'PSO failed'
  print '(a,2f10.5)', 'weights: ', result%weights
  print '(a,es12.4)', 'SSE: ', result%objective
end program pso_weights
