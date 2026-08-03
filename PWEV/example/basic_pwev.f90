program basic_pwev
  use pwev
  implicit none
  real(dp) :: data(72)
  type(pwev_control) :: control
  type(pwev_result) :: result
  integer :: i, status
  do i = 1, size(data)
    data(i) = 100.0_dp + 0.12_dp * real(i, dp) + 2.5_dp * sin(0.22_dp * real(i, dp))
  end do
  control%garch_max_iterations = 180
  control%mem_max_iterations = 120
  control%mem_random_starts = 3
  control%pso_iterations = 200
  control%pso_population = 30
  call pwev_fit(data, 0.85_dp, result, status, control)
  if (status /= PWEV_SUCCESS) error stop trim(result%message)
  print '(a,4f10.5)', 'weights: ', result%weights
  print '(a,f10.4)', 'test ensemble RMSE: ', result%accuracy(5, 10)
end program basic_pwev
