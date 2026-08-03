program test_pso
  use pwev
  implicit none
  real(dp) :: x(40, 3), actual(40)
  type(pwev_control) :: control
  type(pwev_pso_result) :: result
  integer :: i
  do i = 1, 40
    x(i, 1) = 0.5_dp + 0.03_dp * real(i, dp)
    x(i, 2) = 1.2_dp + sin(0.2_dp * real(i, dp))
    x(i, 3) = 0.8_dp + cos(0.13_dp * real(i, dp))
  end do
  actual = 0.2_dp * x(:, 1) + 0.7_dp * x(:, 2) + 0.1_dp * x(:, 3)
  control%pso_iterations = 400
  control%pso_population = 50
  control%random_seed = 17
  call pso_ensemble_weights(actual, x, control, result)
  if (result%status /= PWEV_SUCCESS) error stop 'PSO failed'
  if (any(result%weights < 0.0_dp) .or. any(result%weights > 1.0_dp)) error stop 'PSO bounds failed'
  if (result%objective > 1.0e-5_dp) error stop 'PSO objective too large'
  if (maxval(abs(result%weights - [0.2_dp, 0.7_dp, 0.1_dp])) > 2.0e-2_dp) &
    error stop 'PSO weights inaccurate'
  print *, 'test_pso: PASS'
end program test_pso
