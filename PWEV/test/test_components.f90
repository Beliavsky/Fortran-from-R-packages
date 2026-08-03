program test_components
  use pwev
  implicit none
  real(dp) :: train_actual(30), test_actual(8), train_models(30, 4), test_models(8, 4)
  type(pwev_control) :: control
  type(pwev_result) :: result
  integer :: i, status
  do i = 1, 30
    train_models(i, 1) = 10.0_dp + 0.1_dp * real(i, dp)
    train_models(i, 2) = 9.0_dp + sin(0.2_dp * real(i, dp))
    train_models(i, 3) = 11.0_dp + cos(0.3_dp * real(i, dp))
    train_models(i, 4) = 8.0_dp + 0.05_dp * real(i, dp)
  end do
  train_actual = 0.15_dp * train_models(:, 1) + 0.55_dp * train_models(:, 2) + &
    0.2_dp * train_models(:, 3) + 0.1_dp * train_models(:, 4)
  do i = 1, 8
    test_models(i, 1) = 13.0_dp + 0.1_dp * real(i, dp)
    test_models(i, 2) = 9.0_dp + sin(0.2_dp * real(i + 30, dp))
    test_models(i, 3) = 11.0_dp + cos(0.3_dp * real(i + 30, dp))
    test_models(i, 4) = 9.5_dp + 0.05_dp * real(i, dp)
  end do
  test_actual = 0.15_dp * test_models(:, 1) + 0.55_dp * test_models(:, 2) + &
    0.2_dp * test_models(:, 3) + 0.1_dp * test_models(:, 4)
  control%pso_iterations = 350
  control%pso_population = 40
  control%random_seed = 91
  control%round_accuracy = .false.
  call pwev_fit_from_components(train_actual, test_actual, train_models, test_models, result, status, control)
  if (status /= PWEV_SUCCESS) error stop 'component workflow failed'
  if (any(shape(result%train_fitted) /= [30, 6])) error stop 'train shape failed'
  if (any(shape(result%test_pred) /= [8, 6])) error stop 'test shape failed'
  if (any(shape(result%accuracy) /= [5, 18])) error stop 'accuracy shape failed'
  if (maxval(abs(result%train_fitted(:, 6) - train_actual)) > 2.0e-2_dp) &
    error stop 'ensemble train fit failed'
  print *, 'test_components: PASS'
end program test_components
