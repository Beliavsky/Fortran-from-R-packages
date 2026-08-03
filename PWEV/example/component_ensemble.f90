program component_ensemble
  use pwev
  implicit none
  real(dp) :: y_train(24), y_test(6), f_train(24, 4), f_test(6, 4)
  type(pwev_control) :: control
  type(pwev_result) :: result
  integer :: i, status
  do i = 1, 24
    f_train(i, :) = [1.0_dp + 0.04_dp * i, 1.5_dp + sin(0.2_dp * i), &
      1.2_dp + cos(0.17_dp * i), 0.8_dp + 0.03_dp * i]
  end do
  y_train = matmul(f_train, [0.1_dp, 0.6_dp, 0.2_dp, 0.1_dp])
  do i = 1, 6
    f_test(i, :) = [2.0_dp + 0.04_dp * i, 1.5_dp + sin(0.2_dp * (i + 24)), &
      1.2_dp + cos(0.17_dp * (i + 24)), 1.5_dp + 0.03_dp * i]
  end do
  y_test = matmul(f_test, [0.1_dp, 0.6_dp, 0.2_dp, 0.1_dp])
  control%pso_iterations = 300
  control%pso_population = 40
  call pwev_fit_from_components(y_train, y_test, f_train, f_test, result, status, control)
  if (status /= PWEV_SUCCESS) error stop trim(result%message)
  print '(a,4f10.5)', 'estimated weights: ', result%weights
end program component_ensemble
