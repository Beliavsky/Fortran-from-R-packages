program test_base_models
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pwev
  implicit none
  real(dp) :: train_data(48), test_data(12), train_models(48, 4), test_models(12, 4)
  type(pwev_control) :: control
  integer :: i, model_status(4)
  do i = 1, 48
    train_data(i) = 100.0_dp + 0.15_dp * real(i, dp) + 2.0_dp * sin(0.27_dp * real(i, dp))
  end do
  do i = 1, 12
    test_data(i) = 107.0_dp + 0.15_dp * real(i, dp) + 1.5_dp * sin(0.31_dp * real(i, dp))
  end do
  control%garch_max_iterations = 150
  control%mem_max_iterations = 100
  control%mem_random_starts = 2
  call fit_pwev_base_models(train_data, test_data, control, train_models, test_models, model_status)
  if (any(model_status /= PWEV_SUCCESS)) error stop 'base-model fit failed'
  if (any(.not. ieee_is_finite(train_models)) .or. any(.not. ieee_is_finite(test_models))) &
    error stop 'non-finite base prediction'
  if (any(train_models(:, 4) <= 0.0_dp) .or. any(test_models(:, 4) <= 0.0_dp)) &
    error stop 'MEM positivity failed'
  print *, 'test_base_models: PASS'
end program test_base_models
