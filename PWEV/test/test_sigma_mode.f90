program test_sigma_mode
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pwev
  implicit none
  real(dp) :: train_data(50), test_data(10), train_models(50, 4), test_models(10, 4)
  type(pwev_control) :: control
  integer :: i, model_status(4)
  do i = 1, 50
    train_data(i) = 5.0_dp + 0.03_dp * real(i, dp) + 0.4_dp * sin(0.4_dp * real(i, dp))
  end do
  do i = 1, 10
    test_data(i) = 6.5_dp + 0.05_dp * real(i, dp)
  end do
  control%garch_output = PWEV_GARCH_SIGMA
  control%mem_oos_mode = PWEV_MEM_RECURSIVE_OOS
  control%garch_max_iterations = 150
  control%mem_max_iterations = 100
  control%mem_random_starts = 2
  call fit_pwev_base_models(train_data, test_data, control, train_models, test_models, model_status)
  if (any(model_status /= PWEV_SUCCESS)) error stop 'sigma-mode base fit failed'
  if (any(train_models(:, 1:3) <= 0.0_dp) .or. any(test_models(:, 1:3) <= 0.0_dp)) &
    error stop 'sigma forecast positivity failed'
  if (any(.not. ieee_is_finite(test_models))) error stop 'sigma forecast finite check failed'
  print *, 'test_sigma_mode: PASS'
end program test_sigma_mode
