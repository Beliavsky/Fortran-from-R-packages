program test_full_workflow
  use pwev
  implicit none
  real(dp) :: data(64)
  type(pwev_control) :: control
  type(pwev_result) :: result
  integer :: i, status
  do i = 1, 64
    data(i) = 80.0_dp + 0.25_dp * real(i, dp) + 3.0_dp * sin(0.21_dp * real(i, dp))
  end do
  control%garch_max_iterations = 150
  control%mem_max_iterations = 100
  control%mem_random_starts = 2
  control%pso_iterations = 150
  control%pso_population = 24
  control%random_seed = 31415
  call pwev_fit(data, 0.8_dp, result, status, control)
  if (status /= PWEV_SUCCESS) error stop 'full PWEV workflow failed'
  if (result%train_size /= 51 .or. result%test_size /= 13) error stop 'split size failed'
  if (size(result%weights) /= 4) error stop 'weight size failed'
  if (any(result%weights < 0.0_dp) .or. any(result%weights > 1.0_dp)) error stop 'weight bounds failed'
  if (any(shape(result%accuracy) /= [5, 18])) error stop 'accuracy table failed'
  print *, 'test_full_workflow: PASS'
end program test_full_workflow
