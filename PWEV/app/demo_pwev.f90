program demo_pwev
  use pwev
  implicit none
  real(dp) :: data(80)
  type(pwev_control) :: control
  type(pwev_result) :: result
  integer :: i, status

  do i = 1, size(data)
    data(i) = 90.0_dp + 0.2_dp * real(i, dp) + 2.0_dp * sin(0.18_dp * real(i, dp)) + &
      0.5_dp * cos(0.47_dp * real(i, dp))
  end do
  control%garch_max_iterations = 200
  control%mem_max_iterations = 140
  control%mem_random_starts = 3
  control%pso_iterations = 250
  control%pso_population = 36
  control%random_seed = 2026
  call pwev_fit(data, 0.85_dp, result, status, control)
  if (status /= PWEV_SUCCESS) error stop trim(result%message)

  print '(a)', 'PWEV modern Fortran demo'
  print '(a,i0,a,i0)', 'training observations: ', result%train_size, ', test observations: ', result%test_size
  print '(a,4f10.5)', 'PSO weights: ', result%weights
  print '(a,f12.5)', 'training SSE: ', result%pso_objective
  print '(a,f10.4)', 'ensemble train RMSE: ', result%accuracy(5, 1)
  print '(a,f10.4)', 'ensemble test RMSE: ', result%accuracy(5, 10)
  print '(a,3f11.4)', 'first test actual/base ensemble: ', result%test_pred(1, 1), &
    result%test_pred(1, 2), result%test_pred(1, 6)
end program demo_pwev
