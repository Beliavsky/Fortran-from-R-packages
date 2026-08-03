program sigma_ensemble
  use pwev
  implicit none
  real(dp) :: proxy(70)
  type(pwev_control) :: control
  type(pwev_result) :: result
  integer :: i, status
  do i = 1, size(proxy)
    proxy(i) = 1.0_dp + 0.3_dp * abs(sin(0.19_dp * real(i, dp))) + 0.002_dp * real(i, dp)
  end do
  control%garch_output = PWEV_GARCH_SIGMA
  control%mem_oos_mode = PWEV_MEM_RECURSIVE_OOS
  control%garch_max_iterations = 180
  control%mem_max_iterations = 120
  control%mem_random_starts = 3
  control%pso_iterations = 200
  control%pso_population = 30
  call pwev_fit(proxy, 0.8_dp, result, status, control)
  if (status /= PWEV_SUCCESS) error stop trim(result%message)
  print '(a,4f10.5)', 'volatility weights: ', result%weights
  print '(a,3f10.5)', 'first ensemble forecasts: ', result%test_pred(1:3, 6)
end program sigma_ensemble
