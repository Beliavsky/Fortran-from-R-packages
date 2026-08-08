program test_current
  use gradient
  use gradient_benchmarks, only : sphere_objective
  implicit none
  type(sqgde_options) :: opt
  type(sqgde_result) :: res

  call get_algo_params(3,opt)
  opt%n_particles = 15
  opt%n_diff = 2
  opt%n_iter = 350
  opt%adapt_scheme = SQGDE_CURRENT
  opt%init_center = 4.0_dp
  opt%init_sd = 2.0_dp
  opt%stop_tol = 1.0e-12_dp
  call seed_rng(1236)
  call optim_sqgde(sphere_objective,opt,res)
  if (res%status /= 0) error stop 'current: optimizer status'
  if (res%weight > 1.0e-4_dp) error stop 'current: poor sphere solution'
  print *, 'PASS test_current', res%weight
end program test_current
