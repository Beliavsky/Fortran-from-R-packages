program test_best
  use gradient
  use gradient_benchmarks, only : shifted_sphere_objective
  implicit none
  type(sqgde_options) :: opt
  type(sqgde_result) :: res

  call get_algo_params(4,opt)
  opt%n_particles = 16
  opt%n_diff = 2
  opt%n_iter = 250
  opt%adapt_scheme = SQGDE_BEST
  opt%init_center = -3.0_dp
  opt%init_sd = 2.0_dp
  opt%stop_tol = 1.0e-11_dp
  call seed_rng(1237)
  call optim_sqgde(shifted_sphere_objective,opt,res)
  if (res%status /= 0) error stop 'best: optimizer status'
  if (res%weight > 1.0e-5_dp) error stop 'best: poor shifted-sphere solution'
  print *, 'PASS test_best', res%weight
end program test_best
