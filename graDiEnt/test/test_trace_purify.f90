program test_trace_purify
  use gradient
  use gradient_benchmarks, only : sphere_objective
  implicit none
  type(sqgde_options) :: opt
  type(sqgde_result) :: res
  integer :: expected_min_evals

  call get_algo_params(2,opt)
  opt%n_particles = 9
  opt%n_diff = 2
  opt%n_iter = 25
  opt%thin = 4
  opt%purify = 5
  opt%return_trace = .true.
  opt%stop_tol = 0.0_dp
  opt%init_center = 1.0_dp
  opt%init_sd = 1.0_dp
  call seed_rng(2222)
  call optim_sqgde(sphere_objective,opt,res)
  if (res%status /= 0) error stop 'trace: optimizer status'
  if (res%trace_count /= 7) error stop 'trace: wrong trace count'
  if (size(res%weights_trace,1) /= 7) error stop 'trace: wrong trace allocation'
  expected_min_evals = opt%n_particles + opt%n_iter*opt%n_particles + &
                       (opt%n_iter/opt%purify)*opt%n_particles
  if (res%evaluations < expected_min_evals) error stop 'purify: evaluation count too small'
  print *, 'PASS test_trace_purify', res%trace_count, res%evaluations
end program test_trace_purify
