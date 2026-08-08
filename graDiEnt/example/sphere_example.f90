program sphere_example
  use gradient
  use gradient_benchmarks, only : sphere_objective
  implicit none
  type(sqgde_options) :: opt
  type(sqgde_result) :: res

  call get_algo_params(5,opt)
  opt%n_particles = 20
  opt%n_diff = 3
  opt%n_iter = 300
  opt%adapt_scheme = SQGDE_RAND
  opt%init_center = 3.0_dp
  opt%init_sd = 2.0_dp
  call seed_rng(2026)
  call optim_sqgde(sphere_objective,opt,res)
  write(*,'(a,es14.6)') 'best objective = ', res%weight
  write(*,'(a,*(f12.6,1x))') 'solution = ', res%solution
end program sphere_example
