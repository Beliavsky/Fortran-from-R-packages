program strategies_example
  use gradient
  use gradient_benchmarks, only : shifted_sphere_objective
  implicit none
  type(sqgde_options) :: opt
  type(sqgde_result) :: res
  integer :: scheme
  character(len=7), parameter :: names(3) = [character(len=7) :: 'rand', 'current', 'best']

  do scheme = SQGDE_RAND, SQGDE_BEST
    call get_algo_params(3,opt)
    opt%n_particles = 15
    opt%n_diff = 2
    opt%n_iter = 250
    opt%adapt_scheme = scheme
    opt%init_center = -2.0_dp
    opt%init_sd = 2.0_dp
    call seed_rng(5000+scheme)
    call optim_sqgde(shifted_sphere_objective,opt,res)
    write(*,'(a7,2x,es14.6)') names(scheme), res%weight
  end do
end program strategies_example
