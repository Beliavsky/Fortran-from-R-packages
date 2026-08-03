! SPDX-License-Identifier: LGPL-3.0-or-later
program global_optimization
  use nloptr_mod
  use nloptr_example_functions, only: multimodal_objective
  implicit none
  type(nloptr_problem) :: problem
  type(nloptr_options) :: options
  type(nloptr_result) :: result

  problem%n = 2
  problem%objective => multimodal_objective
  allocate(problem%lower(2), problem%upper(2))
  problem%lower = -3.0_dp
  problem%upper = 3.0_dp
  options = nl_opts('NLOPT_GN_DIRECT_L')
  options%population = 80
  options%maxeval = 4000

  call nloptr(problem, [2.5_dp, -2.0_dp], options, result)
  write(*, '(a,2f14.8)') 'solution: ', result%solution
  write(*, '(a,es14.6)') 'objective: ', result%objective
end program global_optimization
