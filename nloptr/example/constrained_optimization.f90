! SPDX-License-Identifier: LGPL-3.0-or-later
program constrained_optimization
  use nloptr_mod
  use nloptr_example_functions, only: tutorial_objective, tutorial_constraints
  implicit none
  type(nloptr_problem) :: problem
  type(nloptr_options) :: options
  type(nloptr_result) :: result

  problem%n = 2
  problem%n_ineq = 2
  problem%objective => tutorial_objective
  problem%inequality => tutorial_constraints
  allocate(problem%lower(2), problem%upper(2))
  problem%lower = [-huge(1.0_dp), 0.0_dp]
  problem%upper = huge(1.0_dp)
  options = nl_opts('NLOPT_LD_MMA')
  options%maxeval = 12000
  options%max_outer = 10

  call nloptr(problem, [1.234_dp, 5.678_dp], options, result)
  write(*, '(a,2f14.8)') 'solution: ', result%solution
  write(*, '(a,es14.6)') 'objective: ', result%objective
  write(*, '(a,es14.6)') 'max violation: ', result%max_constraint
end program constrained_optimization
