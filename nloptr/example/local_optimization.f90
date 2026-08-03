! SPDX-License-Identifier: LGPL-3.0-or-later
program local_optimization
  use nloptr_mod
  use nloptr_example_functions, only: rosenbrock_objective
  implicit none
  type(nloptr_problem) :: problem
  type(nloptr_options) :: options
  type(nloptr_result) :: result

  problem%n = 2
  problem%objective => rosenbrock_objective
  allocate(problem%lower(2), problem%upper(2))
  problem%lower = -huge(1.0_dp)
  problem%upper = huge(1.0_dp)
  options = nl_opts('NLOPT_LD_LBFGS')
  options%maxeval = 5000

  call nloptr(problem, [-1.2_dp, 1.0_dp], options, result)
  write(*, '(a,2f14.8)') 'solution: ', result%solution
  write(*, '(a,es14.6)') 'objective: ', result%objective
  write(*, '(a,i0)') 'evaluations: ', result%evaluations
end program local_optimization
