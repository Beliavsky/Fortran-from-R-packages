! SPDX-License-Identifier: LGPL-3.0-or-later
program demo_nloptr
  use nloptr_mod
  use nloptr_example_functions, only: rosenbrock_objective
  implicit none
  type(nloptr_problem) :: problem
  type(nloptr_options) :: options
  type(nloptr_result) :: gradient_result, simplex_result

  problem%n = 2
  problem%objective => rosenbrock_objective
  allocate(problem%lower(2), problem%upper(2))
  problem%lower = -huge(1.0_dp)
  problem%upper = huge(1.0_dp)
  options = nl_opts()
  options%maxeval = 5000

  call lbfgs(problem, [-1.2_dp, 1.0_dp], options, gradient_result)
  call neldermead(problem, [-1.2_dp, 1.0_dp], options, simplex_result)

  write(*, '(a)') 'nloptr-fortran demonstration'
  write(*, '(a,2f12.7,a,es12.4)') 'BFGS:        ', gradient_result%solution, &
    ' f=', gradient_result%objective
  write(*, '(a,2f12.7,a,es12.4)') 'Nelder-Mead: ', simplex_result%solution, &
    ' f=', simplex_result%objective
end program demo_nloptr
