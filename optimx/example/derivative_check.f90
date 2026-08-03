! SPDX-License-Identifier: GPL-2.0-only
program derivative_check_example
  use optimx_mod
  use optimx_example_functions, only: rosenbrock_callback
  implicit none
  type(optimx_problem) :: problem
  type(derivative_check) :: check
  call initialize_problem(problem,2)
  problem%objective => rosenbrock_callback
  problem%has_gradient=.true.;problem%has_hessian=.true.
  call grchk(problem,[0.8_dp,0.7_dp],check)
  write(*,'(a,es12.4)')'maximum relative gradient error: ',check%max_error
end program derivative_check_example
