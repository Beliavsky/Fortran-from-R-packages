! SPDX-License-Identifier: GPL-2.0-only
program demo_optimx
  use optimx_mod
  use optimx_example_functions, only: rosenbrock_callback
  implicit none
  type(optimx_problem) :: problem
  type(optimx_result) :: answer
  call initialize_problem(problem,2)
  problem%objective => rosenbrock_callback
  problem%has_gradient=.true.;problem%has_hessian=.true.
  call optimr(problem,[-1.2_dp,1.0_dp],'Rvmmin',result=answer)
  write(*,'(a,2f12.7)')'Rosenbrock solution: ',answer%par
  write(*,'(a,es12.4)')'objective: ',answer%value
  write(*,'(a,l1,1x,l1)')'KKT checks: ',answer%kkt1,answer%kkt2
end program demo_optimx
