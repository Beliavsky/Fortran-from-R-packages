! SPDX-License-Identifier: GPL-2.0-only
program bounded_problem
  use optimx_mod
  use optimx_example_functions, only: quadratic_callback
  implicit none
  type(optimx_problem) :: problem
  type(optimx_result) :: answer
  call initialize_problem(problem,2)
  problem%objective => quadratic_callback
  problem%has_gradient=.true.;problem%has_hessian=.true.
  problem%lower=[0.0_dp,-1.0_dp];problem%upper=[0.5_dp,1.5_dp]
  call rvmminb(problem,[-2.0_dp,4.0_dp],result=answer)
  write(*,'(a,2f12.6)')'solution: ',answer%par
  write(*,'(a,es12.4)')'objective: ',answer%value
end program bounded_problem
