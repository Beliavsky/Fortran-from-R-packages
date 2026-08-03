! SPDX-License-Identifier: GPL-2.0-only
program rosenbrock_methods
  use optimx_mod
  use optimx_example_functions, only: rosenbrock_callback
  implicit none
  type(optimx_problem) :: problem
  type(optimx_multi_result) :: answer
  character(len=16) :: methods(4)
  integer :: i
  call initialize_problem(problem,2)
  problem%objective => rosenbrock_callback
  problem%has_gradient=.true.;problem%has_hessian=.true.
  methods=[character(len=16)::'Rvmmin','Rcgmin','Nelder-Mead','snewton']
  call optimx(problem,[-1.2_dp,1.0_dp],methods,result=answer)
  do i=1,size(answer%runs)
    write(*,'(a16,1x,es12.4,2(1x,f10.6))')trim(answer%runs(i)%method), &
      answer%runs(i)%value,answer%runs(i)%par
  end do
end program rosenbrock_methods
