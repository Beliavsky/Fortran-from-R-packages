! SPDX-License-Identifier: GPL-2.0-only
program test_solvers
  use optimx_mod
  use optimx_example_functions, only: rosenbrock_callback, quadratic_callback
  implicit none
  type(optimx_problem) :: problem
  type(optimx_control) :: control
  type(optimx_result) :: result
  real(dp) :: x0(2)

  call initialize_problem(problem, 2)
  problem%objective => rosenbrock_callback
  problem%has_gradient = .true.
  problem%has_hessian = .true.
  control = ctrldefault(2)
  control%maxit = 1000
  control%maxfeval = 20000
  x0 = [-1.2_dp, 1.0_dp]

  call rvmmin(problem, x0, control, result)
  call check(result%value < 1.0e-10_dp, 'Rvmmin Rosenbrock value')
  call check(maxval(abs(result%par - 1.0_dp)) < 2.0e-5_dp, 'Rvmmin Rosenbrock parameters')

  call rcgmin(problem, x0, control, result)
  call check(result%value < 1.0e-7_dp, 'Rcgmin Rosenbrock value')

  call snewton(problem, x0, control, result)
  call check(result%value < 1.0e-10_dp, 'snewton Rosenbrock value')

  problem%objective => quadratic_callback
  problem%lower = [0.0_dp, -1.0_dp]
  problem%upper = [0.5_dp, 1.5_dp]
  x0 = [-2.0_dp, 4.0_dp]
  call rvmminb(problem, x0, control, result)
  call check(maxval(abs(result%par - [0.5_dp, 1.5_dp])) < 2.0e-6_dp, 'bounded optimum')

  write(*,'(a)') 'test_solvers: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a,a)') 'FAIL: ', label
      error stop 1
    end if
  end subroutine check
end program test_solvers
