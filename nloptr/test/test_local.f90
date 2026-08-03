! SPDX-License-Identifier: LGPL-3.0-or-later
program test_local
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
  options = nl_opts()
  options%maxeval = 5000
  options%xtol_rel = 1.0e-9_dp
  options%ftol_rel = 1.0e-12_dp

  call lbfgs(problem, [-1.2_dp, 1.0_dp], options, result)
  call check(maxval(abs(result%solution - [1.0_dp, 1.0_dp])) < 2.0e-5_dp, 'lbfgs solution')
  call check(result%objective < 1.0e-9_dp, 'lbfgs objective')

  call neldermead(problem, [-1.2_dp, 1.0_dp], options, result)
  call check(maxval(abs(result%solution - [1.0_dp, 1.0_dp])) < 2.0e-4_dp, 'nelder-mead solution')
  call check(result%objective < 1.0e-7_dp, 'nelder-mead objective')

  call bobyqa(problem, [-1.2_dp, 1.0_dp], options, result)
  call check(result%objective < 2.0e-4_dp, 'bobyqa adapted objective')

  print '(a)', 'test_local: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,a)') 'FAIL: ', label
      error stop 1
    end if
  end subroutine check
end program test_local
