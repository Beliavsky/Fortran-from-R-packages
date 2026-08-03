! SPDX-License-Identifier: LGPL-3.0-or-later
program test_global_and_api
  use nloptr_mod
  use nloptr_example_functions, only: multimodal_objective, rosenbrock_objective
  implicit none
  type(nloptr_problem) :: problem
  type(nloptr_options) :: options, defaults
  type(nloptr_result) :: result

  defaults = nloptr_get_default_options()
  call check(defaults%maxeval > 0, 'default options')
  problem%n = 2
  problem%objective => multimodal_objective
  allocate(problem%lower(2), problem%upper(2))
  problem%lower = -3.0_dp
  problem%upper = 3.0_dp
  options = nl_opts()
  options%population = 80
  options%maxeval = 4000
  options%local_algorithm = 'NLOPT_LD_LBFGS'
  call direct(problem, [2.5_dp, -2.0_dp], options, result)
  call check(result%objective < -1.99_dp, 'global direct adapted objective')
  call check(maxval(abs(result%solution)) < 2.0e-3_dp, 'global solution')
  call check(is_nloptr(result), 'is_nloptr')

  problem%objective => rosenbrock_objective
  options%maxeval = 5000
  call mlsl(problem, [-1.2_dp, 1.0_dp], options, result)
  call check(result%objective < 1.0e-7_dp, 'mlsl wrapper')

  print '(a)', 'test_global_and_api: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,a)') 'FAIL: ', label
      write(*, '(a,es20.10)') 'objective: ', result%objective
      if (allocated(result%solution)) write(*, '(a,*(es16.8,1x))') 'solution: ', result%solution
      error stop 1
    end if
  end subroutine check
end program test_global_and_api
