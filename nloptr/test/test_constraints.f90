! SPDX-License-Identifier: LGPL-3.0-or-later
program test_constraints
  use nloptr_mod
  use nloptr_example_functions, only: tutorial_objective, tutorial_constraints, &
    quadratic_objective, equality_constraint
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
  options = nl_opts()
  options%maxeval = 12000
  options%max_outer = 10
  options%constraint_tol = 2.0e-5_dp
  options%xtol_rel = 1.0e-8_dp
  call mma(problem, [1.234_dp, 5.678_dp], options, result)
  call check(maxval(abs(result%solution - [1.0_dp / 3.0_dp, 8.0_dp / 27.0_dp])) < 3.0e-3_dp, &
    'tutorial constrained solution')
  call check(result%max_constraint < 1.0e-4_dp, 'tutorial feasibility')

  nullify(problem%inequality)
  problem%n_ineq = 0
  problem%n_eq = 1
  problem%objective => quadratic_objective
  problem%equality => equality_constraint
  problem%lower = -10.0_dp
  problem%upper = 10.0_dp
  call slsqp(problem, [0.0_dp, 0.0_dp], options, result)
  call check(maxval(abs(result%solution - [0.0_dp, 1.0_dp])) < 2.0e-3_dp, 'equality solution')
  call check(result%max_constraint < 1.0e-5_dp, 'equality feasibility')

  print '(a)', 'test_constraints: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,a)') 'FAIL: ', label
      write(*, '(a,2es20.10)') 'solution: ', result%solution
      write(*, '(a,es20.10)') 'constraint: ', result%max_constraint
      error stop 1
    end if
  end subroutine check
end program test_constraints
