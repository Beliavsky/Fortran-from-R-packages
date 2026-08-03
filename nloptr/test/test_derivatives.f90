! SPDX-License-Identifier: LGPL-3.0-or-later
program test_derivatives
  use nloptr_mod
  use nloptr_example_functions, only: rosenbrock_objective, tutorial_constraints, quadratic_value_only
  implicit none
  real(dp) :: x(2), grad(2), analytic(2), jac(2, 2), values(2), dummy_jac(2, 2)
  type(derivative_check_result) :: result
  integer :: status

  x = [-1.2_dp, 1.0_dp]
  call nl_grad(x, rosenbrock_objective, grad, status)
  call rosenbrock_objective(x, values(1), analytic, .true., status)
  call check(maxval(abs(grad - analytic)) < 2.0e-5_dp, 'gradient')

  call nl_jacobian(x, 2, tutorial_constraints, jac, status)
  call tutorial_constraints(x, values, dummy_jac, .true., status)
  call check(maxval(abs(jac - dummy_jac)) < 2.0e-5_dp, 'jacobian')

  call check_derivatives(x, rosenbrock_objective, analytic, result, 1.0e-4_dp)
  call check(result%n_warnings == 0, 'derivative checker')
  analytic(1) = analytic(1) + 1.0_dp
  call check_derivatives(x, rosenbrock_objective, analytic, result, 1.0e-4_dp)
  call check(result%n_warnings == 1, 'derivative checker warning')

  call nl_grad([0.0_dp, 0.0_dp], quadratic_value_only, grad, status)
  call check(maxval(abs(grad - [-2.0_dp, -4.0_dp])) < 2.0e-5_dp, 'value-only callback')

  print '(a)', 'test_derivatives: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*, '(a,a)') 'FAIL: ', label
      error stop 1
    end if
  end subroutine check
end program test_derivatives
