! SPDX-License-Identifier: GPL-2.0-or-later
program test_constraints_and_fixed
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, &
    set_equality_constraints, set_fixed, max_lik
  implicit none

  type(maxlik_problem) :: constrained_problem, fixed_problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  real(dp) :: a(1, 2), b(1), start(2)
  integer :: status

  start = [0.0_dp, 0.0_dp]
  call initialize_problem(constrained_problem, 2, objective)
  constrained_problem%gradient => gradient
  constrained_problem%hessian => hessian
  a(1, :) = [1.0_dp, 1.0_dp]
  b = 0.0_dp
  call set_equality_constraints(constrained_problem, a, b, status)
  call assert_true(status == 0, 'set equality constraints')
  control%iterlim = 100
  control%constraint_tol = 1.0e-7_dp
  control%constraint_max_outer = 16
  call max_lik(constrained_problem, start, result, 'nr', control)
  call assert_true(result%constraint_violation <= 1.0e-7_dp, 'constraint violation')
  call assert_close(result%estimate(1), 1.5_dp, 2.0e-5_dp, 'constrained x1')
  call assert_close(result%estimate(2), -1.5_dp, 2.0e-5_dp, 'constrained x2')

  call initialize_problem(fixed_problem, 2, objective)
  fixed_problem%gradient => gradient
  fixed_problem%hessian => hessian
  call set_fixed(fixed_problem, [2])
  start = [-4.0_dp, 5.0_dp]
  call max_lik(fixed_problem, start, result, 'bfgs', control)
  call assert_true(result%converged, 'fixed-parameter fit converged')
  call assert_close(result%estimate(1), 2.0_dp, 1.0e-6_dp, 'free parameter')
  call assert_close(result%estimate(2), 5.0_dp, 0.0_dp, 'fixed parameter')
  call assert_close(result%covariance(2, 2), 0.0_dp, 0.0_dp, 'fixed covariance')
  print '(a)', 'test_constraints_and_fixed: PASS'

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    value = -0.5_dp * ((x(1) - 2.0_dp)**2 + (x(2) + 1.0_dp)**2)
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    g = [-(x(1) - 2.0_dp), -(x(2) + 1.0_dp)]
    status = 0
  end subroutine gradient

  subroutine hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    h = 0.0_dp
    h(1, 1) = -1.0_dp + 0.0_dp * x(1)
    h(2, 2) = -1.0_dp
    status = 0
  end subroutine hessian

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAIL: ' // label
      if (allocated(result%estimate)) print '(a,2f14.6)', 'estimate:', result%estimate
      print '(a,es14.5)', 'constraint violation:', result%constraint_violation
      print '(a,i0,2a)', 'code=', result%code, ' message=', trim(result%message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual - expected) <= tolerance, label)
  end subroutine assert_close

end program test_constraints_and_fixed
