! SPDX-License-Identifier: GPL-2.0-or-later
program test_derivatives
  use maxlik, only: dp, maxlik_problem, derivative_comparison, initialize_problem, compare_derivatives
  implicit none

  type(maxlik_problem) :: problem
  type(derivative_comparison) :: comparison
  real(dp) :: x(2)

  x = [0.4_dp, -0.7_dp]
  call initialize_problem(problem, 2, objective)
  problem%gradient => gradient
  problem%hessian => hessian
  call compare_derivatives(problem, x, comparison, 2.0e-5_dp)

  call assert_true(comparison%status == 0, 'derivative comparison status')
  call assert_true(comparison%gradient_ok, 'analytic gradient')
  call assert_true(comparison%hessian_ok, 'analytic Hessian')
  call assert_close(comparison%max_gradient_error, 0.0_dp, 2.0e-5_dp, 'gradient error')
  call assert_close(comparison%max_hessian_error, 0.0_dp, 2.0e-4_dp, 'Hessian error')
  print '(a)', 'test_derivatives: PASS'

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    value = -0.5_dp * ((x(1) - 1.0_dp)**2 + 2.0_dp * (x(2) + 2.0_dp)**2)
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    g = [-(x(1) - 1.0_dp), -2.0_dp * (x(2) + 2.0_dp)]
    status = 0
  end subroutine gradient

  subroutine hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    h = 0.0_dp
    h(1, 1) = -1.0_dp + 0.0_dp * x(1)
    h(2, 2) = -2.0_dp
    status = 0
  end subroutine hessian

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      print '(a)', 'FAIL: ' // label
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call assert_true(abs(actual - expected) <= tolerance, label)
  end subroutine assert_close

end program test_derivatives
