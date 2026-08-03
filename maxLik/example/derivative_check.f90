! SPDX-License-Identifier: GPL-2.0-or-later
program derivative_check
  use maxlik, only: dp, maxlik_problem, derivative_comparison, initialize_problem, compare_derivatives
  implicit none

  type(maxlik_problem) :: problem
  type(derivative_comparison) :: report

  call initialize_problem(problem, 2, objective)
  problem%gradient => gradient
  problem%hessian => hessian
  call compare_derivatives(problem, [0.3_dp, -0.4_dp], report)

  print '(a,es12.4)', 'maximum gradient error:', report%max_gradient_error
  print '(a,es12.4)', 'maximum Hessian error: ', report%max_hessian_error

contains

  subroutine objective(x, value, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    value = -x(1)**2 - 2.0_dp * x(2)**2 + x(1) * x(2)
    status = 0
  end subroutine objective

  subroutine gradient(x, g, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(out) :: status
    g = [-2.0_dp * x(1) + x(2), -4.0_dp * x(2) + x(1)]
    status = 0
  end subroutine gradient

  subroutine hessian(x, h, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(out) :: status
    h = reshape([-2.0_dp + 0.0_dp * x(1), 1.0_dp, 1.0_dp, -4.0_dp], [2, 2])
    status = 0
  end subroutine hessian

end program derivative_check
