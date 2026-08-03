! SPDX-License-Identifier: GPL-2.0-or-later
program linear_constraint
  use maxlik, only: dp, maxlik_problem, maxlik_control, maxlik_result, initialize_problem, &
    set_equality_constraints, max_lik
  implicit none

  type(maxlik_problem) :: problem
  type(maxlik_control) :: control
  type(maxlik_result) :: result
  real(dp) :: a(1, 2), b(1)
  integer :: status

  call initialize_problem(problem, 2, objective)
  problem%gradient => gradient
  problem%hessian => hessian
  a(1, :) = [1.0_dp, 1.0_dp]
  b = 0.0_dp
  call set_equality_constraints(problem, a, b, status)
  control%constraint_tol = 1.0e-7_dp
  call max_lik(problem, [0.0_dp, 0.0_dp], result, 'nr', control)

  print '(a,2f12.6)', 'constrained estimate:', result%estimate
  print '(a,es12.4)', 'constraint violation:', result%constraint_violation

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

end program linear_constraint
