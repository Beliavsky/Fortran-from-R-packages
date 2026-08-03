! SPDX-License-Identifier: LGPL-3.0-or-later
module nloptr_example_functions
  use nloptr_kinds, only: dp
  implicit none
  private
  public :: rosenbrock_objective, quadratic_objective, tutorial_objective
  public :: tutorial_constraints, equality_constraint, multimodal_objective
  public :: quadratic_value_only
contains
  subroutine rosenbrock_objective(x, value, gradient, need_gradient, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(inout) :: gradient(:)
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    value = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2
    gradient = 0.0_dp
    if (need_gradient) then
      gradient(1) = -2.0_dp * (1.0_dp - x(1)) - 400.0_dp * x(1) * (x(2) - x(1)**2)
      gradient(2) = 200.0_dp * (x(2) - x(1)**2)
    end if
    status = 0
  end subroutine rosenbrock_objective

  subroutine quadratic_objective(x, value, gradient, need_gradient, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(inout) :: gradient(:)
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    value = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
    gradient = 0.0_dp
    if (need_gradient) gradient = [2.0_dp * (x(1) - 1.0_dp), 2.0_dp * (x(2) - 2.0_dp)]
    status = 0
  end subroutine quadratic_objective


  subroutine quadratic_value_only(x, value, gradient, need_gradient, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(inout) :: gradient(:)
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    value = (x(1) - 1.0_dp)**2 + (x(2) - 2.0_dp)**2
    if (.not. need_gradient) gradient = 0.0_dp
    status = 0
  end subroutine quadratic_value_only

  subroutine tutorial_objective(x, value, gradient, need_gradient, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(inout) :: gradient(:)
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    if (x(2) <= 0.0_dp) then
      value = huge(1.0_dp)
      gradient = 0.0_dp
      status = -5
      return
    end if
    value = sqrt(x(2))
    gradient = 0.0_dp
    if (need_gradient) gradient(2) = 0.5_dp / sqrt(x(2))
    status = 0
  end subroutine tutorial_objective

  subroutine tutorial_constraints(x, values, jacobian, need_jacobian, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: values(:)
    real(dp), intent(inout) :: jacobian(:, :)
    logical, intent(in) :: need_jacobian
    integer, intent(out) :: status
    real(dp), parameter :: a(2) = [2.0_dp, -1.0_dp]
    real(dp), parameter :: b(2) = [0.0_dp, 1.0_dp]
    integer :: i
    values = (a * x(1) + b)**3 - x(2)
    jacobian = 0.0_dp
    if (need_jacobian) then
      do i = 1, 2
        jacobian(i, 1) = 3.0_dp * a(i) * (a(i) * x(1) + b(i))**2
        jacobian(i, 2) = -1.0_dp
      end do
    end if
    status = 0
  end subroutine tutorial_constraints

  subroutine equality_constraint(x, values, jacobian, need_jacobian, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: values(:)
    real(dp), intent(inout) :: jacobian(:, :)
    logical, intent(in) :: need_jacobian
    integer, intent(out) :: status
    values(1) = x(1) + x(2) - 1.0_dp
    jacobian = 0.0_dp
    if (need_jacobian) jacobian(1, :) = 1.0_dp
    status = 0
  end subroutine equality_constraint

  subroutine multimodal_objective(x, value, gradient, need_gradient, status)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), intent(inout) :: gradient(:)
    logical, intent(in) :: need_gradient
    integer, intent(out) :: status
    value = 0.2_dp * (x(1)**2 + x(2)**2) - cos(3.0_dp * x(1)) - cos(3.0_dp * x(2))
    gradient = 0.0_dp
    if (need_gradient) then
      gradient(1) = 0.4_dp * x(1) + 3.0_dp * sin(3.0_dp * x(1))
      gradient(2) = 0.4_dp * x(2) + 3.0_dp * sin(3.0_dp * x(2))
    end if
    status = 0
  end subroutine multimodal_objective
end module nloptr_example_functions
