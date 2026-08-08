! SPDX-License-Identifier: MIT
module gradient_benchmarks
  use gradient_kinds, only : dp
  implicit none
  private
  public :: sphere_objective, rosenbrock_objective, shifted_sphere_objective
contains
  function sphere_objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum(x*x)
  end function sphere_objective

  function shifted_sphere_objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = sum((x-1.0_dp)**2)
  end function shifted_sphere_objective

  function rosenbrock_objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    integer :: i
    f = 0.0_dp
    do i = 1, size(x)-1
      f = f + 100.0_dp*(x(i+1)-x(i)*x(i))**2 + (1.0_dp-x(i))**2
    end do
  end function rosenbrock_objective
end module gradient_benchmarks
