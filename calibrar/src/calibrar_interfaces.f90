! SPDX-License-Identifier: GPL-2.0-only
module calibrar_interfaces
  use calibrar_kinds, only : dp
  implicit none
  private
  public :: scalar_objective, gradient_callback, vector_objective

  abstract interface
    function scalar_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function scalar_objective

    subroutine gradient_callback(x, g)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
    end subroutine gradient_callback

    subroutine vector_objective(x, f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: f(:)
    end subroutine vector_objective
  end interface
end module calibrar_interfaces
