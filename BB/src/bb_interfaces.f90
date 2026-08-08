! SPDX-License-Identifier: GPL-3.0-only
module bb_interfaces
  use bb_kinds, only: dp
  implicit none
  private

  public :: bb_scalar_fn, bb_gradient_fn, bb_vector_fn, bb_projection_fn

  abstract interface
    function bb_scalar_fn(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function bb_scalar_fn

    subroutine bb_gradient_fn(x, gradient)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: gradient(:)
    end subroutine bb_gradient_fn

    subroutine bb_vector_fn(x, value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: value(:)
    end subroutine bb_vector_fn

    subroutine bb_projection_fn(x, projected, ok)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: projected(:)
      logical, intent(out) :: ok
    end subroutine bb_projection_fn
  end interface
end module bb_interfaces
