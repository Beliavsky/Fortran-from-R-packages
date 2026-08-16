! SPDX-License-Identifier: GPL-2.0-only
module rsolnp_callbacks
   use rsolnp_kinds, only : dp
   implicit none
   private

   public :: objective_callback, gradient_callback
   public :: vector_callback, jacobian_callback

   abstract interface
      subroutine objective_callback(x, value, data)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value
         class(*), intent(in), optional :: data
      end subroutine objective_callback

      subroutine gradient_callback(x, gradient, data)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: gradient(:)
         class(*), intent(in), optional :: data
      end subroutine gradient_callback

      subroutine vector_callback(x, value, data)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value(:)
         class(*), intent(in), optional :: data
      end subroutine vector_callback

      subroutine jacobian_callback(x, jacobian, data)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: jacobian(:, :)
         class(*), intent(in), optional :: data
      end subroutine jacobian_callback
   end interface

end module rsolnp_callbacks
