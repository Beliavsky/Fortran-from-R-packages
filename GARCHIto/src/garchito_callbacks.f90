! SPDX-License-Identifier: GPL-3.0-only
module garchito_callbacks
   use garchito_kinds, only : dp
   implicit none
   private

   public :: objective_callback, projection_callback

   abstract interface
      subroutine objective_callback(x, value, data)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: value
         class(*), intent(in) :: data
      end subroutine objective_callback

      subroutine projection_callback(x, data)
         import dp
         real(dp), intent(inout) :: x(:)
         class(*), intent(in) :: data
      end subroutine projection_callback
   end interface
end module garchito_callbacks
