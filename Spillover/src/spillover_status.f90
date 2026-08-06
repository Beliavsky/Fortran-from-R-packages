! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_status
   implicit none
   private

   integer, parameter, public :: spillover_success = 0
   integer, parameter, public :: spillover_invalid_argument = 1
   integer, parameter, public :: spillover_singular_matrix = 2
   integer, parameter, public :: spillover_not_positive_definite = 3
   integer, parameter, public :: spillover_allocation_error = 4
   integer, parameter, public :: spillover_iteration_limit = 5

   public :: set_status

contains

   subroutine set_status(info, message, code, text)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message
      integer, intent(in) :: code
      character(len=*), intent(in) :: text

      if (present(info)) info = code
      if (present(message)) message = text
   end subroutine set_status

end module spillover_status
