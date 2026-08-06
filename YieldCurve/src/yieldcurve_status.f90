! SPDX-License-Identifier: GPL-2.0-or-later
module yieldcurve_status
   implicit none
   private

   integer, parameter, public :: yc_success = 0
   integer, parameter, public :: yc_invalid_argument = 1
   integer, parameter, public :: yc_dimension_error = 2
   integer, parameter, public :: yc_rank_deficient = 3
   integer, parameter, public :: yc_no_solution = 4

   public :: set_yc_status

contains

   subroutine set_yc_status(code, text, stat, message)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: stat
      character(len=*), intent(out), optional :: message

      if (present(stat)) stat = code
      if (present(message)) message = text
   end subroutine set_yc_status

end module yieldcurve_status
