! SPDX-License-Identifier: GPL-2.0-or-later
module fincal_status
   implicit none
   private

   integer, parameter, public :: fincal_ok = 0
   integer, parameter, public :: fincal_invalid_input = 1
   integer, parameter, public :: fincal_size_mismatch = 2
   integer, parameter, public :: fincal_no_root = 3
   integer, parameter, public :: fincal_insufficient_inventory = 4
   integer, parameter, public :: fincal_nonfinite_result = 5
   integer, parameter, public :: fincal_weights_not_unit = 6

   public :: fincal_status_message
contains
   pure function fincal_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message

      select case (status)
      case (fincal_ok)
         message = 'success'
      case (fincal_invalid_input)
         message = 'invalid input'
      case (fincal_size_mismatch)
         message = 'array sizes do not conform'
      case (fincal_no_root)
         message = 'no root found in the search interval'
      case (fincal_insufficient_inventory)
         message = 'inventory is insufficient for the requested sale'
      case (fincal_nonfinite_result)
         message = 'calculation produced a nonfinite result'
      case (fincal_weights_not_unit)
         message = 'portfolio weights do not sum to one'
      case default
         message = 'unknown status'
      end select
   end function fincal_status_message
end module fincal_status
