! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_status
   implicit none
   integer, parameter, public :: icsnp_ok = 0
   integer, parameter, public :: icsnp_invalid_input = 1
   integer, parameter, public :: icsnp_singular = 2
   integer, parameter, public :: icsnp_iteration_limit = 3
   integer, parameter, public :: icsnp_numerical_error = 4
contains
   pure function icsnp_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message
      select case (status)
      case (icsnp_ok)
         message = "success"
      case (icsnp_invalid_input)
         message = "invalid input"
      case (icsnp_singular)
         message = "singular or non-positive-definite matrix"
      case (icsnp_iteration_limit)
         message = "iteration limit reached"
      case default
         message = "numerical error"
      end select
   end function icsnp_status_message
end module icsnp_status
