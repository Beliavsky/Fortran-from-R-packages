! SPDX-License-Identifier: GPL-3.0-only
module smoots_status
   implicit none
   private
   integer, parameter, public :: sm_ok = 0
   integer, parameter, public :: sm_invalid_input = 1
   integer, parameter, public :: sm_singular = 2
   integer, parameter, public :: sm_iteration_limit = 3
   integer, parameter, public :: sm_fit_failed = 4
   integer, parameter, public :: sm_allocation_failed = 5
   public :: smoots_status_message
contains
   pure function smoots_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=64) :: message
      select case (status)
      case (sm_ok); message = 'success'
      case (sm_invalid_input); message = 'invalid input'
      case (sm_singular); message = 'singular linear system'
      case (sm_iteration_limit); message = 'iteration limit reached'
      case (sm_fit_failed); message = 'model fitting failed'
      case (sm_allocation_failed); message = 'allocation failed'
      case default; message = 'unknown status'
      end select
   end function smoots_status_message
end module smoots_status
