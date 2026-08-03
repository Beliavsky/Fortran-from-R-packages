! SPDX-License-Identifier: GPL-2.0-or-later
module fints_status
   implicit none
   private
   integer, parameter, public :: fints_ok = 0
   integer, parameter, public :: fints_invalid_input = 1
   integer, parameter, public :: fints_singular = 2
   integer, parameter, public :: fints_nonstationary = 3
   integer, parameter, public :: fints_iteration_limit = 4
   integer, parameter, public :: fints_no_data = 5
   integer, parameter, public :: fints_numerical_failure = 6
   integer, parameter, public :: fints_io_error = 7
   public :: fints_status_message
contains
   pure function fints_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message

      select case (status)
      case (fints_ok)
         message = 'ok'
      case (fints_invalid_input)
         message = 'invalid input'
      case (fints_singular)
         message = 'singular system'
      case (fints_nonstationary)
         message = 'nonstationary autoregressive model'
      case (fints_iteration_limit)
         message = 'iteration limit reached'
      case (fints_no_data)
         message = 'no usable data'
      case (fints_numerical_failure)
         message = 'numerical failure'
      case (fints_io_error)
         message = 'input/output error'
      case default
         message = 'unknown status'
      end select
   end function fints_status_message
end module fints_status
