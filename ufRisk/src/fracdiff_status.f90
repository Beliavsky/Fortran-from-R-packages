! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_status
   implicit none
   private

   integer, parameter, public :: fd_ok = 0
   integer, parameter, public :: fd_invalid_input = 1
   integer, parameter, public :: fd_gamma_error = 2
   integer, parameter, public :: fd_optimization_failed = 3
   integer, parameter, public :: fd_iteration_limit = 4
   integer, parameter, public :: fd_singular_hessian = 5
   integer, parameter, public :: fd_unstable_ar = 6
   integer, parameter, public :: fd_insufficient_data = 7

   public :: fracdiff_status_message

contains

   pure function fracdiff_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message

      select case (status)
      case (fd_ok)
         message = "ok"
      case (fd_invalid_input)
         message = "invalid input"
      case (fd_gamma_error)
         message = "gamma function or fractional recursion failure"
      case (fd_optimization_failed)
         message = "optimization failed"
      case (fd_iteration_limit)
         message = "optimization iteration limit reached"
      case (fd_singular_hessian)
         message = "singular or indefinite Hessian"
      case (fd_unstable_ar)
         message = "autoregressive polynomial is not stationary"
      case (fd_insufficient_data)
         message = "insufficient data"
      case default
         message = "unknown status"
      end select
   end function fracdiff_status_message

end module fracdiff_status
