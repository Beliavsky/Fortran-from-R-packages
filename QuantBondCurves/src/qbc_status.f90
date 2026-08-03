! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_status
   implicit none
   private
   integer, parameter, public :: qbc_success = 0
   integer, parameter, public :: qbc_invalid_argument = 1
   integer, parameter, public :: qbc_size_mismatch = 2
   integer, parameter, public :: qbc_no_convergence = 3
   integer, parameter, public :: qbc_singular = 4
   integer, parameter, public :: qbc_infeasible = 5
   public :: qbc_status_message
contains
   pure function qbc_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=80) :: message
      select case (status)
      case (qbc_success); message = 'success'
      case (qbc_invalid_argument); message = 'invalid argument'
      case (qbc_size_mismatch); message = 'size mismatch'
      case (qbc_no_convergence); message = 'iteration did not converge'
      case (qbc_singular); message = 'singular numerical problem'
      case (qbc_infeasible); message = 'infeasible problem'
      case default; message = 'unknown status'
      end select
   end function qbc_status_message
end module qbc_status
