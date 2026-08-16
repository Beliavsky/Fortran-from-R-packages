! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_status
   implicit none
   private
   integer, parameter, public :: pracma_ok = 0
   integer, parameter, public :: pracma_invalid_argument = 1
   integer, parameter, public :: pracma_dimension_mismatch = 2
   integer, parameter, public :: pracma_singular = 3
   integer, parameter, public :: pracma_not_converged = 4
   integer, parameter, public :: pracma_not_bracketed = 5
   integer, parameter, public :: pracma_nonfinite = 6
   integer, parameter, public :: pracma_max_iterations = 7
   integer, parameter, public :: pracma_not_positive_definite = 8
   integer, parameter, public :: pracma_infeasible = 9
   integer, parameter, public :: pracma_unsupported = 10
   public :: pracma_status_message
contains
   pure function pracma_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message
      select case (status)
      case (pracma_ok);                    message = 'success'
      case (pracma_invalid_argument);      message = 'invalid argument'
      case (pracma_dimension_mismatch);    message = 'dimension mismatch'
      case (pracma_singular);              message = 'singular matrix or system'
      case (pracma_not_converged);         message = 'iteration did not converge'
      case (pracma_not_bracketed);         message = 'root is not bracketed'
      case (pracma_nonfinite);             message = 'nonfinite value encountered'
      case (pracma_max_iterations);        message = 'maximum iterations reached'
      case (pracma_not_positive_definite); message = 'matrix is not positive definite'
      case (pracma_infeasible);            message = 'problem is infeasible'
      case (pracma_unsupported);           message = 'operation is not supported'
      case default;                        message = 'unknown status'
      end select
   end function pracma_status_message
end module pracma_status
