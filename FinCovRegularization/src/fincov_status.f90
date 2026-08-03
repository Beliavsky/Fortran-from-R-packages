! SPDX-License-Identifier: GPL-2.0-only
module fincov_status
   implicit none
   private

   integer, parameter, public :: fincov_ok = 0
   integer, parameter, public :: fincov_invalid_input = 1
   integer, parameter, public :: fincov_size_mismatch = 2
   integer, parameter, public :: fincov_singular_matrix = 3
   integer, parameter, public :: fincov_no_convergence = 4
   integer, parameter, public :: fincov_not_positive_definite = 5
   integer, parameter, public :: fincov_allocation_failure = 6

   public :: fincov_status_message
contains
   pure function fincov_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message

      select case (status)
      case (fincov_ok)
         message = "ok"
      case (fincov_invalid_input)
         message = "invalid input"
      case (fincov_size_mismatch)
         message = "size mismatch"
      case (fincov_singular_matrix)
         message = "singular or numerically rank-deficient matrix"
      case (fincov_no_convergence)
         message = "iteration did not converge"
      case (fincov_not_positive_definite)
         message = "matrix is not positive definite"
      case (fincov_allocation_failure)
         message = "allocation failure"
      case default
         message = "unknown status"
      end select
   end function fincov_status_message
end module fincov_status
