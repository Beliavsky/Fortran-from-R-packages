! SPDX-License-Identifier: GPL-3.0-only
module fingraph_status
   implicit none
   private
   integer, parameter, public :: fg_ok = 0
   integer, parameter, public :: fg_invalid_input = 1
   integer, parameter, public :: fg_size_mismatch = 2
   integer, parameter, public :: fg_singular_matrix = 3
   integer, parameter, public :: fg_no_convergence = 4
   integer, parameter, public :: fg_not_positive_definite = 5
   integer, parameter, public :: fg_allocation_failure = 6
   integer, parameter, public :: fg_infeasible = 7
   public :: fg_status_message
contains
   pure function fg_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message
      select case (status)
      case (fg_ok); message = "ok"
      case (fg_invalid_input); message = "invalid input"
      case (fg_size_mismatch); message = "size mismatch"
      case (fg_singular_matrix); message = "singular or rank-deficient matrix"
      case (fg_no_convergence); message = "iteration did not converge"
      case (fg_not_positive_definite); message = "matrix is not positive definite"
      case (fg_allocation_failure); message = "allocation failure"
      case (fg_infeasible); message = "infeasible constraints"
      case default; message = "unknown status"
      end select
   end function fg_status_message
end module fingraph_status
