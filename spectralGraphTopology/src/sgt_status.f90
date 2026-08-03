! SPDX-License-Identifier: GPL-3.0-only
module sgt_status
   implicit none
   private
   integer, parameter, public :: sgt_ok = 0
   integer, parameter, public :: sgt_invalid_input = 1
   integer, parameter, public :: sgt_size_mismatch = 2
   integer, parameter, public :: sgt_singular_matrix = 3
   integer, parameter, public :: sgt_no_convergence = 4
   integer, parameter, public :: sgt_not_positive_definite = 5
   integer, parameter, public :: sgt_allocation_failure = 6
   integer, parameter, public :: sgt_infeasible = 7
   public :: sgt_status_message
contains
   pure function sgt_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message
      select case (status)
      case (sgt_ok); message = "ok"
      case (sgt_invalid_input); message = "invalid input"
      case (sgt_size_mismatch); message = "size mismatch"
      case (sgt_singular_matrix); message = "singular or rank-deficient matrix"
      case (sgt_no_convergence); message = "iteration did not converge"
      case (sgt_not_positive_definite); message = "matrix is not positive definite"
      case (sgt_allocation_failure); message = "allocation failure"
      case (sgt_infeasible); message = "infeasible constraints"
      case default; message = "unknown status"
      end select
   end function sgt_status_message
end module sgt_status
