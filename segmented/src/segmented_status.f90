! SPDX-License-Identifier: GPL-2.0-or-later
module segmented_status
  implicit none
  private
  integer, parameter, public :: SEG_SUCCESS = 0
  integer, parameter, public :: SEG_INVALID_ARGUMENT = 1
  integer, parameter, public :: SEG_SINGULAR = 2
  integer, parameter, public :: SEG_MAX_ITER = 3
  integer, parameter, public :: SEG_NONFINITE = 4
  integer, parameter, public :: SEG_DIMENSION_ERROR = 5
  integer, parameter, public :: SEG_INFEASIBLE_BREAKPOINT = 6
  integer, parameter, public :: SEG_NLME_ERROR = 7
  public :: segmented_status_message
contains
  pure function segmented_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message
    select case (status)
    case (SEG_SUCCESS)
      message = 'success'
    case (SEG_INVALID_ARGUMENT)
      message = 'invalid argument'
    case (SEG_SINGULAR)
      message = 'singular regression system'
    case (SEG_MAX_ITER)
      message = 'maximum iterations reached'
    case (SEG_NONFINITE)
      message = 'nonfinite value encountered'
    case (SEG_DIMENSION_ERROR)
      message = 'inconsistent dimensions'
    case (SEG_INFEASIBLE_BREAKPOINT)
      message = 'breakpoint outside admissible range'
    case (SEG_NLME_ERROR)
      message = 'mixed-effects fit failed'
    case default
      message = 'unknown status'
    end select
  end function segmented_status_message
end module segmented_status
