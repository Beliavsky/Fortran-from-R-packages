! SPDX-License-Identifier: GPL-3.0-only
module rumidas_status
  implicit none
  private
  integer, parameter, public :: RUMIDAS_SUCCESS = 0
  integer, parameter, public :: RUMIDAS_INVALID_INPUT = 1
  integer, parameter, public :: RUMIDAS_DIMENSION_ERROR = 2
  integer, parameter, public :: RUMIDAS_INVALID_PARAMETER = 3
  integer, parameter, public :: RUMIDAS_NUMERICAL_ERROR = 4
  integer, parameter, public :: RUMIDAS_OPTIMIZATION_ERROR = 5
  public :: rumidas_status_message
contains
  pure function rumidas_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=96) :: message
    select case (status)
    case (RUMIDAS_SUCCESS)
      message = 'success'
    case (RUMIDAS_INVALID_INPUT)
      message = 'invalid input'
    case (RUMIDAS_DIMENSION_ERROR)
      message = 'inconsistent array dimensions'
    case (RUMIDAS_INVALID_PARAMETER)
      message = 'invalid model parameter'
    case (RUMIDAS_NUMERICAL_ERROR)
      message = 'numerical failure'
    case (RUMIDAS_OPTIMIZATION_ERROR)
      message = 'optimization failed'
    case default
      message = 'unknown status'
    end select
  end function rumidas_status_message
end module rumidas_status
