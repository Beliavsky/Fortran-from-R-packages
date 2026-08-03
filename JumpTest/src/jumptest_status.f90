! SPDX-License-Identifier: MIT
module jumptest_status
  implicit none
  private

  integer, parameter, public :: JT_SUCCESS = 0
  integer, parameter, public :: JT_INVALID_ARGUMENT = 1
  integer, parameter, public :: JT_INVALID_DIMENSION = 2
  integer, parameter, public :: JT_NONFINITE_INPUT = 3
  integer, parameter, public :: JT_DEGENERATE_SAMPLE = 4
  integer, parameter, public :: JT_NUMERICAL_FAILURE = 5

  public :: status_message

contains

  pure function status_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message

    select case (status)
    case (JT_SUCCESS)
      message = 'success'
    case (JT_INVALID_ARGUMENT)
      message = 'invalid argument'
    case (JT_INVALID_DIMENSION)
      message = 'invalid dimension'
    case (JT_NONFINITE_INPUT)
      message = 'nonfinite input'
    case (JT_DEGENERATE_SAMPLE)
      message = 'degenerate sample'
    case (JT_NUMERICAL_FAILURE)
      message = 'numerical failure'
    case default
      message = 'unknown status'
    end select
  end function status_message

end module jumptest_status
