! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_status
  implicit none
  private
  integer, parameter, public :: MIXTOOLS_SUCCESS = 0
  integer, parameter, public :: MIXTOOLS_INVALID_ARGUMENT = 1
  integer, parameter, public :: MIXTOOLS_DIMENSION_ERROR = 2
  integer, parameter, public :: MIXTOOLS_SINGULAR = 3
  integer, parameter, public :: MIXTOOLS_NOT_POSITIVE_DEFINITE = 4
  integer, parameter, public :: MIXTOOLS_NOT_CONVERGED = 5
  integer, parameter, public :: MIXTOOLS_NUMERICAL_ERROR = 6
  integer, parameter, public :: MIXTOOLS_UNSUPPORTED = 7
  public :: mixtools_status_message
contains
  pure function mixtools_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message
    select case (status)
    case (MIXTOOLS_SUCCESS); message = "success"
    case (MIXTOOLS_INVALID_ARGUMENT); message = "invalid argument"
    case (MIXTOOLS_DIMENSION_ERROR); message = "dimension error"
    case (MIXTOOLS_SINGULAR); message = "singular matrix"
    case (MIXTOOLS_NOT_POSITIVE_DEFINITE); message = "matrix is not positive definite"
    case (MIXTOOLS_NOT_CONVERGED); message = "iteration limit reached"
    case (MIXTOOLS_NUMERICAL_ERROR); message = "numerical error"
    case (MIXTOOLS_UNSUPPORTED); message = "unsupported option"
    case default; message = "unknown status"
    end select
  end function mixtools_status_message
end module mixtools_status
