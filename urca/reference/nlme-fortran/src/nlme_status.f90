! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_status
  implicit none
  private
  integer, parameter, public :: NLME_SUCCESS = 0
  integer, parameter, public :: NLME_INVALID_ARGUMENT = 1
  integer, parameter, public :: NLME_SINGULAR = 2
  integer, parameter, public :: NLME_NOT_POSITIVE_DEFINITE = 3
  integer, parameter, public :: NLME_MAX_ITER = 4
  integer, parameter, public :: NLME_NONFINITE = 5
  integer, parameter, public :: NLME_CALLBACK_ERROR = 6
  integer, parameter, public :: NLME_DIMENSION_ERROR = 7
  public :: nlme_status_message
contains
  pure function nlme_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message
    select case (status)
    case (NLME_SUCCESS)
      message = 'success'
    case (NLME_INVALID_ARGUMENT)
      message = 'invalid argument'
    case (NLME_SINGULAR)
      message = 'singular system'
    case (NLME_NOT_POSITIVE_DEFINITE)
      message = 'matrix is not positive definite'
    case (NLME_MAX_ITER)
      message = 'maximum iterations reached'
    case (NLME_NONFINITE)
      message = 'nonfinite value encountered'
    case (NLME_CALLBACK_ERROR)
      message = 'model callback failed'
    case (NLME_DIMENSION_ERROR)
      message = 'inconsistent dimensions'
    case default
      message = 'unknown status'
    end select
  end function nlme_status_message
end module nlme_status
