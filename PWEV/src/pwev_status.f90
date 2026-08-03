! SPDX-License-Identifier: GPL-3.0-only
module pwev_status
  implicit none
  private
  integer, parameter, public :: PWEV_SUCCESS = 0
  integer, parameter, public :: PWEV_INVALID_INPUT = 1
  integer, parameter, public :: PWEV_MODEL_FAILURE = 2
  integer, parameter, public :: PWEV_OPTIMIZER_FAILURE = 3
  integer, parameter, public :: PWEV_NUMERICAL_ERROR = 4
  integer, parameter, public :: PWEV_GARCH_MEAN = 1
  integer, parameter, public :: PWEV_GARCH_SIGMA = 2
  integer, parameter, public :: PWEV_MEM_UPSTREAM_OOS = 1
  integer, parameter, public :: PWEV_MEM_RECURSIVE_OOS = 2
  public :: pwev_status_message
contains
  pure function pwev_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=64) :: message
    select case (status)
    case (PWEV_SUCCESS)
      message = 'success'
    case (PWEV_INVALID_INPUT)
      message = 'invalid input'
    case (PWEV_MODEL_FAILURE)
      message = 'one or more base models failed'
    case (PWEV_OPTIMIZER_FAILURE)
      message = 'ensemble optimizer failed'
    case (PWEV_NUMERICAL_ERROR)
      message = 'numerical error'
    case default
      message = 'unknown status'
    end select
  end function pwev_status_message
end module pwev_status
