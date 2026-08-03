! SPDX-License-Identifier: MIT
module mfgarch_status
  implicit none
  private

  integer, parameter, public :: mfgarch_success = 0
  integer, parameter, public :: mfgarch_invalid_argument = 1
  integer, parameter, public :: mfgarch_dimension_error = 2
  integer, parameter, public :: mfgarch_numerical_error = 3
  integer, parameter, public :: mfgarch_not_converged = 4
  integer, parameter, public :: mfgarch_singular_matrix = 5

  public :: mfgarch_message

contains

  pure function mfgarch_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message

    select case (status)
    case (mfgarch_success)
      message = 'success'
    case (mfgarch_invalid_argument)
      message = 'invalid argument'
    case (mfgarch_dimension_error)
      message = 'dimension mismatch'
    case (mfgarch_numerical_error)
      message = 'numerical error'
    case (mfgarch_not_converged)
      message = 'optimizer did not converge'
    case (mfgarch_singular_matrix)
      message = 'singular matrix'
    case default
      message = 'unknown status'
    end select
  end function mfgarch_message

end module mfgarch_status
