! SPDX-License-Identifier: GPL-3.0-only
module tscopula_status
  implicit none
  private
  integer, parameter, public :: tsc_success = 0
  integer, parameter, public :: tsc_invalid_input = 1
  integer, parameter, public :: tsc_numerical_failure = 2
  integer, parameter, public :: tsc_nonstationary = 3
  integer, parameter, public :: tsc_unsupported = 4
  integer, parameter, public :: tsc_not_converged = 5

  type, public :: tsc_error
    integer :: code = tsc_success
    character(len=:), allocatable :: message
  contains
    procedure :: ok => error_ok
  end type tsc_error

  public :: set_error, clear_error
contains
  pure logical function error_ok(this) result(ok)
    class(tsc_error), intent(in) :: this
    ok = this%code == tsc_success
  end function error_ok

  subroutine clear_error(error)
    type(tsc_error), intent(out) :: error
    error%code = tsc_success
    error%message = 'success'
  end subroutine clear_error

  subroutine set_error(error, code, message)
    type(tsc_error), intent(out) :: error
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    error%code = code
    error%message = message
  end subroutine set_error
end module tscopula_status
