! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_status
  implicit none
  private
  integer, parameter, public :: tvmvp_ok = 0
  integer, parameter, public :: tvmvp_invalid_input = 1
  integer, parameter, public :: tvmvp_singular = 2
  integer, parameter, public :: tvmvp_not_positive_definite = 3
  integer, parameter, public :: tvmvp_convergence_failure = 4
  integer, parameter, public :: tvmvp_insufficient_data = 5

  type, public :: tvmvp_error
    integer :: code = tvmvp_ok
    character(len=256) :: message = ''
  contains
    procedure :: failed => tvmvp_error_failed
  end type tvmvp_error

  public :: set_error, clear_error
contains
  pure logical function tvmvp_error_failed(self)
    class(tvmvp_error), intent(in) :: self
    tvmvp_error_failed = self%code /= tvmvp_ok
  end function tvmvp_error_failed

  subroutine clear_error(err)
    type(tvmvp_error), intent(out) :: err
    err%code = tvmvp_ok
    err%message = ''
  end subroutine clear_error

  subroutine set_error(err, code, message)
    type(tvmvp_error), intent(out) :: err
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    err%code = code
    err%message = message
  end subroutine set_error
end module tvmvp_status
