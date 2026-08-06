! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_status
    implicit none
    private

    integer, parameter, public :: ltsa_success = 0
    integer, parameter, public :: ltsa_invalid_input = 1
    integer, parameter, public :: ltsa_not_positive_definite = 2
    integer, parameter, public :: ltsa_singular = 3
    integer, parameter, public :: ltsa_nonstationary = 4
    integer, parameter, public :: ltsa_dh_condition_failed = 5
    public :: set_error

    type, public :: ltsa_error
        integer :: code = ltsa_success
        character(len=160) :: message = ''
    contains
        procedure :: ok => ltsa_error_ok
    end type ltsa_error

contains

    logical function ltsa_error_ok(self) result(ok)
        class(ltsa_error), intent(in) :: self
        ok = self%code == ltsa_success
    end function ltsa_error_ok

    subroutine set_error(error, code, message)
        type(ltsa_error), intent(out) :: error
        integer, intent(in) :: code
        character(len=*), intent(in) :: message
        error%code = code
        error%message = message
    end subroutine set_error

end module ltsa_status
