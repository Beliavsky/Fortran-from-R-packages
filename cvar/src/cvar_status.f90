! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of cvar 0.6 by Georgi N. Boshnakov.
module cvar_status
    implicit none
    private

    integer, parameter, public :: cvar_ok = 0
    integer, parameter, public :: cvar_invalid_probability = 1
    integer, parameter, public :: cvar_invalid_scale = 2
    integer, parameter, public :: cvar_invalid_sample = 3
    integer, parameter, public :: cvar_bracket_failure = 4
    integer, parameter, public :: cvar_nonconvergence = 5
    integer, parameter, public :: cvar_invalid_model = 6
    integer, parameter, public :: cvar_allocation_failure = 7

    public :: cvar_status_message

contains

    pure function cvar_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=:), allocatable :: message

        select case (status)
        case (cvar_ok)
            message = "success"
        case (cvar_invalid_probability)
            message = "probability must be strictly between zero and one"
        case (cvar_invalid_scale)
            message = "scale or slope must be positive"
        case (cvar_invalid_sample)
            message = "sample must contain at least one finite value"
        case (cvar_bracket_failure)
            message = "unable to bracket the requested quantile"
        case (cvar_nonconvergence)
            message = "numerical method did not converge"
        case (cvar_invalid_model)
            message = "invalid GARCH(1,1) model"
        case (cvar_allocation_failure)
            message = "memory allocation failed"
        case default
            message = "unknown status"
        end select
    end function cvar_status_message

end module cvar_status
