! SPDX-License-Identifier: MIT
! Translated from the computational core of the R package trust 0.1-9.
module trust_types
    use trust_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: trust_step_newton = 1
    integer, parameter, public :: trust_step_easy_easy = 2
    integer, parameter, public :: trust_step_hard_easy = 3
    integer, parameter, public :: trust_step_hard_hard = 4

    integer, parameter, public :: trust_ok = 0
    integer, parameter, public :: trust_err_bad_input = 1
    integer, parameter, public :: trust_err_initial_objective = 2
    integer, parameter, public :: trust_err_objective = 3
    integer, parameter, public :: trust_err_eigensolver = 4

    type, public :: trust_options
        real(dp) :: rinit = 1.0_dp
        real(dp) :: rmax = 100.0_dp
        integer :: iterlim = 100
        real(dp) :: fterm = sqrt(epsilon(1.0_dp))
        real(dp) :: mterm = sqrt(epsilon(1.0_dp))
        logical :: minimize = .true.
        logical :: save_history = .false.
        real(dp), allocatable :: parscale(:)
    end type trust_options

    type, public :: trust_history
        integer :: n = 0
        real(dp), allocatable :: argument(:, :)
        real(dp), allocatable :: argument_try(:, :)
        integer, allocatable :: step_type(:)
        logical, allocatable :: accepted(:)
        real(dp), allocatable :: radius(:)
        real(dp), allocatable :: step_norm(:)
        real(dp), allocatable :: rho(:)
        real(dp), allocatable :: value(:)
        real(dp), allocatable :: value_try(:)
        real(dp), allocatable :: predicted_difference(:)
    end type trust_history

    type, public :: trust_result
        real(dp), allocatable :: argument(:)
        real(dp), allocatable :: gradient(:)
        real(dp), allocatable :: hessian(:, :)
        real(dp) :: value = huge(1.0_dp)
        logical :: converged = .false.
        integer :: iterations = 0
        integer :: status = trust_ok
        character(len=160) :: message = ''
        type(trust_history) :: history
    end type trust_result

    abstract interface
        subroutine trust_objective(x, value, gradient, hessian, status)
            import :: dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: value
            real(dp), intent(out) :: gradient(:)
            real(dp), intent(out) :: hessian(:, :)
            integer, intent(out) :: status
        end subroutine trust_objective
    end interface
    public :: trust_objective

end module trust_types
