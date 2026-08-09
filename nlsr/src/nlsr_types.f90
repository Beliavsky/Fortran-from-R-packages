! SPDX-License-Identifier: GPL-2.0-only
module nlsr_types
    use nlsr_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: nlsr_ok = 0
    integer, parameter, public :: nlsr_max_feval = 1
    integer, parameter, public :: nlsr_max_jeval = 2
    integer, parameter, public :: nlsr_no_progress = 3
    integer, parameter, public :: nlsr_bad_input = 4
    integer, parameter, public :: nlsr_callback_error = 5
    integer, parameter, public :: nlsr_singular = 6

    integer, parameter, public :: jac_forward = 1
    integer, parameter, public :: jac_backward = 2
    integer, parameter, public :: jac_central = 3
    integer, parameter, public :: jac_richardson = 4

    type, public :: nlsr_control
        integer :: femax = 10000
        integer :: jemax = 5000
        real(dp) :: lamda = 1.0e-4_dp
        real(dp) :: laminc = 10.0_dp
        real(dp) :: lamdec = 4.0_dp
        integer :: nbtlim = 6
        real(dp) :: ndstep = 1.0e-7_dp
        real(dp) :: offset = 100.0_dp
        real(dp) :: phi = 1.0_dp
        real(dp) :: psi = 0.0_dp
        logical :: rofftest = .true.
        logical :: smallsstest = .true.
        real(dp) :: stepredn = 0.0_dp
        real(dp) :: scale_offset = 1.0_dp
        integer :: jacobian_method = jac_central
    end type nlsr_control

    type, public :: nlsr_result
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: jacobian(:,:)
        real(dp), allocatable :: lower(:)
        real(dp), allocatable :: upper(:)
        real(dp), allocatable :: weights(:)
        integer, allocatable :: bdmask(:)
        logical, allocatable :: masked(:)
        real(dp) :: ssquares = huge(1.0_dp)
        real(dp) :: roff = huge(1.0_dp)
        real(dp) :: lamda = 1.0e-4_dp
        integer :: feval = 0
        integer :: jeval = 0
        integer :: status = nlsr_bad_input
        logical :: converged = .false.
    end type nlsr_result

    abstract interface
        subroutine residual_fn(par, residual, ierr)
            import dp
            real(dp), intent(in) :: par(:)
            real(dp), intent(out) :: residual(:)
            integer, intent(out) :: ierr
        end subroutine residual_fn

        subroutine jacobian_fn(par, jac, ierr)
            import dp
            real(dp), intent(in) :: par(:)
            real(dp), intent(out) :: jac(:,:)
            integer, intent(out) :: ierr
        end subroutine jacobian_fn

        subroutine weights_fn(par, residual, weights, ierr)
            import dp
            real(dp), intent(in) :: par(:), residual(:)
            real(dp), intent(out) :: weights(:)
            integer, intent(out) :: ierr
        end subroutine weights_fn
    end interface

    public :: residual_fn, jacobian_fn, weights_fn
end module nlsr_types
