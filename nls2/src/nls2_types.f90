! SPDX-License-Identifier: GPL-2.0-only
module nls2_types
    use nls2_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: nls2_ok = 0
    integer, parameter, public :: nls2_maxiter = 1
    integer, parameter, public :: nls2_singular = 2
    integer, parameter, public :: nls2_bad_input = 3
    integer, parameter, public :: nls2_model_error = 4
    integer, parameter, public :: nls2_no_finite_start = 5

    type, public :: nls_control
        integer :: maxiter = 50
        real(dp) :: tol = 1.0e-5_dp
        real(dp) :: min_factor = 1.0_dp / 1024.0_dp
        real(dp) :: scale_offset = 0.0_dp
        logical :: central_diff = .false.
        logical :: warn_only = .false.
        real(dp) :: diff_step = sqrt(epsilon(1.0_dp))
    end type nls_control

    type, public :: nls_result
        real(dp), allocatable :: par(:)
        real(dp), allocatable :: linear_par(:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: covariance(:,:)
        real(dp) :: rss = huge(1.0_dp)
        real(dp) :: sigma = huge(1.0_dp)
        real(dp) :: fin_tol = huge(1.0_dp)
        integer :: iterations = 0
        integer :: evaluations = 0
        integer :: start_index = 0
        integer :: status = nls2_bad_input
        logical :: converged = .false.
    end type nls_result

    type, public :: nls2_search_result
        type(nls_result) :: best
        type(nls_result), allocatable :: fits(:)
        real(dp), allocatable :: starts(:,:)
        real(dp), allocatable :: start_rss(:)
        integer :: n_candidates = 0
        integer :: status = nls2_bad_input
    end type nls2_search_result

    abstract interface
        subroutine nls_model(x, par, yhat, ierr)
            import dp
            real(dp), intent(in) :: x(:,:)
            real(dp), intent(in) :: par(:)
            real(dp), intent(out) :: yhat(:)
            integer, intent(out) :: ierr
        end subroutine nls_model

        subroutine nls_jacobian(x, par, jac, ierr)
            import dp
            real(dp), intent(in) :: x(:,:)
            real(dp), intent(in) :: par(:)
            real(dp), intent(out) :: jac(:,:)
            integer, intent(out) :: ierr
        end subroutine nls_jacobian

        subroutine plinear_basis(x, theta, basis, ierr)
            import dp
            real(dp), intent(in) :: x(:,:)
            real(dp), intent(in) :: theta(:)
            real(dp), intent(out) :: basis(:,:)
            integer, intent(out) :: ierr
        end subroutine plinear_basis
    end interface

    public :: nls_model, nls_jacobian, plinear_basis

end module nls2_types
