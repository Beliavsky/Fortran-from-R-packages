! Upstream license declaration: GPL (version unspecified)
module limsolve_types
    use limsolve_kinds, only: dp
    implicit none
    private

    integer, parameter, public :: LS_SUCCESS = 0
    integer, parameter, public :: LS_INFEASIBLE = 1
    integer, parameter, public :: LS_UNBOUNDED = 2
    integer, parameter, public :: LS_SINGULAR = 3
    integer, parameter, public :: LS_MAXITER = 4
    integer, parameter, public :: LS_INVALID = 5
    integer, parameter, public :: LS_NUMERICAL = 6

    type, public :: solve_result
        real(dp), allocatable :: x(:)
        real(dp) :: residual_norm = 0.0_dp
        real(dp) :: solution_norm = 0.0_dp
        logical :: is_error = .false.
        integer :: status = LS_INVALID
        integer :: numiter = 0
        real(dp), allocatable :: unconstrained_solution(:)
        real(dp), allocatable :: covariance(:,:)
        integer :: rank_eq = 0
        integer :: rank_app = 0
    contains
        procedure :: succeeded => solve_result_succeeded
    end type solve_result

    type, public :: resolution_result
        real(dp), allocatable :: row(:)
        real(dp), allocatable :: col(:)
        integer :: nsolvable = 0
    end type resolution_result

    type, public :: range_result
        real(dp), allocatable :: range(:,:) ! (n,2): min,max
        real(dp), allocatable :: central(:)
        real(dp), allocatable :: all_x(:,:)
        integer :: status = LS_INVALID
    end type range_result

    type, public :: sample_result
        real(dp), allocatable :: x(:,:)
        real(dp), allocatable :: q(:,:)
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: jump(:)
        real(dp) :: accepted_ratio = 0.0_dp
        integer :: status = LS_INVALID
    end type sample_result

contains

    logical function solve_result_succeeded(this)
        class(solve_result), intent(in) :: this
        solve_result_succeeded = this%status == LS_SUCCESS .and. .not. this%is_error
    end function solve_result_succeeded

end module limsolve_types
