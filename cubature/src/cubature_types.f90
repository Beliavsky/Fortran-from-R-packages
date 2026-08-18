! SPDX-License-Identifier: GPL-3.0-or-later
module cubature_types
    use cubature_kinds, only : dp, i8
    implicit none
    private

    integer, parameter, public :: ERROR_INDIVIDUAL = 0
    integer, parameter, public :: ERROR_PAIRED     = 1
    integer, parameter, public :: ERROR_L2         = 2
    integer, parameter, public :: ERROR_L1         = 3
    integer, parameter, public :: ERROR_LINF       = 4

    integer, parameter, public :: CUBATURE_SUCCESS = 0
    integer, parameter, public :: CUBATURE_MAXEVAL = 1
    integer, parameter, public :: CUBATURE_BADARG  = 2
    integer, parameter, public :: CUBATURE_FAILURE = 3

    type, public :: cubature_result
        real(dp), allocatable :: integral(:)
        real(dp), allocatable :: error(:)
        real(dp), allocatable :: prob(:)
        integer(i8) :: evaluations = 0_i8
        integer :: return_code = CUBATURE_SUCCESS
        integer :: nregions = 0
    contains
        procedure :: converged => result_converged
    end type cubature_result

    type, public :: vegas_options
        integer(i8) :: min_eval = 0_i8
        integer(i8) :: max_eval = 1000000_i8
        integer :: seed = 12345
        integer :: nstart = 1000
        integer :: nincrease = 500
        integer :: nbatch = 1000
        integer :: grid_no = 0
        integer :: nbins = 32
        integer :: max_iter = 64
    end type vegas_options

    type, public :: suave_options
        integer(i8) :: min_eval = 0_i8
        integer(i8) :: max_eval = 1000000_i8
        integer :: seed = 0
        integer :: nnew = 1000
        integer :: nmin = 50
        real(dp) :: flatness = 50.0_dp
        integer :: max_regions = 4096
    end type suave_options

    type, public :: divonne_options
        integer(i8) :: min_eval = 0_i8
        integer(i8) :: max_eval = 1000000_i8
        integer :: seed = 0
        integer :: key1 = 47
        integer :: key2 = 1
        integer :: key3 = 1
        integer :: max_pass = 5
        real(dp) :: border = 0.0_dp
        real(dp) :: max_chisq = 10.0_dp
        real(dp) :: min_deviation = 0.25_dp
        integer :: max_regions = 4096
    end type divonne_options

    type, public :: cuhre_options
        integer(i8) :: min_eval = 0_i8
        integer(i8) :: max_eval = 1000000_i8
        integer :: key = 0
    end type cuhre_options

    abstract interface
        subroutine cubature_integrand(x, value)
            import dp
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: value(:)
        end subroutine cubature_integrand

        subroutine cubature_integrand_v(x, value)
            import dp
            real(dp), intent(in) :: x(:, :)
            real(dp), intent(out) :: value(:, :)
        end subroutine cubature_integrand_v
    end interface
    public :: cubature_integrand, cubature_integrand_v, errors_converged

contains

    logical function result_converged(self, abs_tol, rel_tol, norm)
        class(cubature_result), intent(in) :: self
        real(dp), intent(in) :: abs_tol, rel_tol
        integer, intent(in), optional :: norm
        integer :: nrm
        nrm = ERROR_INDIVIDUAL
        if (present(norm)) nrm = norm
        if (.not. allocated(self%integral) .or. .not. allocated(self%error)) then
            result_converged = .false.
        else
            result_converged = errors_converged(self%integral, self%error, abs_tol, rel_tol, nrm)
        end if
    end function result_converged

    logical function errors_converged(val, err, abs_tol, rel_tol, norm)
        real(dp), intent(in) :: val(:), err(:), abs_tol, rel_tol
        integer, intent(in) :: norm
        real(dp) :: nv, ne
        integer :: j

        if (size(val) /= size(err)) then
            errors_converged = .false.
            return
        end if
        if (size(val) == 0) then
            errors_converged = .true.
            return
        end if

        select case (norm)
        case (ERROR_INDIVIDUAL)
            errors_converged = all(err <= max(abs_tol, rel_tol * abs(val)))
        case (ERROR_PAIRED)
            errors_converged = .true.
            j = 1
            do while (j <= size(val))
                if (j == size(val)) then
                    nv = abs(val(j)); ne = err(j)
                else
                    nv = sqrt(val(j) * val(j) + val(j + 1) * val(j + 1))
                    ne = sqrt(err(j) * err(j) + err(j + 1) * err(j + 1))
                end if
                if (ne > max(abs_tol, rel_tol * nv)) then
                    errors_converged = .false.
                    exit
                end if
                j = j + 2
            end do
        case (ERROR_L2)
            nv = sqrt(sum(val * val)); ne = sqrt(sum(err * err))
            errors_converged = ne <= max(abs_tol, rel_tol * nv)
        case (ERROR_L1)
            nv = sum(abs(val)); ne = sum(abs(err))
            errors_converged = ne <= max(abs_tol, rel_tol * nv)
        case (ERROR_LINF)
            nv = maxval(abs(val)); ne = maxval(abs(err))
            errors_converged = ne <= max(abs_tol, rel_tol * nv)
        case default
            errors_converged = .false.
        end select
    end function errors_converged

end module cubature_types
