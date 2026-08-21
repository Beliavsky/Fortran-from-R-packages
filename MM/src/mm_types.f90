! Computational translation of the R package MM 1.7-0.
! Upstream license: GPL-2. This translation is GPL-2.
module mm_types
    use mm_kinds, only : dp
    implicit none
    private

    type, public :: paras_type
        real(dp), allocatable :: vals(:)
        character(len=64), allocatable :: pnames(:)
    end type paras_type

    type, public :: mb_type
        integer, allocatable :: counts(:,:)
        integer, allocatable :: m(:)
        character(len=64), allocatable :: pnames(:)
    end type mb_type

    type, public :: suffstats_type
        integer :: y_total = 0
        real(dp) :: nobs = 0.0_dp
        real(dp), allocatable :: row_sums(:)
        real(dp), allocatable :: cross_prods(:,:)
    end type suffstats_type

    type, public :: gunter_type
        integer, allocatable :: tbl(:,:)
        integer, allocatable :: d(:)
    end type gunter_type

    type, public :: gunter_mb_type
        integer, allocatable :: tbl(:,:)
        integer, allocatable :: d(:)
        integer, allocatable :: m(:)
    end type gunter_mb_type

    type, public :: glm_fit_type
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: fitted(:)
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: iterations = 0
        logical :: converged = .false.
    end type glm_fit_type

    type, public :: mm_fit_type
        type(paras_type) :: parameters
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: iterations = 0
        logical :: converged = .false.
    end type mm_fit_type

end module mm_types
