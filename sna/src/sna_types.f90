! sna-fortran: computational translation of the R sna package.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_types
    use sna_kinds, only : dp
    implicit none
    private

    type, public :: geodist_result
        real(dp), allocatable :: distance(:,:)
        real(dp), allocatable :: counts(:,:)
    end type geodist_result

    type, public :: component_result
        integer, allocatable :: membership(:)
        integer, allocatable :: csize(:)
        integer :: n_components = 0
    end type component_result

    type, public :: clique_result
        integer, allocatable :: count_by_size(:)
        integer, allocatable :: vertex_count(:,:)
        integer, allocatable :: comembership(:,:,:)
    end type clique_result

    type, public :: path_census_result
        real(dp), allocatable :: count(:)
        real(dp), allocatable :: vertex_count(:,:)
        real(dp), allocatable :: dyad_count(:,:,:)
    end type path_census_result

    type, public :: regression_result
        real(dp), allocatable :: coef(:)
        real(dp), allocatable :: se(:)
        real(dp), allocatable :: statistic(:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residual(:)
        real(dp) :: sigma2 = 0.0_dp
        real(dp) :: loglik = 0.0_dp
        integer :: nobs = 0
        integer :: rank = 0
        logical :: converged = .false.
    end type regression_result

    type, public :: qap_result
        real(dp) :: observed = 0.0_dp
        real(dp), allocatable :: simulated(:)
        real(dp) :: p_lower = 0.0_dp
        real(dp) :: p_upper = 0.0_dp
        real(dp) :: p_two_sided = 0.0_dp
    end type qap_result

    type, public :: brokerage_result
        real(dp), allocatable :: raw(:,:)
        real(dp), allocatable :: aggregate(:,:)
        real(dp), allocatable :: expected(:,:), sd(:,:), z(:,:)
        real(dp), allocatable :: global_expected(:), global_sd(:), global_z(:)
        real(dp), allocatable :: group_expected(:,:), group_sd(:,:)
        integer, allocatable :: class_ids(:), class_sizes(:)
    end type brokerage_result

end module sna_types
