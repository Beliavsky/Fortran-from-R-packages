! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

module dirmult_types
    use iso_fortran_env, only : real64
    implicit none
    private

    integer, parameter, public :: dp = real64

    type, public :: dirmult_fit_type
        real(dp) :: loglik = 0.0_dp
        integer :: iterations = 0
        real(dp), allocatable :: gamma(:)
        real(dp), allocatable :: pi(:)
        real(dp) :: theta = 0.0_dp
        logical :: converged = .false.
        integer :: info = 0
    end type dirmult_fit_type

    type, public :: mom_result_type
        real(dp) :: theta = 0.0_dp
        real(dp) :: se = 0.0_dp
        integer :: info = 0
    end type mom_result_type

    type, public :: dirmult_summary_type
        real(dp), allocatable :: mle(:)
        real(dp), allocatable :: se_mle(:)
        real(dp), allocatable :: mom(:)
        real(dp), allocatable :: se_mom(:)
        integer :: info = 0
    end type dirmult_summary_type

    type, public :: profile_fit_type
        real(dp) :: loglik = 0.0_dp
        integer :: iterations = 0
        real(dp), allocatable :: gamma(:)
        real(dp), allocatable :: pi(:)
        real(dp) :: theta = 0.0_dp
        real(dp) :: lambda = 0.0_dp
        logical :: converged = .false.
        integer :: info = 0
    end type profile_fit_type

    type, public :: profile_grid_type
        real(dp), allocatable :: theta(:)
        real(dp), allocatable :: loglik(:)
        logical, allocatable :: success(:)
        integer :: info = 0
    end type profile_grid_type

    type, public :: count_table_type
        integer, allocatable :: x(:,:)
    end type count_table_type

    type, public :: real_vector_type
        real(dp), allocatable :: value(:)
    end type real_vector_type

    type, public :: equal_theta_table_type
        real(dp), allocatable :: gamma(:)
        real(dp), allocatable :: pi(:)
        real(dp) :: theta = 0.0_dp
        real(dp) :: lambda = 0.0_dp
    end type equal_theta_table_type

    type, public :: equal_theta_fit_type
        real(dp) :: loglik = 0.0_dp
        integer :: iterations = 0
        type(equal_theta_table_type), allocatable :: table(:)
        real(dp) :: common_gamma_sum = 0.0_dp
        logical :: converged = .false.
        integer :: info = 0
    end type equal_theta_fit_type

    type, public :: sim_pop_result_type
        real(dp) :: theta = 0.0_dp
        real(dp), allocatable :: pi(:)
        integer, allocatable :: data(:,:)
        integer :: info = 0
    end type sim_pop_result_type

    type, public :: null_test_result_type
        integer, allocatable :: simulated(:,:,:)
        real(dp), allocatable :: mle_theta(:)
        real(dp), allocatable :: dm_loglik(:)
        real(dp), allocatable :: mom(:)
        real(dp), allocatable :: mn_loglik(:)
        logical, allocatable :: converged(:)
        integer :: info = 0
    end type null_test_result_type

end module dirmult_types
