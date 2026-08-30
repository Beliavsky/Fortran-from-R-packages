! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_types
    use classint_kinds, only: dp
    use e1071, only: bclust_model
    implicit none
    private

    type, public :: classint_options
        real(dp), allocatable :: fixed_breaks(:)
        real(dp), allocatable :: sd_m(:)
        integer :: quantile_type = 7
        integer :: box_quantile_type = 7
        real(dp) :: headtails_threshold = 0.4_dp
        real(dp) :: box_iqr_mult = 1.5_dp
        logical :: box_legacy = .false.
        character(len=16) :: hclust_method = "complete"
        integer :: kmeans_iter_max = 100
        integer :: kmeans_nstart = 1
        integer :: seed = 1
        logical :: sample_large_fisher_jenks = .true.
        integer :: large_n = 3000
        real(dp) :: sample_proportion = 0.1_dp
        integer :: dpih_level = 2
        integer :: dpih_gridsize = 401
        character(len=8) :: dpih_scale = "minim"
        logical :: dpih_has_range = .false.
        real(dp) :: dpih_range(2) = 0.0_dp
        logical :: dpih_truncate = .true.
        integer :: bclust_iter_base = 10
        integer :: bclust_minsize = 0
        character(len=16) :: bclust_distance = "euclidean"
        character(len=16) :: bclust_hclust = "average"
        character(len=16) :: bclust_base_method = "kmeans"
        integer :: bclust_base_centers = 20
        logical :: bclust_final_kmeans = .false.
        logical :: bclust_resample = .true.
        integer :: bclust_maxcluster = 20
    end type classint_options

    type, public :: hclust_1d_model
        integer :: n = 0
        integer, allocatable :: parent(:)
        real(dp), allocatable :: node_mean(:)
        character(len=16) :: method = "complete"
    end type hclust_1d_model

    type, public :: class_intervals
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: breaks(:)
        character(len=16) :: style = ""
        character(len=5) :: interval_closure = "left"
        integer :: nobs = 0
        logical :: sampled = .false.
        real(dp), allocatable :: fisher_stats(:, :)
        logical :: has_hclust = .false.
        type(hclust_1d_model) :: hclust
        logical :: has_bclust = .false.
        type(bclust_model) :: bclust
    end type class_intervals

    type, public :: jenks_indices
        integer :: n_classes = 0
        real(dp) :: goodness_of_fit = 0.0_dp
        real(dp) :: tabular_accuracy = 0.0_dp
        real(dp) :: overview_accuracy = 0.0_dp
        logical :: has_overview = .false.
    end type jenks_indices
end module classint_types
