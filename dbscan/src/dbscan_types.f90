module dbscan_types
    use dbscan_kinds, only : dp
    implicit none
    private

    type, public :: knn_result
        integer :: k = 0
        integer, allocatable :: id(:, :)
        real(dp), allocatable :: dist(:, :)
        integer, allocatable :: shared(:, :)
    end type knn_result

    type, public :: frnn_result
        real(dp) :: eps = 0.0_dp
        integer :: n = 0
        integer, allocatable :: offset(:)
        integer, allocatable :: id(:)
        real(dp), allocatable :: dist(:)
    end type frnn_result

    type, public :: clustering_result
        integer, allocatable :: cluster(:)
        integer :: min_pts = 0
        integer :: k = 0
        integer :: kt = 0
        real(dp) :: eps = 0.0_dp
    end type clustering_result

    type, public :: optics_result
        integer :: min_pts = 0
        real(dp) :: eps = huge(1.0_dp)
        real(dp) :: eps_cl = -1.0_dp
        real(dp) :: xi = -1.0_dp
        integer, allocatable :: order(:)
        integer, allocatable :: predecessor(:)
        integer, allocatable :: cluster(:)
        integer, allocatable :: xi_start(:)
        integer, allocatable :: xi_end(:)
        real(dp), allocatable :: reachdist(:)
        real(dp), allocatable :: coredist(:)
    end type optics_result

    type, public :: dendrogram_result
        integer, allocatable :: merge(:, :)
        integer, allocatable :: order(:)
        real(dp), allocatable :: height(:)
    end type dendrogram_result

    type, public :: hdbscan_result
        integer :: min_pts = 0
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: membership_prob(:)
        real(dp), allocatable :: outlier_scores(:)
        real(dp), allocatable :: cluster_scores(:)
        real(dp), allocatable :: core_dist(:)
        real(dp), allocatable :: mst(:, :)
        type(dendrogram_result) :: hierarchy
    end type hdbscan_result

    type, public :: dbcv_result
        real(dp) :: score = -1.0_dp
        integer :: n = 0
        integer :: d = 0
        integer, allocatable :: cluster_size(:)
        real(dp), allocatable :: dsc(:)
        real(dp), allocatable :: dspc(:, :)
        real(dp), allocatable :: validity(:)
    end type dbcv_result
end module dbscan_types
