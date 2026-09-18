module tclust_types
    use tclust_kinds, only : dp
    implicit none
    private

    type, public :: tclust_result
        integer :: n = 0
        integer :: p = 0
        integer :: k = 0
        integer :: code = 0
        real(dp) :: alpha = 0.0_dp
        real(dp) :: obj = -huge(1.0_dp)
        real(dp) :: nlogl = huge(1.0_dp)
        real(dp) :: restr_fact = 12.0_dp
        real(dp) :: cshape = 1.0e10_dp
        real(dp) :: unrestr_fact = 1.0_dp
        real(dp) :: clacla = huge(1.0_dp)
        real(dp) :: mixmix = huge(1.0_dp)
        real(dp) :: mixcla = huge(1.0_dp)
        character(len=8) :: restriction = 'eigen'
        character(len=4) :: opt = 'HARD'
        logical :: equal_weights = .false.
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: size(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: cov(:, :, :)
        real(dp), allocatable :: posterior(:, :)
        real(dp), allocatable :: mah(:)
        real(dp), allocatable :: x(:, :)
    end type tclust_result

    type, public :: tkmeans_result
        integer :: n = 0
        integer :: p = 0
        integer :: k = 0
        integer :: code = 0
        real(dp) :: alpha = 0.0_dp
        real(dp) :: obj = huge(1.0_dp)
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: size(:)
        real(dp), allocatable :: weights(:)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: x(:, :)
    end type tkmeans_result

    type, public :: rlg_result
        integer :: n = 0
        integer :: p = 0
        integer :: k = 0
        real(dp) :: alpha = 0.0_dp
        real(dp) :: obj = huge(1.0_dp)
        integer, allocatable :: dimensions(:)
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: centers(:, :)
        real(dp), allocatable :: basis(:, :, :)
        real(dp), allocatable :: x(:, :)
    end type rlg_result

    type, public :: simulation_result
        real(dp), allocatable :: x(:, :)
        integer, allocatable :: true_cluster(:)
    end type simulation_result

    type, public :: rand_index_result
        real(dp) :: adjusted_rand = 0.0_dp
        real(dp) :: rand = 0.0_dp
        real(dp) :: mirkin = 0.0_dp
        real(dp) :: hubert = 0.0_dp
    end type rand_index_result

    type, public :: fm_index_result
        real(dp) :: adjusted = 0.0_dp
        real(dp) :: raw = 0.0_dp
        real(dp) :: expected = 0.0_dp
        real(dp) :: variance = 0.0_dp
    end type fm_index_result

    type, public :: discr_fact_result
        real(dp) :: threshold = 0.0_dp
        real(dp) :: ylimmin = 0.0_dp
        integer, allocatable :: assignment(:)
        integer, allocatable :: second_assignment(:)
        real(dp), allocatable :: likelihood(:)
        real(dp), allocatable :: second_likelihood(:)
        real(dp), allocatable :: factor(:)
        real(dp), allocatable :: mean_factor(:)
    end type discr_fact_result

    type, public :: ctlcurves_result
        integer, allocatable :: k_values(:)
        real(dp), allocatable :: alpha_values(:)
        real(dp), allocatable :: obj(:, :)
        real(dp), allocatable :: min_weights(:, :)
        real(dp), allocatable :: unrestr_fact(:, :)
        integer, allocatable :: doubtful(:, :)
        integer, allocatable :: outlying(:, :)
    end type ctlcurves_result

    type, public :: tclust_ic_result
        integer :: n = 0
        integer, allocatable :: k_values(:)
        real(dp), allocatable :: c_values(:)
        real(dp), allocatable :: clacla(:, :)
        real(dp), allocatable :: mixmix(:, :)
        real(dp), allocatable :: mixcla(:, :)
        integer, allocatable :: idxcla(:, :, :)
        integer, allocatable :: idxmix(:, :, :)
    end type tclust_ic_result

    type, public :: tclust_ic_solution_result
        integer :: nsolutions = 0
        integer, allocatable :: best_k(:)
        real(dp), allocatable :: best_c(:)
        logical, allocatable :: spurious(:)
        real(dp), allocatable :: ari_best(:, :)
        real(dp), allocatable :: ari_mix(:, :)
        real(dp), allocatable :: ari_cla(:, :)
    end type tclust_ic_solution_result

end module tclust_types
