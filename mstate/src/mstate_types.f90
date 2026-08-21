module mstate_types
    use mstate_kinds, only : dp
    implicit none
    private

    type, public :: transition_map
        integer :: nstate = 0
        integer :: ntrans = 0
        integer, allocatable :: trans(:, :) ! 0 means unavailable
        integer, allocatable :: from(:)
        integer, allocatable :: to(:)
    end type transition_map

    type, public :: msdata_type
        integer :: n = 0
        integer, allocatable :: id(:), from(:), to(:), trans(:), status(:)
        real(dp), allocatable :: tstart(:), tstop(:), time(:)
    end type msdata_type

    type, public :: etmdata_type
        integer :: n = 0
        integer, allocatable :: id(:), from(:), to(:)
        real(dp), allocatable :: entry(:), exit(:)
    end type etmdata_type

    type, public :: hazard_type
        integer :: nt = 0
        integer :: ntrans = 0
        real(dp), allocatable :: time(:)
        real(dp), allocatable :: haz(:, :)       ! cumulative hazard, (nt,ntrans)
        real(dp), allocatable :: varhaz(:, :, :) ! cumulative covariances, (nt,ntrans,ntrans)
    end type hazard_type

    type, public :: probtrans_type
        integer :: nt = 0
        integer :: nstate = 0
        real(dp), allocatable :: time(:)
        real(dp), allocatable :: p(:, :, :)  ! (nt,nstate,nstate)
        real(dp), allocatable :: se(:, :, :) ! same shape, optional/zero if unavailable
    end type probtrans_type

    type, public :: redrank_result
        integer :: rank = 0
        integer :: niter = 0
        integer :: df = 0
        real(dp) :: loglik = 0.0_dp
        logical :: converged = .false.
        real(dp), allocatable :: alpha(:, :)
        real(dp), allocatable :: gamma(:, :)
        real(dp), allocatable :: beta(:, :)
        real(dp), allocatable :: beta2(:)
        real(dp), allocatable :: alphax(:, :)
    end type redrank_result


    type, public :: cuminc_result
        integer :: nt = 0
        integer :: ncause = 0
        real(dp), allocatable :: time(:)
        real(dp), allocatable :: surv(:), se_surv(:)
        real(dp), allocatable :: cif(:, :), se_cif(:, :)
    end type cuminc_result

    type, public :: crprep_type
        integer :: n = 0
        integer :: nkeep = 0
        logical :: has_truncation = .false.
        integer, allocatable :: id(:), status(:), strata(:), count(:), failcode(:)
        real(dp), allocatable :: tstart(:), tstop(:), weight_cens(:), weight_trunc(:)
        real(dp), allocatable :: keep(:, :)
    end type crprep_type

    type, public :: censor_distribution
        integer :: n = 0
        real(dp), allocatable :: time(:), haz(:), surv(:)
    end type censor_distribution

    type, public :: simulated_msdata_type
        integer :: n = 0
        integer, allocatable :: id(:), from(:), to(:), status(:), trans(:)
        real(dp), allocatable :: tstart(:), tstop(:), duration(:)
    end type simulated_msdata_type


    type, public :: relative_bootstrap_type
        integer :: nt = 0
        integer :: ntrans = 0
        integer :: norig = 0
        integer :: b = 0
        integer, allocatable :: nvalid(:), original_nvalid(:)
        logical, allocatable :: valid_rep(:)
        logical :: has_original = .false.
        real(dp), allocatable :: haz(:, :, :)      ! (nt,ntrans,b), NaN when unavailable
        real(dp), allocatable :: varhaz(:, :)      ! pointwise bootstrap variances
        real(dp), allocatable :: original_haz(:, :, :)
        real(dp), allocatable :: original_varhaz(:, :)
    end type relative_bootstrap_type

    type, public :: relative_msfit_type
        character(len=9) :: variance_mode = 'fixed'
        logical :: has_bootstrap = .false.
        type(hazard_type) :: fit
        type(hazard_type) :: bootstrap_fit
        type(transition_map) :: trans
        type(relative_bootstrap_type) :: bootstrap
        integer, allocatable :: link(:, :), nlink(:)
        real(dp), allocatable :: population_haz(:, :)
    end type relative_msfit_type

    type, public :: markov_test_result
        integer :: transition = 0
        integer :: from_state = 0
        integer :: to_state = 0
        integer :: b = 0
        integer :: nsub = 0
        character(len=8) :: dist = 'poisson'
        integer, allocatable :: qualset(:), nobs_grid(:)
        real(dp), allocatable :: grid(:)
        real(dp), allocatable :: beta(:), beta_vcov(:, :)
        real(dp), allocatable :: zbar(:, :)
        real(dp), allocatable :: est_cov(:, :, :)
        real(dp), allocatable :: obs_chisq_trace(:)
        real(dp), allocatable :: n_wb_trace(:, :, :)
        real(dp), allocatable :: nch_wb_trace(:, :)
        real(dp), allocatable :: est_quant(:, :, :)
        real(dp), allocatable :: orig_stat(:), p_stat_wb(:)
        real(dp) :: orig_ch_stat = 0.0_dp
        real(dp) :: p_ch_stat_wb = 0.0_dp
    end type markov_test_result

end module mstate_types
