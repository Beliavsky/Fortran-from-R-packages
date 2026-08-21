module circstats_types
    use circstats_kinds, only: dp
    implicit none
    private

    type, public :: circ_dispersion_result
        integer :: n = 0
        real(dp) :: r = 0.0_dp
        real(dp) :: rbar = 0.0_dp
        real(dp) :: variance = 0.0_dp
    end type

    type, public :: circ_summary_result
        integer :: n = 0
        real(dp) :: mean_dir = 0.0_dp
        real(dp) :: rho = 0.0_dp
    end type

    type, public :: trig_moment_result
        real(dp) :: mu = 0.0_dp
        real(dp) :: rho = 0.0_dp
        real(dp) :: cosine = 0.0_dp
        real(dp) :: sine = 0.0_dp
    end type

    type, public :: test_result
        real(dp) :: statistic = 0.0_dp
        real(dp) :: p_value = -1.0_dp
        real(dp) :: critical = -1.0_dp
        real(dp) :: p_lower = -1.0_dp
        real(dp) :: p_upper = -1.0_dp
        integer :: df = 0
        logical :: reject = .false.
    end type

    type, public :: circ_cor_result
        real(dp) :: r = 0.0_dp
        real(dp) :: statistic = 0.0_dp
        real(dp) :: p_value = -1.0_dp
    end type

    type, public :: vm_fit_result
        real(dp) :: mu = 0.0_dp
        real(dp) :: kappa = 0.0_dp
    end type

    type, public :: wrapped_cauchy_fit_result
        real(dp) :: mu = 0.0_dp
        real(dp) :: rho = 0.0_dp
        integer :: iterations = 0
        logical :: converged = .false.
    end type

    type, public :: change_point_result
        integer :: n = 0
        real(dp) :: rho = 0.0_dp
        real(dp) :: rmax = 0.0_dp
        real(dp) :: rave = 0.0_dp
        real(dp) :: tmax = 0.0_dp
        real(dp) :: tave = 0.0_dp
        integer :: k_r = 0
        integer :: k_t = 0
    end type

    type, public :: rao_homogeneity_result
        type(test_result) :: polar
        type(test_result) :: dispersion
    end type

    type, public :: circ_reg_result
        real(dp) :: rho = 0.0_dp
        real(dp) :: a_k = 0.0_dp
        real(dp) :: kappa = 0.0_dp
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        real(dp), allocatable :: coef(:,:)
        real(dp) :: pvalues(2) = 1.0_dp
        logical :: higher_order_significant = .false.
    end type

    type, public :: vm_bootstrap_result
        real(dp) :: mu_ci(2) = 0.0_dp
        real(dp) :: kappa_ci(2) = 0.0_dp
        real(dp), allocatable :: mu_reps(:)
        real(dp), allocatable :: kappa_reps(:)
    end type
end module circstats_types
