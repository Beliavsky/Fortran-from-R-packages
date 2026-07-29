! SPDX-License-Identifier: MIT
module bayesianou_types
  use bayesianou_kinds, only : dp, status_ok
  implicit none
  private

  integer, parameter, public :: level_canonical = 1
  integer, parameter, public :: level_both_full = 2
  integer, parameter, public :: level_both_lean = 3
  integer, parameter, public :: level_n1_lean = 4

  type, public :: ou_level_flags
    logical :: cubic = .true.
    logical :: sv = .true.
    logical :: student_t = .true.
    logical :: hierarchy = .true.
  end type ou_level_flags

  type, public :: ou_level_spec_type
    type(ou_level_flags) :: level1
    type(ou_level_flags) :: level2
  end type ou_level_spec_type

  type, public :: zscore_result
    real(dp), allocatable :: mz(:,:)
    real(dp), allocatable :: mu(:)
    real(dp), allocatable :: sd(:)
    integer :: status = status_ok
  end type zscore_result

  type, public :: ou_priors
    real(dp) :: sigma_delta = 0.002_dp
    real(dp) :: beta1_mean = 0.0_dp
    real(dp) :: beta1_sd = 0.5_dp
    real(dp) :: nu_shape = 2.0_dp
    real(dp) :: nu_rate = 0.1_dp
    real(dp) :: rho_mean = 0.7_dp
    real(dp) :: rho_sd = 0.2_dp
    real(dp) :: sigma_phi_meas_sd = 0.5_dp
    real(dp) :: sigma_k_recon = 0.10_dp
  end type ou_priors

  type, public :: ou_options
    integer :: n_levels = 1
    real(dp) :: train_frac = 0.70_dp
    logical :: fit_full = .false.
    integer :: chains = 4
    integer :: iterations = 2000
    integer :: warmup = 1000
    integer :: thin = 1
    integer :: seed = 1234
    logical :: com_in_mean = .true.
    logical :: hard_sum_zero = .true.
    logical :: orthogonalize_tmg = .true.
    logical :: factor_from_x = .true.
    logical :: use_train_loadings = .true.
    logical :: theta_separation_hard = .false.
    logical :: k_reconstruction = .false.
    real(dp) :: sigma_phi_meas_fixed = -1.0_dp
    real(dp) :: kappa_cap = 2.0_dp
    real(dp) :: proposal_scale = 0.045_dp
    type(ou_level_spec_type) :: level_spec = ou_level_spec_type( &
      ou_level_flags(.true.,.true.,.true.,.true.), &
      ou_level_flags(.false.,.false.,.false.,.true.))
    type(ou_priors) :: priors
  end type ou_options

  type, public :: ou_input
    real(dp), allocatable :: y(:,:), x(:,:), tmg(:), com(:,:), capital(:,:)
    real(dp), allocatable :: gprime(:), value_anchor(:,:)
    real(dp), allocatable :: k_cost(:,:), k_hat(:,:)
  end type ou_input

  type, public :: ou_summary
    real(dp), allocatable :: theta(:), kappa(:), a3(:), beta0(:)
    real(dp), allocatable :: alpha(:), rho(:), sigma_eta(:)
    real(dp), allocatable :: kappa_p(:), mu_const(:), sigma_p(:), a3_p(:)
    real(dp) :: beta1 = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: nu = 10.0_dp
    real(dp) :: m1 = 0.0_dp
    real(dp) :: m_v = 0.0_dp
    real(dp) :: sigma_phi_meas = 0.5_dp
  end type ou_summary

  type, public :: ou_draws
    real(dp), allocatable :: theta(:,:), kappa(:,:), a3(:,:), beta0(:,:)
    real(dp), allocatable :: alpha(:,:), rho(:,:), sigma_eta(:,:)
    real(dp), allocatable :: kappa_p(:,:), mu_const(:,:), sigma_p(:,:), a3_p(:,:)
    real(dp), allocatable :: beta1(:), gamma(:), nu(:), m1(:), m_v(:), sigma_phi_meas(:)
    real(dp), allocatable :: phi(:,:,:), h(:,:,:)
    integer, allocatable :: chain(:)
  end type ou_draws

  type, public :: loo_result
    real(dp) :: elpd_loo = 0.0_dp
    real(dp) :: se_elpd = 0.0_dp
    real(dp) :: p_loo = 0.0_dp
    real(dp) :: looic = 0.0_dp
    real(dp), allocatable :: pointwise(:)
    real(dp), allocatable :: pareto_k(:)
  end type loo_result

  type, public :: oos_metric
    integer :: horizon = 1
    real(dp) :: rmse = 0.0_dp
    real(dp) :: mae = 0.0_dp
    integer :: n_obs = 0
  end type oos_metric

  type, public :: ou_diagnostics
    real(dp), allocatable :: rhat(:), ess(:)
    real(dp) :: rhat_max = 1.0_dp
    real(dp) :: rhat_share = 0.0_dp
    integer :: divergences = 0
    real(dp) :: acceptance_rate = 0.0_dp
    type(loo_result) :: loo
    type(oos_metric), allocatable :: oos(:)
  end type ou_diagnostics

  type, public :: ou_fit_result
    type(ou_options) :: options
    type(ou_summary) :: summary
    type(ou_draws) :: draws
    type(ou_diagnostics) :: diagnostics
    type(zscore_result) :: zy, zx
    real(dp), allocatable :: ztmg(:), common_factor(:), com_wmean(:), com_wsd(:)
    real(dp), allocatable :: phi_median(:,:), h_median(:,:)
    integer :: t_train = 0
    integer :: t_lik = 0
    integer :: status = status_ok
    character(len=256) :: message = ''
  end type ou_fit_result

  type, public :: stability_result
    real(dp), allocatable :: intervals(:,:)
    logical :: stable = .false.
    real(dp) :: probability = 0.0_dp
  end type stability_result

  type, public :: rubin_result
    real(dp), allocatable :: estimate(:), total_sd(:), within_var(:), between_var(:)
    real(dp), allocatable :: df(:), fmi(:), lo(:), hi(:)
  end type rubin_result

  type, public :: ou_mi_result
    type(ou_fit_result), allocatable :: fits(:)
    type(rubin_result) :: pooled
    integer :: status = status_ok
  end type ou_mi_result

end module bayesianou_types
