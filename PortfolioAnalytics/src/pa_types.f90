! SPDX-License-Identifier: GPL-3.0-only
module pa_types
  use pa_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: obj_min_variance = 1
  integer, parameter, public :: obj_max_return = 2
  integer, parameter, public :: obj_quadratic_utility = 3
  integer, parameter, public :: obj_max_sharpe = 4
  integer, parameter, public :: obj_min_es = 5
  integer, parameter, public :: obj_max_starr = 6
  integer, parameter, public :: obj_risk_parity = 7
  integer, parameter, public :: obj_min_semideviation = 8
  integer, parameter, public :: obj_min_drawdown = 9
  integer, parameter, public :: obj_min_concentration = 10
  integer, parameter, public :: obj_min_csm = 11
  integer, parameter, public :: obj_max_csm_ratio = 12
  integer, parameter, public :: obj_min_eqs = 13

  integer, parameter, public :: opt_auto = 0
  integer, parameter, public :: opt_projected_gradient = 1
  integer, parameter, public :: opt_differential_evolution = 2
  integer, parameter, public :: opt_random_search = 3

  integer, parameter, public :: pa_success = 0
  integer, parameter, public :: pa_invalid_input = 1
  integer, parameter, public :: pa_infeasible = 2
  integer, parameter, public :: pa_max_iterations = 3

  type, public :: portfolio_constraints
    real(dp), allocatable :: min_weight(:)
    real(dp), allocatable :: max_weight(:)
    real(dp) :: min_sum = 1.0_dp
    real(dp) :: max_sum = 1.0_dp
    real(dp), allocatable :: group_a(:,:)
    real(dp), allocatable :: group_lower(:)
    real(dp), allocatable :: group_upper(:)
    real(dp), allocatable :: factor_loadings(:,:)
    real(dp), allocatable :: factor_lower(:)
    real(dp), allocatable :: factor_upper(:)
    real(dp), allocatable :: initial_weights(:)
    real(dp) :: turnover_limit = -1.0_dp
    real(dp) :: diversification_min = -1.0_dp
    real(dp) :: return_target = -huge(1.0_dp)
    real(dp) :: leverage_limit = -1.0_dp
    real(dp), allocatable :: transaction_cost(:)
    real(dp) :: transaction_cost_limit = -1.0_dp
    integer :: max_positions = -1
    integer :: max_long = -1
    integer :: max_short = -1
    real(dp) :: position_tolerance = 1.0e-10_dp
  end type portfolio_constraints

  type, public :: portfolio_options
    integer :: objective = obj_min_variance
    integer :: optimizer = opt_auto
    real(dp) :: alpha = 0.05_dp
    real(dp) :: risk_aversion = 1.0_dp
    real(dp) :: risk_free = 0.0_dp
    real(dp) :: target_return = 0.0_dp
    real(dp) :: turnover_aversion = 0.0_dp
    real(dp) :: concentration_aversion = 0.0_dp
    real(dp) :: penalty_scale = 1.0e6_dp
    integer :: max_iterations = 1000
    integer :: population_size = 80
    integer :: random_portfolios = 5000
    integer :: seed = 12345
    real(dp) :: tolerance = 1.0e-8_dp
    logical :: local_refine = .true.
  end type portfolio_options

  type, public :: portfolio_result
    real(dp), allocatable :: weights(:)
    real(dp), allocatable :: risk_contribution(:)
    real(dp) :: objective_value = huge(1.0_dp)
    real(dp) :: expected_return = 0.0_dp
    real(dp) :: risk = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: expected_shortfall = 0.0_dp
    real(dp) :: turnover = 0.0_dp
    real(dp) :: hhi = 0.0_dp
    real(dp) :: diversification = 0.0_dp
    logical :: converged = .false.
    logical :: feasible = .false.
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: status = pa_invalid_input
    character(len=160) :: message = ''
  end type portfolio_result

  type, public :: frontier_result
    real(dp), allocatable :: weights(:,:)
    real(dp), allocatable :: target_return(:)
    real(dp), allocatable :: expected_return(:)
    real(dp), allocatable :: risk(:)
    logical, allocatable :: feasible(:)
  end type frontier_result

  type, public :: entropy_result
    real(dp), allocatable :: probabilities(:)
    real(dp), allocatable :: dual(:)
    logical :: converged = .false.
    integer :: iterations = 0
    real(dp) :: objective = huge(1.0_dp)
    real(dp) :: max_violation = huge(1.0_dp)
  end type entropy_result

  type, public :: factor_model_result
    real(dp), allocatable :: means(:)
    real(dp), allocatable :: loadings(:,:)
    real(dp), allocatable :: factors(:,:)
    real(dp), allocatable :: residuals(:,:)
    real(dp), allocatable :: eigenvalues(:)
    integer :: nobs = 0
    integer :: nassets = 0
    integer :: nfactors = 0
  end type factor_model_result

  type, public :: rebalancing_result
    real(dp), allocatable :: weights(:,:)
    real(dp), allocatable :: portfolio_return(:)
    real(dp), allocatable :: wealth(:)
    real(dp), allocatable :: turnover(:)
    real(dp), allocatable :: transaction_cost(:)
  end type rebalancing_result

end module pa_types
