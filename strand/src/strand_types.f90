! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand_types
  use strand_kinds, only : dp
  implicit none
  private

  type, public :: lp_result
    real(dp), allocatable :: x(:)
    real(dp) :: objective = 0.0_dp
    integer :: iterations = 0
    logical :: optimal = .false.
    logical :: feasible = .false.
    logical :: unbounded = .false.
    character(len=160) :: message = ''
  end type lp_result

  type, public :: strategy_spec
    character(len=48) :: name = ''
    real(dp) :: capital = 0.0_dp
    real(dp) :: ideal_long_weight = 1.0_dp
    real(dp) :: ideal_short_weight = 1.0_dp
    real(dp) :: target_long_weight = 0.0_dp
    real(dp) :: target_short_weight = 0.0_dp
    logical :: has_target_weights = .false.
    real(dp) :: position_limit_pct_adv = 30.0_dp
    real(dp) :: position_limit_pct_lmv = 1.0_dp
    real(dp) :: position_limit_pct_smv = 1.0_dp
    real(dp) :: trading_limit_pct_adv = 5.0_dp
  end type strategy_spec

  type, public :: factor_constraint
    character(len=64) :: name = ''
    integer :: strategy = 0
    real(dp), allocatable :: values(:)
    real(dp) :: lower_weight = -huge(1.0_dp)
    real(dp) :: upper_weight = huge(1.0_dp)
  end type factor_constraint

  type, public :: category_constraint
    character(len=64) :: name = ''
    integer :: strategy = 0
    integer, allocatable :: category(:)
    real(dp) :: lower_weight = -huge(1.0_dp)
    real(dp) :: upper_weight = huge(1.0_dp)
  end type category_constraint

  type, public :: optimizer_config
    type(strategy_spec), allocatable :: strategies(:)
    real(dp) :: turnover_limit = -1.0_dp
    character(len=16) :: target_weight_policy = 'full'
    real(dp) :: max_weight_change = -1.0_dp
    real(dp), allocatable :: loosening_sequence(:)
  end type optimizer_config

  type, public :: simulation_config
    real(dp) :: fill_rate_pct_volume = 4.0_dp
    real(dp) :: transaction_cost_pct = 0.0_dp
    real(dp) :: financing_cost_pct = 0.0_dp
    real(dp) :: force_trim_factor = -1.0_dp
    logical :: force_exit_non_investable = .false.
  end type simulation_config

  type, public :: optimization_result
    real(dp), allocatable :: trade_nmv(:, :)
    integer, allocatable :: order_shares(:, :)
    integer, allocatable :: order_shares_joint(:)
    real(dp), allocatable :: order_nmv_joint(:)
    real(dp), allocatable :: max_pos_long(:, :)
    real(dp), allocatable :: max_pos_short(:, :)
    real(dp), allocatable :: max_order_gmv(:, :)
    real(dp), allocatable :: target_long_weight(:)
    real(dp), allocatable :: target_short_weight(:)
    real(dp), allocatable :: loosened_fraction(:)
    real(dp) :: objective = 0.0_dp
    integer :: iterations = 0
    logical :: success = .false.
    character(len=160) :: message = ''
  end type optimization_result

  type, public :: exposure_result
    real(dp), allocatable :: factor(:, :)
    real(dp), allocatable :: category(:, :, :)
    integer, allocatable :: category_level(:, :)
    integer, allocatable :: category_count(:)
  end type exposure_result

  type, public :: day_result
    integer, allocatable :: start_shares(:, :)
    integer, allocatable :: order_shares(:, :)
    integer, allocatable :: market_order_shares(:, :)
    integer, allocatable :: transfer_order_shares(:, :)
    integer, allocatable :: market_fill_shares(:, :)
    integer, allocatable :: transfer_fill_shares(:, :)
    integer, allocatable :: end_shares(:, :)
    real(dp), allocatable :: start_nmv(:, :)
    real(dp), allocatable :: end_nmv(:, :)
    real(dp), allocatable :: position_pnl(:, :)
    real(dp), allocatable :: trade_costs(:, :)
    real(dp), allocatable :: financing_costs(:, :)
    real(dp), allocatable :: gross_pnl(:, :)
    real(dp), allocatable :: net_pnl(:, :)
    real(dp), allocatable :: fill_rate(:)
    logical :: success = .false.
    character(len=160) :: message = ''
  end type day_result

  type, public :: simulation_result
    type(day_result), allocatable :: day(:)
    integer, allocatable :: final_shares(:, :)
    real(dp), allocatable :: gross_pnl(:, :)
    real(dp), allocatable :: net_pnl(:, :)
    real(dp), allocatable :: end_gmv(:, :)
    real(dp), allocatable :: end_nmv(:, :)
    real(dp), allocatable :: turnover(:, :)
    logical :: success = .false.
    character(len=160) :: message = ''
  end type simulation_result

  type, public :: performance_stats
    real(dp) :: total_pnl = 0.0_dp
    real(dp) :: total_return = 0.0_dp
    real(dp) :: annualized_return = 0.0_dp
    real(dp) :: annualized_volatility = 0.0_dp
    real(dp) :: annualized_sharpe = 0.0_dp
    real(dp) :: max_drawdown = 0.0_dp
    real(dp) :: average_gmv = 0.0_dp
    real(dp) :: average_nmv = 0.0_dp
    real(dp) :: average_turnover = 0.0_dp
    real(dp) :: holding_period_months = 0.0_dp
  end type performance_stats

end module strand_types
