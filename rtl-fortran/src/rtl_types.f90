! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_types
  use rtl_kinds, only: dp
  implicit none
  private

  type, public :: status_t
    logical :: ok = .true.
    character(len=256) :: message = "ok"
  end type status_t

  type, public :: option_result
    real(dp) :: price = 0.0_dp
    real(dp) :: delta = 0.0_dp
    real(dp) :: gamma = 0.0_dp
    real(dp) :: vega = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: rho = 0.0_dp
    type(status_t) :: status
  end type option_result

  type, public :: spread_option_result
    real(dp) :: price = 0.0_dp
    real(dp) :: delta_f1 = 0.0_dp
    real(dp) :: delta_f2 = 0.0_dp
    real(dp) :: gamma_f1 = 0.0_dp
    real(dp) :: gamma_f2 = 0.0_dp
    real(dp) :: gamma_cross = 0.0_dp
    real(dp) :: vega_1 = 0.0_dp
    real(dp) :: vega_2 = 0.0_dp
    real(dp) :: theta = 0.0_dp
    real(dp) :: rho = 0.0_dp
    logical :: monitoring_used = .false.
    type(status_t) :: status
  end type spread_option_result

  type, public :: crr_tree_result
    real(dp) :: price = 0.0_dp
    real(dp), allocatable :: asset(:, :)
    real(dp), allocatable :: option(:, :)
    type(status_t) :: status
  end type crr_tree_result

  type, public :: path_result
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: values(:, :)
    type(status_t) :: status
  end type path_result

  type, public :: ou_fit_result
    real(dp) :: theta = 0.0_dp
    real(dp) :: mu = 0.0_dp
    real(dp) :: sigma = 0.0_dp
    real(dp) :: half_life_periods = 0.0_dp
    type(status_t) :: status
  end type ou_fit_result

  type, public :: multivariate_result
    real(dp), allocatable :: mean(:)
    real(dp), allocatable :: sd(:)
    real(dp), allocatable :: correlation(:, :)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: simulations(:, :)
    logical :: covariance_adjusted = .false.
    type(status_t) :: status
  end type multivariate_result

  type, public :: bond_result
    real(dp) :: price = 0.0_dp
    real(dp) :: duration = 0.0_dp
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: cash_flow(:)
    real(dp), allocatable :: discount_factor(:)
    real(dp), allocatable :: present_value(:)
    real(dp), allocatable :: duration_contribution(:)
    type(status_t) :: status
  end type bond_result

  type, public :: npv_result
    real(dp) :: value = 0.0_dp
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: cash_flow(:)
    real(dp), allocatable :: discount_factor(:)
    real(dp), allocatable :: present_value(:)
    type(status_t) :: status
  end type npv_result

  type, public :: irs_result
    real(dp) :: present_value = 0.0_dp
    real(dp) :: duration = 0.0_dp
    integer, allocatable :: payment_dates(:)
    real(dp), allocatable :: time(:)
    real(dp), allocatable :: fixed_leg(:)
    real(dp), allocatable :: floating_leg(:)
    real(dp), allocatable :: net(:)
    type(status_t) :: status
  end type irs_result

  type, public :: commodity_weight_result
    integer :: days_first = 0
    integer :: days_second = 0
    real(dp) :: first_weight = 0.0_dp
    type(status_t) :: status
  end type commodity_weight_result

  type, public :: frontier_result
    real(dp), allocatable :: weights(:, :)
    real(dp), allocatable :: expected_return(:)
    real(dp), allocatable :: risk(:)
    real(dp), allocatable :: sharpe(:)
    integer :: minimum_risk_index = 0
    integer :: maximum_sharpe_index = 0
    type(status_t) :: status
  end type frontier_result

  type, public :: lp_result
    real(dp) :: objective = 0.0_dp
    real(dp), allocatable :: solution(:)
    integer :: iterations = 0
    type(status_t) :: status
  end type lp_result

  type, public :: refinery_result
    real(dp) :: profit = 0.0_dp
    real(dp), allocatable :: slate(:)
    real(dp), allocatable :: margin(:)
    type(status_t) :: status
  end type refinery_result

  type, public :: beta_result
    real(dp), allocatable :: all(:)
    real(dp), allocatable :: bull(:)
    real(dp), allocatable :: bear(:)
    type(status_t) :: status
  end type beta_result

  type, public :: trade_stats_result
    real(dp) :: cumulative_return = 0.0_dp
    real(dp) :: annualized_return = 0.0_dp
    real(dp) :: annualized_sd = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: omega = 0.0_dp
    real(dp) :: fraction_winning = 0.0_dp
    real(dp) :: fraction_in_market = 0.0_dp
    integer :: maximum_drawdown_length = 0
    real(dp) :: maximum_drawdown = 0.0_dp
    type(status_t) :: status
  end type trade_stats_result

  type, public :: strategy_result
    real(dp), allocatable :: ret_close_close(:)
    real(dp), allocatable :: ret_open_close(:)
    real(dp), allocatable :: ret_close_open(:)
    real(dp), allocatable :: short_average(:)
    real(dp), allocatable :: long_average(:)
    integer, allocatable :: signal(:)
    integer, allocatable :: trade(:)
    integer, allocatable :: position(:)
    real(dp), allocatable :: strategy_return(:)
    real(dp), allocatable :: cumulative_equity(:)
    type(status_t) :: status
  end type strategy_result

end module rtl_types
