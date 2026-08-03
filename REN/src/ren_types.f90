! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_types
  use ren_kinds, only : dp
  implicit none
  private

  integer, parameter, public :: ren_success = 0
  integer, parameter, public :: ren_invalid_argument = 1
  integer, parameter, public :: ren_dimension_error = 2
  integer, parameter, public :: ren_numerical_error = 3
  integer, parameter, public :: ren_dependency_error = 4

  integer, parameter, public :: ren_method_count = 13
  character(len=24), parameter, public :: ren_method_names(ren_method_count) = [character(len=24) :: &
    'MV', 'JM', 'TRP_min', 'APP-Ridge', 'Fan.et.al_JM', 'FZY', 'EW', 'LW', &
    'TRP_clu', 'TZ', 'SW', 'SW_clu', 'SW_shrink']

  type, public :: asset_group
    integer, allocatable :: index(:)
  end type asset_group

  type, public :: cluster_result
    type(asset_group), allocatable :: groups(:)
    real(dp), allocatable :: merge_correlation(:)
    integer :: status = ren_success
  end type cluster_result

  type, public :: prepared_data_type
    real(dp), allocatable :: x(:, :)
    integer, allocatable :: month(:)
    integer, allocatable :: count(:)
    integer, allocatable :: date(:)
    integer, allocatable :: retained_columns(:)
    integer :: status = ren_success
  end type prepared_data_type

  type, public :: analysis_options
    integer :: cluster_repetitions = 100
    integer :: stochastic_samples = 1000
    integer :: random_seed = 1729
    real(dp) :: variance_tolerance = 1.0e-14_dp
    logical :: legacy_weight_timing = .true.
  end type analysis_options

  type, public :: analysis_result
    integer :: status = ren_success
    character(len=24), allocatable :: method(:)
    integer, allocatable :: month_date(:)
    integer, allocatable :: return_date(:)
    real(dp), allocatable :: weights(:, :, :)
    real(dp), allocatable :: turnover(:, :)
    real(dp), allocatable :: gross_returns(:, :)
    real(dp), allocatable :: cumulative_return(:, :)
    real(dp), allocatable :: wealth_index(:, :)
    real(dp), allocatable :: cumulative_turnover(:, :)
    real(dp), allocatable :: turnover_mean(:)
    real(dp), allocatable :: sharpe_ratio(:)
    real(dp), allocatable :: volatility(:)
    real(dp), allocatable :: max_drawdown(:)
    real(dp), allocatable :: vw_weights(:, :)
    real(dp), allocatable :: vw_turnover(:)
    real(dp), allocatable :: vw_gross_returns(:)
    real(dp), allocatable :: vw_cumulative_return(:)
    real(dp), allocatable :: vw_wealth_index(:)
    real(dp) :: vw_to_mean = 0.0_dp
    real(dp) :: vw_sharpe_ratio = 0.0_dp
    real(dp) :: vw_volatility = 0.0_dp
    real(dp) :: vw_max_drawdown = 0.0_dp
  end type analysis_result

  public :: ren_status_message
contains
  pure function ren_status_message(status) result(message)
    integer, intent(in) :: status
    character(len=:), allocatable :: message
    select case (status)
    case (ren_success)
      message = 'success'
    case (ren_invalid_argument)
      message = 'invalid argument'
    case (ren_dimension_error)
      message = 'dimension error'
    case (ren_numerical_error)
      message = 'numerical error'
    case (ren_dependency_error)
      message = 'dependency error'
    case default
      message = 'unknown status'
    end select
  end function ren_status_message
end module ren_types
