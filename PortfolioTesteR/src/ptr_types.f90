! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_types
  use ptr_kinds, only : dp
  implicit none
  private

  type, public :: backtest_result
    real(dp), allocatable :: portfolio_value(:)
    real(dp), allocatable :: returns(:)
    real(dp), allocatable :: cash(:)
    real(dp), allocatable :: shares(:,:)
    real(dp), allocatable :: executed_weights(:,:)
    real(dp), allocatable :: turnover(:)
    real(dp) :: initial_capital = 0.0_dp
    real(dp) :: total_return = 0.0_dp
    real(dp) :: annualized_return = 0.0_dp
    real(dp) :: annualized_volatility = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: max_drawdown = 0.0_dp
    integer :: n_transactions = 0
    integer :: first_trade = 0
    logical :: bankrupt = .false.
  end type backtest_result

  type, public :: performance_result
    real(dp) :: total_return = 0.0_dp
    real(dp) :: annualized_return = 0.0_dp
    real(dp) :: annualized_volatility = 0.0_dp
    real(dp) :: sharpe = 0.0_dp
    real(dp) :: sortino = 0.0_dp
    real(dp) :: max_drawdown = 0.0_dp
    real(dp) :: calmar = 0.0_dp
    real(dp) :: omega = 0.0_dp
    real(dp) :: var = 0.0_dp
    real(dp) :: cvar = 0.0_dp
    real(dp) :: win_rate = 0.0_dp
    real(dp) :: turnover = 0.0_dp
  end type performance_result

  type, public :: linear_model
    real(dp), allocatable :: coef(:)
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: scale(:)
    real(dp) :: lambda = 0.0_dp
    logical :: fitted = .false.
  end type linear_model

  type, public :: grid_result
    real(dp), allocatable :: params(:,:)
    real(dp), allocatable :: scores(:)
    integer :: best_index = 0
    real(dp) :: best_score = -huge(1.0_dp)
  end type grid_result

  type, public :: walk_forward_result
    integer, allocatable :: is_start(:), is_end(:), oos_start(:), oos_end(:)
    integer, allocatable :: best_index(:)
    real(dp), allocatable :: chosen_params(:,:)
    real(dp), allocatable :: is_score(:), oos_score(:), oos_return(:)
    real(dp), allocatable :: stitched_value(:)
    integer, allocatable :: stitched_index(:)
  end type walk_forward_result

end module ptr_types
