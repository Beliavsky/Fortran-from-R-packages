! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_optimization
  use ptr_kinds, only : dp
  use ptr_types, only : backtest_result, grid_result
  use ptr_backtest, only : run_backtest
  use ptr_utils, only : is_finite
  implicit none
  private
  public :: strategy_builder, metric_callback, run_param_grid, metric_sharpe

  abstract interface
    subroutine strategy_builder(prices, params, weights, status)
      import dp
      real(dp), intent(in) :: prices(:,:), params(:)
      real(dp), allocatable, intent(out) :: weights(:,:)
      integer, intent(out) :: status
    end subroutine strategy_builder

    real(dp) function metric_callback(backtest)
      import dp, backtest_result
      type(backtest_result), intent(in) :: backtest
    end function metric_callback
  end interface

contains

  real(dp) function metric_sharpe(backtest)
    type(backtest_result), intent(in) :: backtest
    metric_sharpe = backtest%sharpe
  end function metric_sharpe

  subroutine run_param_grid(prices, params, builder, result, metric, &
      initial_capital, cost_bps, frequency)
    real(dp), intent(in) :: prices(:,:), params(:,:)
    procedure(strategy_builder) :: builder
    type(grid_result), intent(out) :: result
    procedure(metric_callback), optional :: metric
    real(dp), intent(in), optional :: initial_capital, cost_bps, frequency
    real(dp), allocatable :: weights(:,:)
    type(backtest_result) :: bt
    real(dp) :: capital, cost, freq, score
    integer :: i, status
    capital = 100000.0_dp; if (present(initial_capital)) capital = initial_capital
    cost = 0.0_dp; if (present(cost_bps)) cost = cost_bps
    freq = 52.0_dp; if (present(frequency)) freq = frequency
    allocate(result%params(size(params,1),size(params,2)), result%scores(size(params,2)))
    result%params = params; result%scores = -huge(1.0_dp)
    result%best_index = 0; result%best_score = -huge(1.0_dp)
    do i = 1, size(params,2)
      call builder(prices, params(:,i), weights, status)
      if (status /= 0) cycle
      if (any(shape(weights) /= shape(prices))) cycle
      call run_backtest(prices, weights, capital, bt, cost_bps=cost, frequency=freq)
      if (present(metric)) then
        score = metric(bt)
      else
        score = metric_sharpe(bt)
      end if
      if (.not. is_finite(score)) cycle
      result%scores(i) = score
      if (result%best_index == 0 .or. score > result%best_score) then
        result%best_index = i; result%best_score = score
      end if
    end do
  end subroutine run_param_grid

end module ptr_optimization
