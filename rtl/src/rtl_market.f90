! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_market
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use rtl_kinds, only: dp
  use rtl_types, only: beta_result, trade_stats_result, strategy_result
  use rtl_stats, only: beta_value, mean_value, sample_sd
  implicit none
  private

  public :: compute_returns, roll_adjust_mask, prompt_beta, trade_stats
  public :: moving_average_strategy, trade_strategy_sma, trade_strategy_dy
  public :: returns, rolladjust, promptBeta, tradeStats, tradeStrategySMA, tradeStrategyDY

  interface returns
    module procedure compute_returns
  end interface returns

  interface rolladjust
    module procedure roll_adjust_mask
  end interface rolladjust

  interface promptBeta
    module procedure prompt_beta
  end interface promptBeta

  interface tradeStats
    module procedure trade_stats
  end interface tradeStats

  interface tradeStrategySMA
    module procedure trade_strategy_sma
  end interface tradeStrategySMA

  interface tradeStrategyDY
    module procedure trade_strategy_dy
  end interface tradeStrategyDY

contains

  function compute_returns(values, period_return, return_type) result(output)
    real(dp), intent(in) :: values(:, :)
    integer, intent(in) :: period_return
    character(len=*), intent(in) :: return_type
    real(dp), allocatable :: output(:, :)
    integer :: n, p, i
    character(len=:), allocatable :: kind_name

    n = size(values, 1)
    p = size(values, 2)
    if (period_return < 1 .or. period_return >= n) then
      allocate(output(0, p))
      return
    end if
    allocate(output(n - period_return, p))
    kind_name = lowercase(trim(return_type))
    select case (kind_name)
    case ("abs")
      output = values(period_return + 1:n, :) - values(1:n - period_return, :)
    case ("rel")
      output = values(period_return + 1:n, :) / values(1:n - period_return, :) - 1.0_dp
    case ("log")
      output = log(values(period_return + 1:n, :) / values(1:n - period_return, :))
    case default
      output = 0.0_dp
    end select
    do i = 1, size(output, 1)
      where (.not. ieee_is_finite(output(i, :))) output(i, :) = 0.0_dp
    end do
  end function compute_returns

  pure function roll_adjust_mask(dates, expiry_dates) result(keep)
    integer, intent(in) :: dates(:), expiry_dates(:)
    logical :: keep(size(dates))
    integer :: i
    keep = .true.
    do i = 2, size(dates)
      if (any(dates(i - 1) == expiry_dates)) keep(i) = .false.
    end do
  end function roll_adjust_mask

  function prompt_beta(return_matrix) result(output)
    real(dp), intent(in) :: return_matrix(:, :)
    type(beta_result) :: output
    logical, allocatable :: bull_mask(:), bear_mask(:)
    integer :: p, i

    p = size(return_matrix, 2)
    if (size(return_matrix, 1) < 2 .or. p < 1) then
      output%status%ok = .false.
      output%status%message = "prompt_beta requires at least two rows and one column"
      return
    end if
    allocate(output%all(p), output%bull(p), output%bear(p))
    allocate(bull_mask(size(return_matrix, 1)), bear_mask(size(return_matrix, 1)))
    bull_mask = return_matrix(:, 1) > 0.0_dp
    bear_mask = return_matrix(:, 1) < 0.0_dp
    do i = 1, p
      output%all(i) = beta_value(return_matrix(:, i), return_matrix(:, 1))
      output%bull(i) = beta_value(return_matrix(:, i), return_matrix(:, 1), bull_mask)
      output%bear(i) = beta_value(return_matrix(:, i), return_matrix(:, 1), bear_mask)
    end do
  end function prompt_beta

  function trade_stats(return_series, risk_free_rate, periods_per_year) result(output)
    real(dp), intent(in) :: return_series(:)
    real(dp), intent(in), optional :: risk_free_rate
    integer, intent(in), optional :: periods_per_year
    type(trade_stats_result) :: output
    real(dp) :: wealth, annual_scale, rf, peak, drawdown, gains, losses, threshold
    integer :: n, scale, i, current_length, nonzero_count, win_count

    n = size(return_series)
    if (n < 1 .or. any(return_series <= -1.0_dp)) then
      output%status%ok = .false.
      output%status%message = "returns must be nonempty and greater than -1"
      return
    end if
    scale = 252
    if (present(periods_per_year)) scale = periods_per_year
    rf = 0.0_dp
    if (present(risk_free_rate)) rf = risk_free_rate
    wealth = product(1.0_dp + return_series)
    output%cumulative_return = wealth - 1.0_dp
    annual_scale = real(scale, dp) / real(n, dp)
    output%annualized_return = wealth**annual_scale - 1.0_dp
    output%annualized_sd = sample_sd(return_series) * sqrt(real(scale, dp))
    if (output%annualized_sd > 0.0_dp) then
      output%sharpe = (output%annualized_return - rf) / output%annualized_sd
    end if
    threshold = rf / real(scale, dp)
    gains = sum(max(return_series - threshold, 0.0_dp))
    losses = sum(max(threshold - return_series, 0.0_dp))
    if (losses > 0.0_dp) output%omega = gains / losses
    nonzero_count = count(abs(return_series) > epsilon(1.0_dp))
    win_count = count(return_series > 0.0_dp)
    if (nonzero_count > 0) output%fraction_winning = real(win_count, dp) / real(nonzero_count, dp)
    output%fraction_in_market = real(nonzero_count, dp) / real(n, dp)

    wealth = 1.0_dp
    peak = 1.0_dp
    current_length = 0
    output%maximum_drawdown = 0.0_dp
    do i = 1, n
      wealth = wealth * (1.0_dp + return_series(i))
      if (wealth >= peak) then
        peak = wealth
        current_length = 0
      else
        current_length = current_length + 1
        output%maximum_drawdown_length = max(output%maximum_drawdown_length, current_length)
      end if
      drawdown = wealth / peak - 1.0_dp
      output%maximum_drawdown = min(output%maximum_drawdown, drawdown)
    end do
  end function trade_stats

  function moving_average_strategy(open_price, close_price, short_window, long_window) result(output)
    real(dp), intent(in) :: open_price(:), close_price(:)
    integer, intent(in) :: short_window, long_window
    type(strategy_result) :: output
    integer :: n, i
    real(dp) :: ret_new, ret_existing, ret_other

    n = min(size(open_price), size(close_price))
    if (n < 2 .or. short_window < 1 .or. long_window < 1) then
      output%status%ok = .false.
      output%status%message = "invalid moving-average strategy inputs"
      return
    end if
    allocate(output%ret_close_close(n), output%ret_open_close(n), output%ret_close_open(n))
    allocate(output%short_average(n), output%long_average(n), output%signal(n))
    allocate(output%trade(n), output%position(n), output%strategy_return(n))
    allocate(output%cumulative_equity(n))
    output%ret_close_close = 0.0_dp
    output%ret_open_close = 0.0_dp
    output%ret_close_open = 0.0_dp
    output%short_average = 0.0_dp
    output%long_average = 0.0_dp
    output%signal = 0
    output%trade = 0
    output%position = 0
    output%strategy_return = 0.0_dp
    output%cumulative_equity = 1.0_dp

    do i = 2, n
      output%ret_close_close(i) = close_price(i) / close_price(i - 1) - 1.0_dp
      output%ret_open_close(i) = (close_price(i) - open_price(i)) / close_price(i)
      output%ret_close_open(i) = open_price(i) / close_price(i - 1) - 1.0_dp
    end do
    call simple_moving_average(close_price(1:n), short_window, output%short_average)
    call simple_moving_average(close_price(1:n), long_window, output%long_average)
    do i = 1, n
      if (i >= max(short_window, long_window)) then
        if (output%short_average(i) > output%long_average(i)) output%signal(i) = 1
        if (output%short_average(i) < output%long_average(i)) output%signal(i) = -1
      end if
      if (i >= 3) output%trade(i) = output%signal(i - 1) - output%signal(i - 2)
      if (i == 1) then
        output%position(i) = output%trade(i)
      else
        output%position(i) = output%position(i - 1) + output%trade(i)
      end if
      ret_new = 0.0_dp
      ret_existing = 0.0_dp
      ret_other = 0.0_dp
      if (output%position(i) == output%trade(i)) then
        ret_new = real(output%position(i), dp) * output%ret_open_close(i)
      end if
      if (output%position(i) /= 0 .and. output%trade(i) == 0) then
        ret_existing = real(output%position(i), dp) * output%ret_close_close(i)
      end if
      if ((output%position(i) - output%trade(i)) /= 0 .and. output%trade(i) /= 0) then
        ret_other = (1.0_dp + output%ret_close_open(i) * &
          real(output%position(i) - output%trade(i), dp)) * &
          (1.0_dp + output%ret_open_close(i) * real(output%position(i), dp)) - 1.0_dp
      end if
      output%strategy_return(i) = ret_new + ret_existing + ret_other
      if (i == 1) then
        output%cumulative_equity(i) = 1.0_dp + output%strategy_return(i)
      else
        output%cumulative_equity(i) = output%cumulative_equity(i - 1) * &
          (1.0_dp + output%strategy_return(i))
      end if
    end do
  end function moving_average_strategy

  function trade_strategy_sma(open_price, close_price, short_window, long_window) result(output)
    real(dp), intent(in) :: open_price(:), close_price(:)
    integer, intent(in) :: short_window, long_window
    type(strategy_result) :: output
    output = moving_average_strategy(open_price, close_price, short_window, long_window)
  end function trade_strategy_sma

  function trade_strategy_dy(open_price, close_price, short_window, long_window) result(output)
    real(dp), intent(in) :: open_price(:), close_price(:)
    integer, intent(in) :: short_window, long_window
    type(strategy_result) :: output
    ! Upstream tradeStrategyDY is identical to tradeStrategySMA and does not use dividend yield.
    output = moving_average_strategy(open_price, close_price, short_window, long_window)
  end function trade_strategy_dy

  pure subroutine simple_moving_average(x, window, average)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: window
    real(dp), intent(out) :: average(:)
    integer :: i
    average = 0.0_dp
    do i = window, size(x)
      average(i) = mean_value(x(i - window + 1:i))
    end do
  end subroutine simple_moving_average

  pure function lowercase(text) result(output)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: output
    integer :: i, code
    output = text
    do i = 1, len(text)
      code = iachar(output(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code + 32)
    end do
  end function lowercase

end module rtl_market
