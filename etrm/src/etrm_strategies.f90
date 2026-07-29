! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm_strategies
   use etrm_kinds, only : dp
   use etrm_math, only : normal_cdf, signum, trade_round
   use etrm_status, only : etrm_ok, etrm_err_size, etrm_err_argument, set_status
   use etrm_types, only : strategy_result
   implicit none
   private

   public :: cppi, dppi, obpi, shpi, slpi

   interface cppi
      module procedure cppi_scalar
      module procedure cppi_vector
   end interface cppi

   interface dppi
      module procedure dppi_scalar
      module procedure dppi_vector
   end interface dppi

contains

   subroutine cppi_scalar(q, f, target_percent, risk_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent, risk_percent
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp) :: rp(1)

      rp(1) = risk_percent
      call cppi_core(q, f, target_percent, rp, result, status, message, transaction_cost, integer_trades)
   end subroutine cppi_scalar

   subroutine cppi_vector(q, f, target_percent, risk_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent, risk_percent(:)
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades

      call cppi_core(q, f, target_percent, risk_percent, result, status, message, &
         transaction_cost, integer_trades)
   end subroutine cppi_vector

   subroutine cppi_core(q, f, target_percent, risk_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent, risk_percent(:)
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp) :: tcost, target_price, gap, desired, cost
      logical :: ints
      integer :: i, n

      call strategy_options(transaction_cost, integer_trades, tcost, ints)
      if (.not. validate_common(q, f, target_percent, tcost, status, message)) return
      n = size(f)
      if (size(risk_percent) /= 1 .and. size(risk_percent) /= n) then
         call set_status(status, message, etrm_err_size, &
            "CPPI risk_percent must have length one or match the price vector")
         return
      end if
      if (any(risk_percent <= 0.0_dp)) then
         call set_status(status, message, etrm_err_argument, "CPPI risk percentages must be positive")
         return
      end if

      call allocate_strategy(result, "CPPI", q, f, tcost, ints, .true.)
      target_price = f(1) * (1.0_dp + target_percent)
      result%target = target_price
      if (size(risk_percent) == 1) then
         result%risk_factor = f(1) * risk_percent(1)
      else
         result%risk_factor = f(1) * risk_percent
      end if

      result%portfolio(1) = f(1)
      gap = directional_gap(q, target_price, result%portfolio(1))
      if (gap > result%risk_factor(1)) then
         desired = 0.0_dp
      else
         desired = (1.0_dp - gap / result%risk_factor(1)) * q
      end if
      result%position(1) = trade_round(desired, ints)
      result%trade(1) = result%position(1)
      result%exposed(1) = q - result%position(1)
      result%hedge(1) = abs(result%position(1) / q)
      cost = (f(1) + tcost * signum(result%trade(1))) * result%trade(1)
      result%portfolio(1) = (f(1) * result%exposed(1) + cost) / q

      do i = 2, n
         gap = directional_gap(q, target_price, result%portfolio(i - 1))
         if (gap > result%risk_factor(i - 1)) then
            desired = 0.0_dp
         else if (gap < 0.0_dp) then
            desired = q
         else
            desired = (1.0_dp - gap / result%risk_factor(i - 1)) * q
         end if
         result%position(i) = trade_round(desired, ints)
         result%trade(i) = result%position(i) - result%position(i - 1)
         result%exposed(i) = q - result%position(i)
         result%hedge(i) = abs(result%position(i) / q)
         cost = cost + (f(i) + tcost * signum(result%trade(i))) * result%trade(i)
         result%portfolio(i) = (f(i) * result%exposed(i) + cost) / q
      end do

      call set_status(status, message, etrm_ok, "")
   end subroutine cppi_core

   subroutine dppi_scalar(q, f, target_percent, risk_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent, risk_percent
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp) :: rp(1)

      rp(1) = risk_percent
      call dppi_core(q, f, target_percent, rp, result, status, message, transaction_cost, integer_trades)
   end subroutine dppi_scalar

   subroutine dppi_vector(q, f, target_percent, risk_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent, risk_percent(:)
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades

      call dppi_core(q, f, target_percent, risk_percent, result, status, message, &
         transaction_cost, integer_trades)
   end subroutine dppi_vector

   subroutine dppi_core(q, f, target_percent, risk_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent, risk_percent(:)
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp) :: tcost, gap, desired, cost
      logical :: ints
      integer :: i, n

      call strategy_options(transaction_cost, integer_trades, tcost, ints)
      if (.not. validate_common(q, f, target_percent, tcost, status, message)) return
      n = size(f)
      if (size(risk_percent) /= 1 .and. size(risk_percent) /= n) then
         call set_status(status, message, etrm_err_size, &
            "DPPI risk_percent must have length one or match the price vector")
         return
      end if
      if (any(risk_percent <= 0.0_dp)) then
         call set_status(status, message, etrm_err_argument, "DPPI risk percentages must be positive")
         return
      end if

      call allocate_strategy(result, "DPPI", q, f, tcost, ints, .true.)
      if (size(risk_percent) == 1) then
         result%risk_factor = f(1) * risk_percent(1)
      else
         result%risk_factor = f * risk_percent
      end if

      result%portfolio(1) = f(1)
      result%target(1) = f(1) * (1.0_dp + target_percent)
      gap = directional_gap(q, result%target(1), result%portfolio(1))
      if (gap > result%risk_factor(1)) then
         desired = 0.0_dp
      else
         desired = (1.0_dp - gap / result%risk_factor(1)) * q
      end if
      result%position(1) = trade_round(desired, ints)
      result%trade(1) = result%position(1)
      result%exposed(1) = q - result%position(1)
      result%hedge(1) = abs(result%position(1) / q)
      cost = (f(1) + tcost * signum(result%trade(1))) * result%trade(1)
      result%portfolio(1) = (f(1) * result%exposed(1) + cost) / q

      do i = 2, n
         gap = directional_gap(q, result%target(i - 1), result%portfolio(i - 1))
         if (gap > result%risk_factor(i - 1)) then
            desired = 0.0_dp
         else if (gap < 0.0_dp) then
            desired = q
         else
            desired = (1.0_dp - gap / result%risk_factor(i - 1)) * q
         end if
         result%position(i) = trade_round(desired, ints)
         result%trade(i) = result%position(i) - result%position(i - 1)
         result%exposed(i) = q - result%position(i)
         result%hedge(i) = abs(result%position(i) / q)
         cost = cost + (f(i) + tcost * signum(result%trade(i))) * result%trade(i)
         result%portfolio(i) = (f(i) * result%exposed(i) + cost) / q
         if (target_percent < 0.0_dp) then
            result%target(i) = max(result%portfolio(i) * (1.0_dp + target_percent), &
               result%target(i - 1))
         else
            result%target(i) = min(result%portfolio(i) * (1.0_dp + target_percent), &
               result%target(i - 1))
         end if
      end do

      call set_status(status, message, etrm_ok, "")
   end subroutine dppi_core

   subroutine obpi(q, f, volatility, days_left, result, status, message, strike_price, &
      interest_rate, trading_days, transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), volatility
      integer, intent(in) :: days_left
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: strike_price, interest_rate, transaction_cost
      integer, intent(in), optional :: trading_days
      logical, intent(in), optional :: integer_trades
      real(dp) :: k, r, tcost, t, d1, d2, premium, cost
      logical :: ints
      integer :: i, n, tdays

      call strategy_options(transaction_cost, integer_trades, tcost, ints)
      n = size(f)
      if (n < 1) then
         call set_status(status, message, etrm_err_size, "Price vector must not be empty")
         return
      end if
      if (abs(q) <= tiny(1.0_dp)) then
         call set_status(status, message, etrm_err_argument, "Volume must be nonzero")
         return
      end if
      if (tcost < 0.0_dp) then
         call set_status(status, message, etrm_err_argument, "Transaction cost cannot be negative")
         return
      end if
      if (volatility <= 0.0_dp) then
         call set_status(status, message, etrm_err_argument, "Volatility must be positive")
         return
      end if
      if (any(f <= 0.0_dp)) then
         call set_status(status, message, etrm_err_argument, "Futures prices must be positive")
         return
      end if
      k = f(1)
      if (present(strike_price)) k = strike_price
      if (k <= 0.0_dp) then
         call set_status(status, message, etrm_err_argument, "Strike price must be positive")
         return
      end if
      r = 0.0_dp
      if (present(interest_rate)) r = interest_rate
      tdays = 250
      if (present(trading_days)) tdays = trading_days
      if (tdays <= 0) then
         call set_status(status, message, etrm_err_argument, "Trading days per year must be positive")
         return
      end if
      if (days_left < n) then
         call set_status(status, message, etrm_err_argument, &
            "days_left must be at least the number of price observations")
         return
      end if

      call allocate_strategy(result, "OBPI", q, f, tcost, ints, .false.)
      result%strike_price = k
      result%annual_volatility = volatility
      result%interest_rate = r
      result%trading_days = tdays

      t = real(days_left, dp) / real(tdays, dp)
      d1 = (log(f(1) / k) + 0.5_dp * volatility**2 * t) / (volatility * sqrt(t))
      d2 = d1 - volatility * sqrt(t)
      if (q > 0.0_dp) then
         premium = exp(-r * t) * (f(1) * normal_cdf(d1) - k * normal_cdf(d2))
         result%target = k + premium
      else
         premium = exp(-r * t) * (k * normal_cdf(-d2) - f(1) * normal_cdf(-d1))
         result%target = k - premium
      end if

      do i = 1, n
         t = real(days_left - i + 1, dp) / real(tdays, dp)
         d1 = (log(f(i) / k) + 0.5_dp * volatility**2 * t) / (volatility * sqrt(t))
         if (q > 0.0_dp) then
            result%position(i) = trade_round(q * normal_cdf(d1) * exp(-r * t), ints)
         else
            result%position(i) = trade_round(q * normal_cdf(-d1) * exp(-r * t), ints)
         end if
      end do

      cost = 0.0_dp
      do i = 1, n
         if (i == 1) then
            result%trade(i) = result%position(i)
         else
            result%trade(i) = result%position(i) - result%position(i - 1)
         end if
         result%exposed(i) = q - result%position(i)
         result%hedge(i) = result%position(i) / q
         cost = cost + (f(i) + signum(result%trade(i)) * tcost) * result%trade(i)
         result%portfolio(i) = (cost + result%exposed(i) * f(i)) / q
      end do

      call set_status(status, message, etrm_ok, "")
   end subroutine obpi

   subroutine slpi(q, f, target_percent, result, status, message, transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp) :: tcost, target_price
      logical :: ints
      integer :: i, crossing, n

      call strategy_options(transaction_cost, integer_trades, tcost, ints)
      if (.not. validate_common(q, f, target_percent, tcost, status, message)) return
      n = size(f)
      call allocate_strategy(result, "SLPI", q, f, tcost, ints, .false.)
      target_price = f(1) * (1.0_dp + target_percent)
      result%target = target_price
      result%position = 0.0_dp
      call recompute_path(q, f, tcost, result)

      crossing = 0
      if (q > 0.0_dp) then
         do i = 1, n
            if (result%portfolio(i) >= target_price) then
               crossing = i
               exit
            end if
         end do
      else
         do i = 1, n
            if (result%portfolio(i) <= target_price) then
               crossing = i
               exit
            end if
         end do
      end if
      if (crossing > 0) then
         result%position(crossing:n) = q
         call recompute_path(q, f, tcost, result)
      end if

      call set_status(status, message, etrm_ok, "")
   end subroutine slpi

   subroutine shpi(q, f, days_left, target_percent, result, status, message, &
      transaction_cost, integer_trades)
      real(dp), intent(in) :: q, f(:), target_percent
      integer, intent(in) :: days_left
      type(strategy_result), intent(out) :: result
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp) :: tcost, target_price
      logical :: ints
      integer :: i, crossing, n

      call strategy_options(transaction_cost, integer_trades, tcost, ints)
      if (.not. validate_common(q, f, target_percent, tcost, status, message)) return
      n = size(f)
      if (days_left < n) then
         call set_status(status, message, etrm_err_argument, &
            "days_left must be at least the number of price observations")
         return
      end if

      call allocate_strategy(result, "SHPI", q, f, tcost, ints, .false.)
      target_price = f(1) * (1.0_dp + target_percent)
      result%target = target_price
      do i = 1, n
         result%position(i) = trade_round(real(i, dp) / real(days_left, dp) * q, ints)
      end do
      call recompute_path(q, f, tcost, result)

      crossing = 0
      if (q > 0.0_dp) then
         do i = 1, n
            if (result%portfolio(i) >= target_price) then
               crossing = i
               exit
            end if
         end do
      else
         do i = 1, n
            if (result%portfolio(i) <= target_price) then
               crossing = i
               exit
            end if
         end do
      end if
      if (crossing > 0) then
         result%position(crossing:n) = q
         call recompute_path(q, f, tcost, result)
      end if

      call set_status(status, message, etrm_ok, "")
   end subroutine shpi

   subroutine allocate_strategy(result, name, q, f, transaction_cost, integer_trades, with_risk)
      type(strategy_result), intent(out) :: result
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: q, f(:), transaction_cost
      logical, intent(in) :: integer_trades, with_risk
      integer :: n

      n = size(f)
      result%name = name
      result%volume = q
      result%transaction_cost = transaction_cost
      result%integer_trades = integer_trades
      allocate(result%market(n), result%trade(n), result%exposed(n), result%position(n), &
         result%hedge(n), result%target(n), result%portfolio(n))
      if (with_risk) allocate(result%risk_factor(n))
      result%market = f
      result%trade = 0.0_dp
      result%exposed = 0.0_dp
      result%position = 0.0_dp
      result%hedge = 0.0_dp
      result%target = 0.0_dp
      result%portfolio = 0.0_dp
      if (with_risk) result%risk_factor = 0.0_dp
   end subroutine allocate_strategy

   subroutine recompute_path(q, f, transaction_cost, result)
      real(dp), intent(in) :: q, f(:), transaction_cost
      type(strategy_result), intent(inout) :: result
      real(dp) :: cost
      integer :: i

      cost = 0.0_dp
      do i = 1, size(f)
         if (i == 1) then
            result%trade(i) = result%position(i)
         else
            result%trade(i) = result%position(i) - result%position(i - 1)
         end if
         result%exposed(i) = q - result%position(i)
         result%hedge(i) = result%position(i) / q
         cost = cost + (f(i) + transaction_cost * signum(result%trade(i))) * result%trade(i)
         result%portfolio(i) = (cost + result%exposed(i) * f(i)) / q
      end do
   end subroutine recompute_path

   subroutine strategy_options(transaction_cost, integer_trades, tcost, ints)
      real(dp), intent(in), optional :: transaction_cost
      logical, intent(in), optional :: integer_trades
      real(dp), intent(out) :: tcost
      logical, intent(out) :: ints

      tcost = 0.0_dp
      if (present(transaction_cost)) tcost = transaction_cost
      ints = .true.
      if (present(integer_trades)) ints = integer_trades
   end subroutine strategy_options

   logical function validate_common(q, f, target_percent, transaction_cost, status, message) result(valid)
      real(dp), intent(in) :: q, f(:), target_percent, transaction_cost
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message

      valid = .false.
      if (size(f) < 1) then
         call set_status(status, message, etrm_err_size, "Price vector must not be empty")
      else if (abs(q) <= tiny(1.0_dp)) then
         call set_status(status, message, etrm_err_argument, "Volume must be nonzero")
      else if (abs(target_percent) <= tiny(1.0_dp)) then
         call set_status(status, message, etrm_err_argument, "Target percentage cannot be zero")
      else if (transaction_cost < 0.0_dp) then
         call set_status(status, message, etrm_err_argument, "Transaction cost cannot be negative")
      else if (q < 0.0_dp .and. target_percent > 0.0_dp) then
         call set_status(status, message, etrm_err_argument, &
            "A seller cannot set a target above the initial market price")
      else if (q > 0.0_dp .and. target_percent < 0.0_dp) then
         call set_status(status, message, etrm_err_argument, &
            "A buyer cannot set a target below the initial market price")
      else
         valid = .true.
      end if
   end function validate_common

   pure real(dp) function directional_gap(q, target, portfolio) result(gap)
      real(dp), intent(in) :: q, target, portfolio
      if (q > 0.0_dp) then
         gap = target - portfolio
      else
         gap = portfolio - target
      end if
   end function directional_gap

end module etrm_strategies
