module pmwr_backtest
   use pmwr_kinds, only : dp
   use pmwr_types, only : backtest_result
   implicit none
   private
   public :: run_backtest

contains

   subroutine run_backtest(close_prices, targets, result, open_prices, signal_mask, rebalance_mask, &
                           convert_weights, trade_at_open, lag, burn_in, fraction, initial_position, &
                           initial_cash, transaction_cost, cashflow, tolerance)
      real(dp), intent(in) :: close_prices(:,:), targets(:,:)
      type(backtest_result), intent(out) :: result
      real(dp), intent(in), optional :: open_prices(:,:)
      logical, intent(in), optional :: signal_mask(:)
      logical, intent(in), optional :: rebalance_mask(:,:)
      logical, intent(in), optional :: convert_weights, trade_at_open
      integer, intent(in), optional :: lag, burn_in
      real(dp), intent(in), optional :: fraction
      real(dp), intent(in), optional :: initial_position(:), initial_cash
      real(dp), intent(in), optional :: transaction_cost(:), cashflow(:)
      real(dp), intent(in), optional :: tolerance

      real(dp), allocatable :: exec_price(:,:), tc(:), cf(:), desired(:), delta(:), init_pos(:)
      logical, allocatable :: sig(:), reb(:,:)
      logical :: weights, at_open
      real(dp) :: frac, cash0, cost, tol
      integer :: n, m, t, j, source, l, b, ntr, k

      n = size(close_prices, 1); m = size(close_prices, 2)
      if (any(shape(targets) /= shape(close_prices))) error stop "run_backtest: target shape mismatch"
      if (n < 1) error stop "run_backtest: no prices"
      allocate(exec_price(n, m)); exec_price = close_prices
      at_open = .true.; if (present(trade_at_open)) at_open = trade_at_open
      if (at_open .and. present(open_prices)) then
         if (any(shape(open_prices) /= shape(close_prices))) error stop "run_backtest: open shape mismatch"
         exec_price = open_prices
      end if
      allocate(sig(n)); sig = .true.
      if (present(signal_mask)) then
         if (size(signal_mask) /= n) error stop "run_backtest: signal mask size"
         sig = signal_mask
      end if
      allocate(reb(n, m)); reb = .true.
      if (present(rebalance_mask)) then
         if (any(shape(rebalance_mask) /= shape(close_prices))) error stop "run_backtest: rebalance shape mismatch"
         reb = rebalance_mask
      end if
      allocate(tc(m)); tc = 0.0_dp
      if (present(transaction_cost)) then
         if (size(transaction_cost) == 1) then
            tc = transaction_cost(1)
         else if (size(transaction_cost) == m) then
            tc = transaction_cost
         else
            error stop "run_backtest: transaction cost size"
         end if
      end if
      allocate(cf(n)); cf = 0.0_dp
      if (present(cashflow)) then
         if (size(cashflow) /= n) error stop "run_backtest: cashflow size"
         cf = cashflow
      end if
      allocate(init_pos(m)); init_pos = 0.0_dp
      if (present(initial_position)) then
         if (size(initial_position) /= m) error stop "run_backtest: initial position size"
         init_pos = initial_position
      end if
      cash0 = 0.0_dp; if (present(initial_cash)) cash0 = initial_cash
      weights = .false.; if (present(convert_weights)) weights = convert_weights
      frac = 1.0_dp; if (present(fraction)) frac = fraction
      l = 1; if (present(lag)) l = lag
      b = 0; if (present(burn_in)) b = burn_in
      tol = 1.0e-8_dp; if (present(tolerance)) tol = tolerance
      if (l < 0) error stop "run_backtest: negative lag"

      allocate(result%position(n, m), result%suggested_position(n, m), result%cash(n), result%wealth(n), &
               result%cumulative_cost(n))
      allocate(desired(m), delta(m))
      result%position = 0.0_dp; result%suggested_position = 0.0_dp
      result%cash = 0.0_dp; result%wealth = 0.0_dp; result%cumulative_cost = 0.0_dp
      result%position(1, :) = init_pos
      result%suggested_position(1, :) = init_pos
      result%cash(1) = cash0 + cf(1)
      result%wealth(1) = dot_product(init_pos, close_prices(1, :)) + result%cash(1)
      result%initial_wealth = result%wealth(1)

      do t = 2, n
         result%suggested_position(t, :) = result%suggested_position(t - 1, :)
         source = t - l
         if (source >= 1 .and. source <= n .and. sig(source)) then
            desired = targets(source, :)
            if (weights) then
               do j = 1, m
                  if (abs(close_prices(t - 1, j)) > epsilon(1.0_dp)) then
                     desired(j) = desired(j) * result%wealth(t - 1) / close_prices(t - 1, j)
                  else
                     desired(j) = 0.0_dp
                  end if
               end do
            end if
            result%suggested_position(t, :) = desired
         end if
         delta = result%suggested_position(t, :) - result%position(t - 1, :)
         if (t <= b) delta = 0.0_dp
         do j = 1, m
            if (.not. reb(t, j) .or. abs(delta(j)) < tol) delta(j) = 0.0_dp
         end do
         delta = frac * delta
         cost = sum(abs(delta) * tc * exec_price(t, :))
         result%position(t, :) = result%position(t - 1, :) + delta
         result%cash(t) = result%cash(t - 1) - dot_product(delta, exec_price(t, :)) - cost + cf(t)
         result%cumulative_cost(t) = result%cumulative_cost(t - 1) + cost
         result%wealth(t) = dot_product(result%position(t, :), close_prices(t, :)) + result%cash(t)
      end do

      ntr = 0
      do t = 1, n
         if (t == 1) then
            delta = result%position(1, :) - init_pos
         else
            delta = result%position(t, :) - result%position(t - 1, :)
         end if
         ntr = ntr + count(abs(delta) >= tol)
      end do
      result%journal%n = ntr
      allocate(result%journal%timestamp(ntr), result%journal%amount(ntr), result%journal%price(ntr), &
               result%journal%instrument(ntr), result%journal%account(ntr))
      k = 0
      do t = 2, n
         delta = result%position(t, :) - result%position(t - 1, :)
         do j = 1, m
            if (abs(delta(j)) >= tol) then
               k = k + 1
               result%journal%timestamp(k) = real(t, dp)
               result%journal%amount(k) = delta(j)
               result%journal%price(k) = exec_price(t, j)
               result%journal%instrument(k) = j
               result%journal%account(k) = 1
            end if
         end do
      end do
      result%journal%n = k
   end subroutine run_backtest

end module pmwr_backtest
