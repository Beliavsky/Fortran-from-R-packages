! SPDX-License-Identifier: LGPL-2.0-or-later
program test_invariants
   use ftrading
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   real(dp) :: x(80), high(80), low(80), close(80), volume(80), open(80), y(80), z(80), matrix(80, 3)
   real(dp) :: drawdown, qnan, inf, small(6), benchmark(5), scalar, huge_series(4)
   real(dp), allocatable :: rolling(:), expected(:), signals(:), durations(:)
   integer, allocatable :: from(:), to(:), indices(:)
   integer :: i, n, j
   procedure(rolling_operation), pointer :: callback
   qnan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   do i = 1, 80
      x(i) = real(i, dp) + sin(real(i, dp))
      close(i) = 50.0_dp + 0.2_dp * real(i, dp) + sin(real(i, dp))
      high(i) = close(i) + 2.0_dp
      low(i) = close(i) - 3.0_dp
      open(i) = close(i) - 0.3_dp
      volume(i) = 1000.0_dp + real(mod(i * 7, 31), dp)
   end do
   call max_draw_down([3.0_dp, 3.0_dp, 1.0_dp, 3.0_dp, 1.0_dp], drawdown, from, to)
   call check(drawdown == 2.0_dp, 'absolute drawdown')
   call check(all(from == [2, 4]) .and. all(to == [3, 5]), 'every tie with latest preceding peak')
   call max_draw_down([1.0_dp, 2.0_dp, 3.0_dp], drawdown, from, to)
   call check(drawdown == 0.0_dp .and. all(from == [1, 2, 3]) .and. all(to == from), 'zero drawdown ties')
   call max_draw_down([1.0_dp, qnan], drawdown, from, to)
   call check(ieee_is_nan(drawdown) .and. size(from) == 0, 'drawdown missing input rejected')
   call max_draw_down([real(dp) ::], drawdown, from, to)
   call check(ieee_is_nan(drawdown) .and. size(to) == 0, 'empty drawdown')
   benchmark = [0.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 5.0_dp]
   scalar = sharpe_ratio(benchmark, scale = 1.0_dp)
   call near([scalar], [1.25_dp / sqrt(8.75_dp / 3.0_dp)], 'Sharpe difference/sample sd')
   call near([sharpe_ratio(benchmark, rates = [0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp], scale = 1.0_dp)], &
      [1.0_dp / sqrt(8.75_dp / 3.0_dp)], 'time-varying risk-free returns')
   call near([sterling_ratio(benchmark)], [5.0_dp], 'Sterling absolute decline')
   call check(ieee_is_nan(sharpe_ratio([1.0_dp, 1.0_dp, 1.0_dp])), 'constant zero Sharpe')
   scalar = sharpe_ratio([1.0_dp, 2.0_dp, 3.0_dp])
   call check(.not. ieee_is_finite(scalar) .and. .not. ieee_is_nan(scalar), 'constant positive Sharpe infinity')
   call check(ieee_is_nan(sharpe_ratio([1.0_dp, 2.0_dp])), 'single return has undefined sample sd')
   call check(ieee_is_nan(sterling_ratio([1.0_dp])), 'singleton Sterling')
   call check(ieee_is_nan(sharpe_ratio(benchmark, rates = [0.0_dp])), 'rates length validation')
   do n = 1, 7
      rolling = roll_fun(x, n, 'mean')
      do i = 1, size(rolling)
         call near([rolling(i)], [sum(x(i:i + n - 1)) / real(n, dp)], 'independent rolling mean')
      end do
      rolling = roll_var(x, n)
      if (n == 1) then
         call check(all(ieee_is_nan(rolling)), 'singleton sample variance')
      else
         do i = 1, size(rolling)
            scalar = sum((x(i:i + n - 1) - sum(x(i:i + n - 1)) / real(n, dp)) ** 2) / real(n - 1, dp)
            call near([rolling(i)], [scalar], 'independent rolling sample variance')
         end do
      end if
   end do
   nullify(callback)
   rolling = roll_fun(x, 5, callback)
   call check(size(rolling) == 0, 'unassociated callback rejected')
   callback => span
   rolling = roll_fun(x, 5, callback)
   expected = roll_fun(x, 5, 'max') - roll_fun(x, 5, 'min')
   call near(rolling, expected, 'pure rolling callback')
   small = [1.0_dp, 2.0_dp, qnan, 4.0_dp, 5.0_dp, 6.0_dp]
   rolling = roll_fun(small, 2, 'mean', trim = .false.)
   call check(ieee_is_nan(rolling(1)) .and. all(ieee_is_nan(rolling(3:4))), 'window-local NA propagation')
   call near(rolling(5:), [4.5_dp, 5.5_dp], 'rolling recovers after interior NA')
   rolling = roll_fun(small, 2, 'mean', na_rm = .true.)
   call near(rolling, [1.5_dp, 3.0_dp, 4.5_dp, 5.5_dp], 'global NA removal before windows')
   rolling = roll_fun(x, 0, 'mean')
   call check(size(rolling) == 0, 'invalid zero window')
   rolling = roll_fun(x, 81, 'max')
   call check(size(rolling) == 0, 'window too long')
   rolling = roll_fun([1.0_dp, inf], 2, 'mean')
   call check(.not. ieee_is_finite(rolling(1)), 'infinite rolling mean')
   y = ema_ta(x, 0.25_dp, startup = 4)
   call near(y(:1), [sum(x(:4)) / 4.0_dp], 'EMA seed is first result')
   do i = 2, 80
      call near(y(i:i), [0.25_dp * x(i) + 0.75_dp * y(i - 1)], 'EMA recursion')
   end do
   call near(ema_ta(x, 3.0_dp), ema_ta(x, 0.5_dp), 'period-to-decay conversion')
   call near(ewma(x, 0.0_dp), spread(sum(x) / 80.0_dp, 1, 80), 'zero decay returns full-series mean')
   call check(all(ieee_is_nan(ema_ta(x, -1.0_dp))), 'invalid negative decay')
   call check(all(ieee_is_nan(ema_ta(x, 0.1_dp, startup = 81))), 'oversized seed')
   call check(all(ieee_is_nan(ema_ta([1.0_dp, 2.0_dp], 0.1_dp))), 'short default startup')
   huge_series = huge(1.0_dp) / 4.0_dp
   rolling = roll_fun(huge_series, 4, 'mean')
   call near(rolling / huge_series(1), [1.0_dp], 'large finite rolling mean')
   huge_series = huge(1.0_dp) * 0.75_dp
   call near(ema_ta(huge_series, 0.0_dp) / huge_series(1), spread(1.0_dp, 1, 4), 'overflow-safe constant mean')
   rolling = roll_var(huge_series, 4)
   call near(rolling, [0.0_dp], 'huge constant sample variance')
   rolling = roll_var([1.0e150_dp, 2.0e150_dp, 3.0e150_dp], 3)
   call near(rolling / 1.0e300_dp, [1.0_dp], 'large finite sample variance')
   y = mom_ta(close, 5)
   call near(y(:5), spread(0.0_dp, 1, 5), 'momentum zeros')
   call near(y(6:), close(6:) - close(:75), 'momentum differences')
   call near(roc_ta(close, 5) * close, y, 'ROC current price denominator')
   call near(wpr_ta(close, high, low, 5), fpk_ta(close, high, low, 5), 'upstream WPR equals FPK')
   y = fpk_ta(close, high, low, 5)
   do i = 1, 80
      j = max(1, i - 5)
      scalar = (close(i) - minval(low(j:i))) / (maxval(high(j:i)) - minval(low(j:i)))
      call near(y(i:i), [scalar], 'FPK partial startup and lag+1')
   end do
   call check(all(y >= 0.0_dp .and. y <= 1.0_dp), 'FPK range')
   y = rsi_ta(close, 5)
   call check(all(y >= 0.0_dp .and. y <= 1.0_dp), 'RSI fractional range')
   call near(rsi_ta(x * 0.0_dp + 1.0_dp, 5), spread(qnan, 1, 80), 'flat RSI undefined')
   y = cdo_ta(close, 3, 7, 4)
   call near(y, macd_ta(close, 3, 7) - cds_ta(close, lag3 = 4), 'CDO preserves signal defaults')
   matrix = bollinger_ta(close, 5)
   call near(matrix(:, 2), close, 'Bollinger middle is raw close')
   rolling = sma(close, 5)
   call near((matrix(5:, 1) + matrix(5:, 3)) / 2.0_dp, rolling, 'Bollinger bands centered on SMA')
   y = adi_ta(high, low, close, volume)
   call near(y(2:) - y(:79), (2.0_dp * close(2:) - high(2:) - low(2:)) / &
      (high(2:) - low(2:)) * volume(2:), 'accumulation increments')
   rolling = obv_ta([1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp], [5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp])
   call near(rolling, [0.0_dp, 0.0_dp, 7.0_dp, -1.0_dp, -1.0_dp], 'OBV exact ties and zero initial value')
   call near(nvi_ta(close, volume) + pvi_ta(close, volume), cumulative(roc_ta(close, 1) / 100.0_dp), &
      'volume-index partition on strictly changing volume')
   y = garmanklass_ta(open, high, low, close)
   z = garmanklass_ta(open * 2.0_dp, high * 2.0_dp, low * 2.0_dp, close * 2.0_dp)
   call check(all(z(2:) > y(2:)), 'raw-price Garman-Klass is scale-dependent upstream')
   call check(all(ieee_is_nan(fpk_ta(close, high(:5), low, 3))), 'mismatched indicator lengths')
   call check(all(ieee_is_nan(chaikinv_ta(high, low, 3, 100))), 'oversized volatility lag')
   call trade_signals([-1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], signals, indices)
   call near(signals, [1.0_dp, -1.0_dp], 'trade signals')
   call check(all(indices == [3, 5]), 'trade signal timestamps')
   durations = trade_lengths([1.0_dp, 4.0_dp, 10.0_dp])
   call near(durations, [3.0_dp, 6.0_dp], 'trade lengths')
   scalar = hit_rate([1.0_dp, -1.0_dp, 0.0_dp, qnan], [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp])
   call near([scalar], [0.833_dp], 'hit rate missing removal and three-digit rounding')
   rolling = mom_ta([real(dp) ::], 1)
   call check(size(rolling) == 0, 'empty momentum is safe')
   call trade_signals([real(dp) ::], signals, indices)
   call check(size(signals) == 0 .and. size(indices) == 0, 'empty trade signals')
   rolling = trade_lengths([real(dp) ::])
   call check(size(rolling) == 0, 'empty timestamp differences')
   call check(ieee_is_nan(hit_rate([real(dp) ::], [real(dp) ::])), 'empty hit rate')
   call near(fpk_ta([1.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], 1), &
      [qnan, qnan], 'zero high-low stochastic range')
   print *, 'fTrading independent invariants passed' 
contains
   subroutine check(ok, label)
      logical, intent(in) :: ok !! True when the independently stated invariant holds.
      character(len=*), intent(in) :: label !! Diagnostic identifying a failed invariant.
      if (.not. ok) then
         print *, 'FAILED: ', label
         error stop 1
      end if
   end subroutine
   subroutine near(actual, wanted, label)
      real(dp), intent(in) :: actual(:) !! Computed vector including possible missing observations.
      real(dp), intent(in) :: wanted(:) !! Independently derived expected vector with matching shape.
      character(len=*), intent(in) :: label !! Diagnostic identifying the comparison.
      integer :: k
      call check(size(actual) == size(wanted), label)
      do k = 1, size(actual)
         if (ieee_is_nan(wanted(k))) then
            call check(ieee_is_nan(actual(k)), label)
         else
            call check(ieee_is_finite(actual(k)), label)
            call check(abs(actual(k) - wanted(k)) <= 2.0e-11_dp * max(1.0_dp, abs(wanted(k))), label)
         end if
      end do
   end subroutine
   pure function span(window) result(value)
      real(dp), intent(in) :: window(:) !! Nonmissing sample window whose range is required.
      real(dp) :: value
      value = maxval(window) - minval(window)
   end function
   pure function cumulative(values) result(sums)
      real(dp), intent(in) :: values(:) !! Sequence to accumulate independently for a test oracle.
      real(dp) :: sums(size(values)), total
      integer :: k
      total = 0.0_dp
      do k = 1, size(values)
         total = total + values(k)
         sums(k) = total
      end do
   end function
end program test_invariants
