! SPDX-License-Identifier: GPL-2.0-or-later
program test_properties
   use ttr
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   real(dp) :: x(40), y(40), z(40), h(40, 3), a(40, 4), b(40, 3), v(40), missing(40)
   integer :: i, nmiss, begin
   logical :: valid
   x = [(real(i, dp), i = 1, 40)]
   v = 1.0_dp
   missing = x
   missing(:3) = nan_value()
   call na_check(missing, 5, nmiss, begin, valid)
   call require(valid .and. nmiss == 3 .and. begin == 8, 'Leading missing indices')
   missing(20) = nan_value()
   call na_check(missing, 5, nmiss, begin, valid)
   call require(.not. valid, 'Interior missing rejected by validation')
   call require(all(ieee_is_nan(sma(missing, 5))), 'Interior missing returns unavailable output')
   call require(all(ieee_is_nan(ema(x, 41))), 'Too-large EMA seed returns unavailable output')
   call require(all(ieee_is_nan(wma(x, 0))), 'Zero window returns unavailable output')
   call require(all(ieee_is_nan(run_cov(x, x(:39), 5))), 'Mismatched covariance shapes')
   call require(size(sma(x(:0), 3)) == 0, 'Empty series preserves empty shape')
   call close_vectors(run_sum(x, 5), 5.0_dp * run_mean(x, 5), 'Sum and mean identity')
   call close_vectors(run_cov(x, x, 5), run_var(x, 5), 'Covariance and variance identity')
   call close_vectors(run_sd(x, 5) ** 2, run_var(x, 5), 'Standard deviation identity')
   call close_vectors(sma(x, 5), sma(3.0_dp * x + 2.0_dp, 5) / 3.0_dp - 2.0_dp / 3.0_dp, &
      'SMA affine equivariance')
   call close_vectors(vwma(x, v, 5), sma(x, 5), 'Constant-volume VWMA equals SMA')
   call close_vectors(vwap(x, v, 5), vwma(x, v, 5), 'VWAP is the VWMA alias')
   call close_vectors(alma(v, 5), sma(v, 5), 'Normalized ALMA preserves constants')
   call close_vectors(wma(v, 5), sma(v, 5), 'Normalized WMA preserves constants')
   y = rsi(x, 5)
   call require(all(abs(y(6:) - 100.0_dp) < 1.0e-12_dp), 'Increasing RSI is 100')
   y = cmo(x, 5)
   call require(all(abs(y(6:) - 100.0_dp) < 1.0e-12_dp), 'Increasing CMO is 100')
   y = run_percent_rank(v, 5)
   call require(all(abs(y(5:) - 0.5_dp) < 1.0e-12_dp), 'All-tie rank is tie multiplier')
   y = run_percent_rank(x, 5)
   call require(all(abs(y(5:) - 0.9_dp) < 1.0e-12_dp), 'Increasing rank is upper half-step')
   y = ema(x, 5)
   call require(all(ieee_is_nan(y(:4))), 'EMA warm-up is unavailable')
   call require(abs(y(5) - 3.0_dp) < 1.0e-12_dp, 'EMA starts at sample mean')
   call require(abs(y(6) - 4.0_dp) < 1.0e-12_dp, 'EMA recurrence has 2/(n+1) decay')
   y = wilder_sum(x, 5)
   call require(abs(y(5) - 13.0_dp) < 1.0e-12_dp, 'Wilder sum uses upstream nontrivial seed')
   missing = x
   missing(:2) = nan_value()
   y = wilder_sum(missing, 1)
   call require(all(ieee_is_nan(y)), 'Wilder n=1 leading-NA upstream propagation quirk')
   y = run_percent_rank(missing, 1, cumulative = .true.)
   call require(abs(y(1) - 0.5_dp) < 1.0e-12_dp, 'Cumulative rank n=1 upstream seed quirk')
   call require(abs(y(2) - 0.25_dp) < 1.0e-12_dp, 'Cumulative rank counts missing prefix in denominator')
   h(:, 1) = x + 2.0_dp
   h(:, 2) = x - 1.0_dp
   h(:, 3) = x
   a = bbands(h, 5)
   call close_vectors(a(:, 2), (a(:, 1) + a(:, 3)) / 2.0_dp, 'Bollinger symmetry')
   b = donchian_channel(h(:, :2), 5)
   call close_vectors(b(:, 2), (b(:, 1) + b(:, 3)) / 2.0_dp, 'Donchian midpoint')
   b = tr(h)
   call require(all(ieee_is_nan(b(1, :))), 'True range needs previous close')
   call require(all(abs(b(2:, 1) - 3.0_dp) < 1.0e-12_dp), 'True range known non-gap result')
   a = atr(h, 5)
   call require(all(ieee_is_nan(a(:5, 2))), 'ATR has extra previous-close warm-up row')
   call require(all(abs(a(6:, 2) - 3.0_dp) < 1.0e-12_dp), 'Constant ATR invariant')
   h = 10.0_dp
   y = clv(h)
   call require(all(y == 0.0_dp), 'Flat-range CLV is zero')
   y = wpr(h, 5)
   call require(all(abs(y(5:) - 0.5_dp) < 1.0e-12_dp), 'Flat-range WPR is one half')
   y = mfi(h, v, 5)
   call require(all(abs(y(6:) - 50.0_dp) < 1.0e-12_dp), 'Unchanged-price MFI is fifty')
   y = obv(x, v)
   call close_vectors(y, x, 'Rising unit-volume OBV equals elapsed count')
   y = obv(41.0_dp - x, v)
   call close_vectors(y, 2.0_dp - x, 'Falling OBV includes first volume positively')
   y = momentum(x, 3)
   call require(all(abs(y(4:) - 3.0_dp) < 1.0e-12_dp), 'Momentum lag difference')
   z = roc(2.0_dp ** (x / 10.0_dp), 5)
   call require(all(abs(z(6:) - log(2.0_dp) / 2.0_dp) < 1.0e-12_dp), 'Log return on geometric prices')
   print *, 'TTR independent edge and algebraic property tests passed'
contains
   subroutine require(condition, label)
      logical, intent(in) :: condition !! Property that must hold for this deterministic test.
      character(*), intent(in) :: label !! Human-readable description printed on failure.
      if (condition) return
      print *, label
      error stop 'TTR property failed'
   end subroutine

   subroutine close_vectors(actual, expected, label)
      real(dp), intent(in) :: actual(:) !! Computed series to compare element by element.
      real(dp), intent(in) :: expected(:) !! Independent formula, including matching NaN padding.
      character(*), intent(in) :: label !! Name of the algebraic property being verified.
      integer :: j
      call require(size(actual) == size(expected), label)
      do j = 1, size(actual)
         if (ieee_is_nan(expected(j))) then
            call require(ieee_is_nan(actual(j)), label)
         else
            call require(.not. ieee_is_nan(actual(j)), label)
            call require(abs(actual(j) - expected(j)) <= 1.0e-10_dp * max(1.0_dp, abs(expected(j))), label)
         end if
      end do
   end subroutine
end program test_properties
