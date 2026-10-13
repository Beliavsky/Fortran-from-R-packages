! SPDX-License-Identifier: LGPL-2.0-or-later
! Derived from fTrading 3042.79 by Diethelm Wuertz and contributors.
module ftrading_indicators
   use r_kinds, only: dp
   use ftrading_base, only: ema_ta, sma, roll_fun, roll_var
   use ttr_base, only: nan_value, divide, ttr_momentum => momentum
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: bias_ta, medprice_ta, typicalprice_ta, wclose_ta, roc_ta, osc_ta, mom_ta
   public :: macd_ta, cds_ta, cdo_ta, vohl_ta, vor_ta, stochastic_ta, fpk_ta, fpd_ta
   public :: spd_ta, apd_ta, wpr_ta, rsi_ta, accel_ta, adi_ta, adoscillator_ta
   public :: bollinger_ta, chaikino_ta, chaikinv_ta, garmanklass_ta, nvi_ta, obv_ta
   public :: pvi_ta, pvtrend_ta, williamsad_ta, williamsr_ta, daily_ta
contains
   pure elemental function medprice_ta(high, low) result(y)
      real(dp), intent(in) :: high !! High price for the same observation; missing values propagate.
      real(dp), intent(in) :: low !! Low price for the same observation; missing values propagate.
      real(dp) :: y
      y = (high + low) / 2.0_dp
   end function

   pure elemental function typicalprice_ta(high, low, close) result(y)
      real(dp), intent(in) :: high !! High price for the same observation; missing values propagate.
      real(dp), intent(in) :: low !! Low price for the same observation; missing values propagate.
      real(dp), intent(in) :: close !! Close price for the same observation; missing values propagate.
      real(dp) :: y
      y = (high + low + close) / 3.0_dp
   end function

   pure elemental function wclose_ta(high, low, close) result(y)
      real(dp), intent(in) :: high !! High price for the same observation; missing values propagate.
      real(dp), intent(in) :: low !! Low price for the same observation; missing values propagate.
      real(dp), intent(in) :: close !! Close price for the same observation; missing values propagate.
      real(dp) :: y
      y = (high + low + 2.0_dp * close) / 4.0_dp
   end function

   pure elemental function vohl_ta(high, low) result(y)
      real(dp), intent(in) :: high !! High price for the same observation; missing values propagate.
      real(dp), intent(in) :: low !! Low price for the same observation; missing values propagate.
      real(dp) :: y
      y = high - low
   end function

   pure elemental function vor_ta(high, low) result(y)
      real(dp), intent(in) :: high !! High price for the same observation; missing values propagate.
      real(dp), intent(in) :: low !! Low price for the same observation; missing values propagate.
      real(dp) :: y
      y = divide(high - low, low)
   end function

   pure elemental function adoscillator_ta(open, high, low, close) result(y)
      real(dp), intent(in) :: open !! Open price for the same observation; missing values propagate.
      real(dp), intent(in) :: high !! High price for the same observation; missing values propagate.
      real(dp), intent(in) :: low !! Low price for the same observation; missing values propagate.
      real(dp), intent(in) :: close !! Close price for the same observation; missing values propagate.
      real(dp) :: y
      y = divide(high - open + close - low, high - low) * 50.0_dp
   end function

   pure function mom_ta(x, lag) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels; differences are zero-padded at the start.
      integer, intent(in), optional :: lag !! Positive differencing lag; default 25, at most size(x).
      real(dp) :: y(size(x))
      integer :: n
      n = 25
      if (present(lag)) n = lag
      y = nan_value()
      if (n < 1 .or. n > size(x)) return
      y = ttr_momentum(x, n)
      y(:n) = 0.0_dp
   end function

   pure function roc_ta(x, lag) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels; denominator is current level, not previous level.
      integer, intent(in), optional :: lag !! Positive differencing lag; default five.
      real(dp) :: y(size(x))
      integer :: n
      n = 5
      if (present(lag)) n = lag
      y = divide(mom_ta(x, n), x)
   end function

   pure function bias_ta(x, lag) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels compared to their startup-seeded EMA.
      integer, intent(in), optional :: lag !! EMA period; default five.
      real(dp) :: y(size(x)), average(size(x))
      integer :: n
      n = 5
      if (present(lag)) n = lag
      average = ema_ta(x, real(n, dp))
      y = divide(x - average, average)
   end function

   pure function osc_ta(x, lag1, lag2) result(y)
      real(dp), intent(in) :: x(:) !! Chronological price or volume observations.
      integer, intent(in), optional :: lag1 !! First EMA period; default 25.
      integer, intent(in), optional :: lag2 !! Second EMA period and denominator; default 65.
      real(dp) :: y(size(x)), average(size(x))
      integer :: n1, n2
      n1 = 25
      n2 = 65
      if (present(lag1)) n1 = lag1
      if (present(lag2)) n2 = lag2
      average = ema_ta(x, real(n2, dp))
      y = divide(ema_ta(x, real(n1, dp)) - average, average)
   end function

   pure function macd_ta(x, lag1, lag2) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels used for the two EMA series.
      integer, intent(in), optional :: lag1 !! Fast EMA period; default 12.
      integer, intent(in), optional :: lag2 !! Slow EMA period; default 26.
      real(dp) :: y(size(x))
      integer :: n1, n2
      n1 = 12
      n2 = 26
      if (present(lag1)) n1 = lag1
      if (present(lag2)) n2 = lag2
      y = ema_ta(x, real(n1, dp)) - ema_ta(x, real(n2, dp))
   end function

   pure function cds_ta(x, lag1, lag2, lag3) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels used for the MACD signal line.
      integer, intent(in), optional :: lag1 !! Fast MACD EMA period; default 12.
      integer, intent(in), optional :: lag2 !! Slow MACD EMA period; default 26.
      integer, intent(in), optional :: lag3 !! Signal smoothing EMA period; default nine.
      real(dp) :: y(size(x))
      integer :: n3
      n3 = 9
      if (present(lag3)) n3 = lag3
      y = ema_ta(macd_ta(x, lag1, lag2), real(n3, dp))
   end function

   pure function cdo_ta(x, lag1, lag2, lag3) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels; signal calculation preserves upstream default-period quirk.
      integer, intent(in), optional :: lag1 !! Fast period for first MACD only; default 12.
      integer, intent(in), optional :: lag2 !! Slow period for first MACD only; default 26.
      integer, intent(in), optional :: lag3 !! Signal smoothing period; default nine.
      real(dp) :: y(size(x))
      y = macd_ta(x, lag1, lag2) - cds_ta(x, lag3 = lag3)
   end function

   pure function fpk_ta(close, high, low, lag) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closing prices.
      real(dp), intent(in) :: high(:) !! Aligned highs; length must match close.
      real(dp), intent(in) :: low(:) !! Aligned lows; length must match close.
      integer, intent(in), optional :: lag !! Number of prior prices, so full window has lag+1 entries; default five.
      real(dp) :: y(size(close)), lowest, highest
      integer :: n, i, j
      y = nan_value()
      if (size(high) /= size(close) .or. size(low) /= size(close)) return
      n = 5
      if (present(lag)) n = lag
      if (n < 1) return
      do i = 1, size(close)
         j = max(1, i - n)
         if (any(ieee_is_nan(high(j:i))) .or. any(ieee_is_nan(low(j:i)))) cycle
         lowest = minval(low(j:i))
         highest = maxval(high(j:i))
         y(i) = divide(close(i) - lowest, highest - lowest)
      end do
   end function

   pure function wpr_ta(close, high, low, lag) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closes; result exactly matches fpkTA, as upstream.
      real(dp), intent(in) :: high(:) !! Aligned highs, same length as close.
      real(dp), intent(in) :: low(:) !! Aligned lows, same length as close.
      integer, intent(in), optional :: lag !! Prior-observation count; default 50.
      real(dp) :: y(size(close))
      integer :: n
      n = 50
      if (present(lag)) n = lag
      y = fpk_ta(close, high, low, n)
   end function

   pure function fpd_ta(close, high, low, lag1, lag2) result(y)
      real(dp), intent(in) :: close(:) !! Chronological close prices, aligned to the other price vectors.
      real(dp), intent(in) :: high(:) !! Chronological high prices, aligned to the other price vectors.
      real(dp), intent(in) :: low(:) !! Chronological low prices, aligned to the other price vectors.
      integer, intent(in), optional :: lag1 !! Stage 1 window/EMA period; default 5.
      integer, intent(in), optional :: lag2 !! Stage 2 window/EMA period; default 3.
      real(dp) :: y(size(close))
      integer :: n
      n = 3
      if (present(lag2)) n = lag2
      y = ema_ta(fpk_ta(close, high, low, lag1), real(n, dp))
   end function

   pure function spd_ta(close, high, low, lag1, lag2, lag3) result(y)
      real(dp), intent(in) :: close(:) !! Chronological close prices, aligned to the other price vectors.
      real(dp), intent(in) :: high(:) !! Chronological high prices, aligned to the other price vectors.
      real(dp), intent(in) :: low(:) !! Chronological low prices, aligned to the other price vectors.
      integer, intent(in), optional :: lag1 !! Stage 1 window/EMA period; default 5.
      integer, intent(in), optional :: lag2 !! Stage 2 window/EMA period; default 3.
      integer, intent(in), optional :: lag3 !! Stage 3 window/EMA period; default 9.
      real(dp) :: y(size(close))
      integer :: n
      n = 9
      if (present(lag3)) n = lag3
      y = ema_ta(fpd_ta(close, high, low, lag1, lag2), real(n, dp))
   end function

   pure function apd_ta(close, high, low, lag1, lag2, lag3, lag4) result(y)
      real(dp), intent(in) :: close(:) !! Chronological close prices, aligned to the other price vectors.
      real(dp), intent(in) :: high(:) !! Chronological high prices, aligned to the other price vectors.
      real(dp), intent(in) :: low(:) !! Chronological low prices, aligned to the other price vectors.
      integer, intent(in), optional :: lag1 !! Stage 1 window/EMA period; default 5.
      integer, intent(in), optional :: lag2 !! Stage 2 window/EMA period; default 3.
      integer, intent(in), optional :: lag3 !! Stage 3 window/EMA period; default 9.
      integer, intent(in), optional :: lag4 !! Stage 4 window/EMA period; default 9.
      real(dp) :: y(size(close))
      integer :: n
      n = 9
      if (present(lag4)) n = lag4
      y = ema_ta(spd_ta(close, high, low, lag1, lag2, lag3), real(n, dp))
   end function

   pure function stochastic_ta(close, high, low, lag1, lag2, lag3, slow) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closing prices.
      real(dp), intent(in) :: high(:) !! High prices aligned to close.
      real(dp), intent(in) :: low(:) !! Low prices aligned to close.
      integer, intent(in), optional :: lag1 !! Prior-observation count for raw K; default five.
      integer, intent(in), optional :: lag2 !! EMA period for raw D; default three.
      integer, intent(in), optional :: lag3 !! EMA period used for both slow columns; default five.
      logical, intent(in), optional :: slow !! True selects slow K and D; default false (fast).
      real(dp) :: y(size(close), 2)
      integer :: n
      y(:, 1) = fpk_ta(close, high, low, lag1)
      y(:, 2) = fpd_ta(close, high, low, lag1, lag2)
      n = 5
      if (present(lag3)) n = lag3
      if (present(slow)) then
         if (slow) then
            y(:, 1) = ema_ta(y(:, 1), real(n, dp))
            y(:, 2) = ema_ta(y(:, 2), real(n, dp))
         end if
      end if
   end function

   pure function rsi_ta(close, lag) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closes; ratio lies in [0,1] when defined.
      integer, intent(in), optional :: lag !! Prior-change count, full width lag+1; default 14.
      real(dp) :: y(size(close)), change(size(close)), gain(size(close))
      integer :: i, j, n
      y = nan_value()
      n = 14
      if (present(lag)) n = lag
      if (n < 1 .or. size(close) < 2) return
      change = mom_ta(close, 1)
      gain = change
      where (gain < 0.0_dp) gain = 0.0_dp
      do i = 1, size(close)
         j = max(1, i - n)
         y(i) = divide(sum(gain(j:i)), sum(abs(change(j:i))))
      end do
      y(1) = y(2)
   end function

   pure function accel_ta(x, n) result(y)
      real(dp), intent(in) :: x(:) !! Chronological levels; first n+1 outputs are unavailable.
      integer, intent(in), optional :: n !! Positive momentum lag; default 12.
      real(dp) :: y(size(x)), changes(size(x))
      integer :: width
      width = 12
      if (present(n)) width = n
      y = nan_value()
      if (width < 1 .or. width >= size(x)) return
      changes = mom_ta(x, width)
      y(width + 2:) = changes(width + 2:) - changes(width + 1:size(x) - 1)
   end function

   pure function adi_ta(high, low, close, volume) result(y)
      real(dp), intent(in) :: high(:) !! Chronological high prices.
      real(dp), intent(in) :: low(:) !! Low prices aligned to high.
      real(dp), intent(in) :: close(:) !! Closing prices aligned to high.
      real(dp), intent(in) :: volume(:) !! Traded volume aligned to high.
      real(dp) :: y(size(high)), running
      integer :: i
      y = nan_value()
      if (size(low) /= size(high) .or. size(close) /= size(high) .or. size(volume) /= size(high)) return
      running = 0.0_dp
      do i = 1, size(high)
         running = running + divide(2.0_dp * close(i) - high(i) - low(i), high(i) - low(i)) * volume(i)
         y(i) = running
      end do
   end function

   pure function bollinger_ta(x, lag, n_sd) result(y)
      real(dp), intent(in) :: x(:) !! Chronological prices; output columns are upper, original x, lower.
      integer, intent(in), optional :: lag !! Complete-window width; default 20.
      real(dp), intent(in), optional :: n_sd !! Sample-standard-deviation multiplier; default two.
      real(dp) :: y(size(x), 3), multiplier
      real(dp), allocatable :: mean(:), deviation(:)
      integer :: n
      y = nan_value()
      y(:, 2) = x
      n = 20
      if (present(lag)) n = lag
      if (n < 1 .or. n > size(x)) return
      multiplier = 2.0_dp
      if (present(n_sd)) multiplier = n_sd
      mean = sma(x, n)
      deviation = multiplier * sqrt(roll_var(x, n))
      y(n:, 1) = mean + deviation
      y(n:, 3) = mean - deviation
   end function

   pure function chaikino_ta(high, low, close, volume, lag1, lag2) result(y)
      real(dp), intent(in) :: high(:) !! Chronological high prices.
      real(dp), intent(in) :: low(:) !! Low prices aligned to high.
      real(dp), intent(in) :: close(:) !! Closing prices aligned to high.
      real(dp), intent(in) :: volume(:) !! Traded volume aligned to high.
      integer, intent(in), optional :: lag1 !! First EMA period; default 10.
      integer, intent(in), optional :: lag2 !! Second EMA period; default three.
      real(dp) :: y(size(high)), accumulation(size(high))
      integer :: n1, n2
      n1 = 10
      n2 = 3
      if (present(lag1)) n1 = lag1
      if (present(lag2)) n2 = lag2
      y = nan_value()
      if (n1 < 1 .or. n2 < 1) return
      accumulation = adi_ta(high, low, close, volume)
      y = ema_ta(accumulation, 2.0_dp / (real(n1, dp) + 1.0_dp)) - &
         ema_ta(accumulation, 2.0_dp / (real(n2, dp) + 1.0_dp))
   end function

   pure function chaikinv_ta(high, low, lag1, lag2) result(y)
      real(dp), intent(in) :: high(:) !! Chronological high prices.
      real(dp), intent(in) :: low(:) !! Aligned low prices, same length as high.
      integer, intent(in), optional :: lag1 !! Range EMA period; default 10.
      integer, intent(in), optional :: lag2 !! Positive change lag, at most series length; default 10.
      real(dp) :: y(size(high)), ranges(size(high))
      integer :: n1, n2
      y = nan_value()
      n1 = 10
      n2 = 10
      if (present(lag1)) n1 = lag1
      if (present(lag2)) n2 = lag2
      if (size(low) /= size(high)) return
      if (n1 < 1 .or. n2 < 1 .or. n2 > size(high)) return
      ranges = ema_ta(high - low, 2.0_dp / (real(n1, dp) + 1.0_dp))
      y(n2 + 1:) = (divide(ranges(n2 + 1:), ranges(:size(high) - n2)) - 1.0_dp) * 100.0_dp
   end function

   pure function garmanklass_ta(open, high, low, close) result(y)
      real(dp), intent(in) :: open(:) !! Strictly positive opening prices; intraday differences remain in price units.
      real(dp), intent(in) :: high(:) !! Aligned highs; upstream raw-price intraday formula is retained.
      real(dp), intent(in) :: low(:) !! Aligned lows; first output is unavailable.
      real(dp), intent(in) :: close(:) !! Aligned strictly positive closes used in overnight log returns.
      real(dp) :: y(size(open)), u, d, cc, overnight, variance
      integer :: i
      y = nan_value()
      if (size(high) /= size(open) .or. size(low) /= size(open) .or. size(close) /= size(open)) return
      do i = 2, size(open)
         if (open(i) <= 0.0_dp .or. close(i - 1) <= 0.0_dp) cycle
         u = high(i) - open(i)
         d = low(i) - open(i)
         cc = close(i) - open(i)
         overnight = (log(open(i)) - log(close(i - 1))) ** 2
         variance = 0.511_dp * (u - d) ** 2 - 0.019_dp * (cc * (u + d) - 2.0_dp * u * d) - 0.383_dp * cc ** 2
         variance = 0.88_dp * variance / 0.808_dp + 0.12_dp * overnight / 0.192_dp
         if (variance >= 0.0_dp) y(i) = sqrt(variance)
      end do
   end function

   pure function volume_index(close, volume, direction) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closing levels.
      real(dp), intent(in) :: volume(:) !! Aligned volumes used to select increasing or decreasing observations.
      integer, intent(in) :: direction !! Negative for falling-volume index; positive for rising-volume index.
      real(dp) :: y(size(close)), changes(size(close)), running
      integer :: i
      y = nan_value()
      if (size(volume) /= size(close) .or. size(close) == 0) return
      changes = roc_ta(close, 1) / 100.0_dp
      running = changes(1) * 0.0_dp
      y(1) = running
      do i = 2, size(close)
         if (real(direction, dp) * (volume(i) - volume(i - 1)) > 0.0_dp) then
            running = running + changes(i)
         else
            running = running + changes(i) * 0.0_dp
         end if
         y(i) = running
      end do
   end function

   pure function nvi_ta(close, volume) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closes; current-price-denominator changes are divided by 100.
      real(dp), intent(in) :: volume(:) !! Aligned volumes selecting the relevant signed volume changes.
      real(dp) :: y(size(close))
      y = volume_index(close, volume, -1)
   end function

   pure function pvi_ta(close, volume) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closes; current-price-denominator changes are divided by 100.
      real(dp), intent(in) :: volume(:) !! Aligned volumes selecting the relevant signed volume changes.
      real(dp) :: y(size(close))
      y = volume_index(close, volume, 1)
   end function

   pure function obv_ta(close, volume) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closes; exact equal prices contribute zero.
      real(dp), intent(in) :: volume(:) !! Aligned volume; first cumulative observation is zero.
      real(dp) :: y(size(close)), running, difference
      integer :: i
      y = nan_value()
      if (size(volume) /= size(close) .or. size(close) == 0) return
      running = 0.0_dp * volume(1)
      y(1) = running
      do i = 2, size(close)
         difference = close(i) - close(i - 1)
         if (ieee_is_nan(difference)) then
            running = nan_value()
         else if (difference == 0.0_dp) then
            running = running + 0.0_dp * volume(i)
         else
            running = running + volume(i) * sign(1.0_dp, difference)
         end if
         y(i) = running
      end do
   end function

   pure function pvtrend_ta(close, volume) result(y)
      real(dp), intent(in) :: close(:) !! Chronological closes; percentage changes use previous close as denominator.
      real(dp), intent(in) :: volume(:) !! Aligned volumes multiplying percentage changes.
      real(dp) :: y(size(close)), running
      integer :: i
      y = nan_value()
      if (size(volume) /= size(close) .or. size(close) == 0) return
      running = 0.0_dp * volume(1)
      y(1) = running
      do i = 2, size(close)
         running = running + volume(i) * (divide(close(i), close(i - 1)) - 1.0_dp) * 100.0_dp
         y(i) = running
      end do
   end function

   pure function williamsad_ta(high, low, close) result(y)
      real(dp), intent(in) :: high(:) !! Chronological highs.
      real(dp), intent(in) :: low(:) !! Aligned lows; same length as high.
      real(dp), intent(in) :: close(:) !! Aligned closes; sign of change selects accumulation contribution.
      real(dp) :: y(size(high)), running, change
      integer :: i
      y = nan_value()
      if (size(low) /= size(high) .or. size(close) /= size(high) .or. size(high) == 0) return
      running = 0.0_dp
      y(1) = running
      do i = 2, size(high)
         change = close(i) - close(i - 1)
         if (ieee_is_nan(change)) then
            running = nan_value()
         else if (change > 0.0_dp) then
            running = running + close(i) - low(i)
         else if (change < 0.0_dp) then
            running = running - (high(i) - close(i))
         end if
         y(i) = running
      end do
   end function

   pure function williamsr_ta(high, low, close, lag) result(y)
      real(dp), intent(in) :: high(:) !! Chronological highs; complete windows unlike fpkTA's startup convention.
      real(dp), intent(in) :: low(:) !! Aligned lows, same length as high.
      real(dp), intent(in) :: close(:) !! Aligned closes; output is percentage distance below rolling high.
      integer, intent(in), optional :: lag !! Complete-window width; default 20.
      real(dp) :: y(size(high))
      real(dp), allocatable :: highest(:), lowest(:)
      integer :: n
      y = nan_value()
      n = 20
      if (present(lag)) n = lag
      if (size(low) /= size(high) .or. size(close) /= size(high)) return
      if (n < 1 .or. n > size(high)) return
      highest = roll_fun(high, n, 'max', trim = .false.)
      lowest = roll_fun(low, n, 'min', trim = .false.)
      y = divide(highest - close, highest - lowest) * 100.0_dp
   end function

   pure function daily_ta(ohlc, indicator, lags, select_column) result(y)
      real(dp), intent(in) :: ohlc(:, :) !! Numeric columns ordered Open, High, Low, Close; at least four columns.
      character(len=*), intent(in), optional :: indicator !! Source .dailyTA selector, e.g. ema, macd, rsi; default ema.
      integer, intent(in), optional :: lags(:) !! Required periods: one normally, two for osc/macd/fpd, three cds/cdo/spd, four apd.
      integer, intent(in), optional :: select_column !! Input column for univariate indicators; default four (Close).
      real(dp) :: y(size(ohlc, 1))
      integer :: column, periods(4), supplied
      character(len=:), allocatable :: which
      y = nan_value()
      if (size(ohlc, 2) < 4) return
      column = 4
      if (present(select_column)) column = select_column
      if (column < 1 .or. column > size(ohlc, 2)) return
      periods = 9
      supplied = 1
      if (present(lags)) then
         supplied = min(4, size(lags))
         periods(:supplied) = lags(:supplied)
      end if
      if (supplied < 1) return
      which = 'ema'
      if (present(indicator)) which = indicator
      select case(which)
      case('osc', 'macd', 'fpd')
         if (supplied < 2) return
      case('cds', 'cdo', 'spd')
         if (supplied < 3) return
      case('apd')
         if (supplied < 4) return
      end select
      select case(which)
      case('ema')
         y = ema_ta(ohlc(:, column), real(periods(1), dp))
      case('bias')
         y = bias_ta(ohlc(:, column), periods(1))
      case('medprice')
         y = medprice_ta(ohlc(:, 2), ohlc(:, 3))
      case('typicalprice')
         y = typicalprice_ta(ohlc(:, 2), ohlc(:, 3), ohlc(:, 4))
      case('wclose')
         y = wclose_ta(ohlc(:, 2), ohlc(:, 3), ohlc(:, 4))
      case('roc')
         y = roc_ta(ohlc(:, column), periods(1))
      case('mom')
         y = mom_ta(ohlc(:, column), periods(1))
      case('osc')
         y = osc_ta(ohlc(:, column), periods(1), periods(2))
      case('macd')
         y = macd_ta(ohlc(:, column), periods(1), periods(2))
      case('cds')
         y = cds_ta(ohlc(:, column), periods(1), periods(2), periods(3))
      case('cdo')
         y = cdo_ta(ohlc(:, column), periods(1), periods(2), periods(3))
      case('vohl')
         y = vohl_ta(ohlc(:, 2), ohlc(:, 3))
      case('vor')
         y = vor_ta(ohlc(:, 2), ohlc(:, 3))
      case('rsi')
         y = rsi_ta(ohlc(:, 4), periods(1))
      case('fpk')
         y = fpk_ta(ohlc(:, 4), ohlc(:, 2), ohlc(:, 3), periods(1))
      case('fpd')
         y = fpd_ta(ohlc(:, 4), ohlc(:, 2), ohlc(:, 3), periods(1), periods(2))
      case('spd')
         y = spd_ta(ohlc(:, 4), ohlc(:, 2), ohlc(:, 3), periods(1), periods(2), periods(3))
      case('apd')
         y = apd_ta(ohlc(:, 4), ohlc(:, 2), ohlc(:, 3), periods(1), periods(2), periods(3), periods(4))
      case('wpr')
         y = wpr_ta(ohlc(:, 4), ohlc(:, 2), ohlc(:, 3), periods(1))
      end select
   end function
end module ftrading_indicators
