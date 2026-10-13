! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from TTR 0.24.4, Joshua M. Ulrich and contributors; see NOTICE.md.
module ttr_indicators
   use ttr_base
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: tr, atr, adx, bbands, cci, clv, cmf, cmo, dpo, dvi, donchian_channel
   public :: emv, gmma, kst, macd, mfi, obv, pbands, rsi, smi, snr, tdi, trix, vhf, wpr
   public :: chaikin_ad, chaikin_volatility, keltner_channels, stoch, ultimate_oscillator, williams_ad
contains
   pure function tr(hlc) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Chronological high, low, close matrix with three columns.
      real(dp) :: z(size(hlc, 1), 3)
      integer :: i
      z = nan_value()
      if (size(hlc, 2) /= 3) return
      do i = 2, size(hlc, 1)
         if (any(ieee_is_nan(hlc(i, :2))) .or. ieee_is_nan(hlc(i - 1, 3))) cycle
         z(i, 2) = max(hlc(i, 1), hlc(i - 1, 3))
         z(i, 3) = min(hlc(i, 2), hlc(i - 1, 3))
         z(i, 1) = z(i, 2) - z(i, 3)
      end do
   end function

   pure function atr(hlc, n, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      integer, intent(in) :: n !! Positive average true range smoothing period.
      character(*), intent(in), optional :: ma_type !! Moving-average selector; default Wilder.
      real(dp) :: z(size(hlc, 1), 4), ranges(size(hlc, 1), 3)
      ranges = tr(hlc)
      z(:, 1) = ranges(:, 1)
      z(:, 2) = moving_average(ranges(:, 1), n, ma_type, 'Wilder')
      z(:, 3:4) = ranges(:, 2:3)
   end function

   pure function adx(hlc, n, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      integer, intent(in) :: n !! Positive directional movement smoothing period.
      character(*), intent(in), optional :: ma_type !! Final DX smoothing; default Wilder.
      real(dp) :: z(size(hlc, 1), 4), dh(size(hlc, 1)), dl(size(hlc, 1))
      real(dp) :: plus(size(hlc, 1)), minus(size(hlc, 1)), ranges(size(hlc, 1), 3), den(size(hlc, 1))
      integer :: i
      z = nan_value()
      if (size(hlc, 2) /= 3) return
      dh = momentum(hlc(:, 1), 1)
      dl = -momentum(hlc(:, 2), 1)
      plus = nan_value()
      minus = nan_value()
      do i = 1, size(dh)
         if (ieee_is_nan(dh(i)) .or. ieee_is_nan(dl(i))) cycle
         plus(i) = 0.0_dp
         minus(i) = 0.0_dp
         if (dh(i) > dl(i) .and. dh(i) > 0.0_dp) plus(i) = dh(i)
         if (dl(i) > dh(i) .and. dl(i) > 0.0_dp) minus(i) = dl(i)
      end do
      ranges = tr(hlc)
      den = wilder_sum(ranges(:, 1), n)
      z(:, 1) = 100.0_dp * divide(wilder_sum(plus, n), den)
      z(:, 2) = 100.0_dp * divide(wilder_sum(minus, n), den)
      z(:, 3) = 100.0_dp * divide(abs(z(:, 1) - z(:, 2)), z(:, 1) + z(:, 2))
      z(:, 4) = moving_average(z(:, 3), n, ma_type, 'Wilder')
   end function

   pure function bbands(hlc, n, sd, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      integer, intent(in) :: n !! Positive band lookback window.
      real(dp), intent(in), optional :: sd !! Population standard-deviation multiplier, default 2.
      character(*), intent(in), optional :: ma_type !! Center smoothing selector, default SMA.
      real(dp) :: z(size(hlc, 1), 4), p(size(hlc, 1)), spread(size(hlc, 1)), mult
      mult = 2.0_dp
      if (present(sd)) mult = sd
      p = prices(hlc)
      spread = mult * run_sd(p, n, sample = .false.)
      z(:, 2) = moving_average(p, n, ma_type)
      z(:, 1) = z(:, 2) - spread
      z(:, 3) = z(:, 2) + spread
      z(:, 4) = divide(p - z(:, 1), z(:, 3) - z(:, 1))
   end function

   pure function cci(hlc, n, constant, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      integer, intent(in) :: n !! Positive commodity channel index window.
      real(dp), intent(in), optional :: constant !! Mean-deviation scaling constant, default 0.015.
      character(*), intent(in), optional :: ma_type !! Center smoothing selector, default SMA.
      real(dp) :: z(size(hlc, 1)), p(size(hlc, 1)), avg(size(hlc, 1)), scale
      scale = 0.015_dp
      if (present(constant)) scale = constant
      p = prices(hlc)
      avg = moving_average(p, n, ma_type)
      z = divide(p - avg, scale * run_mad(p, n, center = avg, median = .false.))
   end function

   pure function clv(hlc) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      real(dp) :: z(size(hlc, 1))
      integer :: i
      z = nan_value()
      if (size(hlc, 2) /= 3) return
      do i = 1, size(z)
         if (any(ieee_is_nan(hlc(i, :)))) cycle
         z(i) = divide(2.0_dp * hlc(i, 3) - hlc(i, 2) - hlc(i, 1), hlc(i, 1) - hlc(i, 2))
         if (.not. ieee_is_finite(z(i))) z(i) = 0.0_dp
      end do
   end function

   pure function cmf(hlc, volume, n) result(z)
      real(dp), intent(in) :: hlc(:, :) !! High, low, close matrix, one row per time.
      real(dp), intent(in) :: volume(:) !! Volume aligned with every price row.
      integer, intent(in) :: n !! Positive Chaikin money flow lookback.
      real(dp) :: z(size(hlc, 1))
      z = nan_value()
      if (size(volume) /= size(z)) return
      z = divide(run_sum(clv(hlc) * volume, n), run_sum(volume, n))
   end function

   pure subroutine up_down(x, up, down)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      real(dp), intent(out) :: up(size(x)) !! Positive one-step changes, zero for declines.
      real(dp), intent(out) :: down(size(x)) !! Absolute negative changes, zero for rises.
      real(dp) :: change(size(x))
      integer :: i
      change = momentum(x, 1)
      up = nan_value()
      down = nan_value()
      do i = 1, size(x)
         if (ieee_is_nan(change(i))) cycle
         up(i) = max(0.0_dp, change(i))
         down(i) = max(0.0_dp, -change(i))
      end do
   end subroutine

   pure function cmo(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n !! Positive Chande momentum window.
      real(dp) :: z(size(x)), up(size(x)), down(size(x))
      call up_down(x, up, down)
      up = run_sum(up, n)
      down = run_sum(down, n)
      z = 100.0_dp * divide(up - down, up + down)
   end function

   pure function rsi(x, n, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n !! Positive relative-strength smoothing period.
      character(*), intent(in), optional :: ma_type !! Smoothing for gains and losses; default Wilder.
      real(dp) :: z(size(x)), up(size(x)), down(size(x))
      call up_down(x, up, down)
      up = moving_average(up, n, ma_type, 'Wilder')
      down = moving_average(down, n, ma_type, 'Wilder')
      z = 100.0_dp * divide(up, up + down)
   end function

   pure function dpo(x, n, shift_by, percent, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n !! Positive detrended-price smoothing period.
      integer, intent(in), optional :: shift_by !! Forward shift of mean; default trunc(n/2+1).
      logical, intent(in), optional :: percent !! Percentage price-to-mean deviation; default false.
      character(*), intent(in), optional :: ma_type !! Center smoothing, default SMA.
      real(dp) :: z(size(x)), avg(size(x))
      integer :: lead
      lead = n / 2 + 1
      if (present(shift_by)) lead = shift_by
      avg = shift(moving_average(x, n, ma_type), -lead)
      z = x - avg
      if (present(percent)) then
         if (percent) z = 100.0_dp * (divide(x, avg) - 1.0_dp)
      end if
   end function

   pure function dvi(x, n, weights, smooth, magnitude, stretch, exact_multiplier) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n !! Positive rolling percent-rank window.
      real(dp), intent(in), optional :: weights(2) !! Magnitude/stretch weights, normalized; default [0.8,0.2].
      integer, intent(in), optional :: smooth !! Initial price averaging window, default 3.
      integer, intent(in), optional :: magnitude(3) !! Magnitude windows, default [5,100,5].
      integer, intent(in), optional :: stretch(3) !! Stretch windows, default [10,100,2].
      real(dp), intent(in), optional :: exact_multiplier !! Percent-rank tie weight, default 1.
      real(dp) :: z(size(x), 3), r(size(x)), b(size(x)), mag(size(x)), str(size(x)), w(2), mult
      integer :: sm, mg(3), st(3), i
      sm = 3
      mg = [5, 100, 5]
      st = [10, 100, 2]
      w = [0.8_dp, 0.2_dp]
      mult = 1.0_dp
      if (present(smooth)) sm = smooth
      if (present(magnitude)) mg = magnitude
      if (present(stretch)) st = stretch
      if (present(weights)) w = weights
      if (present(exact_multiplier)) mult = exact_multiplier
      z = nan_value()
      if (sum(w) == 0.0_dp) return
      w = w / sum(w)
      r = divide(x, sma(x, sm)) - 1.0_dp
      mag = sma((sma(r, mg(1)) + sma(r, mg(2)) / 10.0_dp) / 2.0_dp, mg(3))
      b = nan_value()
      do i = 2, size(x)
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(x(i - 1))) cycle
         b(i) = -1.0_dp
         if (x(i) > x(i - 1)) b(i) = 1.0_dp
      end do
      str = sma((run_sum(b, st(1)) + run_sum(b, st(2)) / 10.0_dp) / 2.0_dp, st(3))
      z(:, 1) = run_percent_rank(mag, n, exact_multiplier = mult)
      z(:, 2) = run_percent_rank(str, n, exact_multiplier = mult)
      z(:, 3) = w(1) * z(:, 1) + w(2) * z(:, 2)
   end function

   pure function donchian_channel(hl, n, include_lag) result(z)
      real(dp), intent(in) :: hl(:, :) !! One close column or two high and low columns.
      integer, intent(in) :: n !! Positive rolling range window.
      logical, intent(in), optional :: include_lag !! Lag the finished channel one row; default false.
      real(dp) :: z(size(hl, 1), 3)
      integer :: j
      z = nan_value()
      if (size(hl, 2) < 1 .or. size(hl, 2) > 2) return
      z(:, 1) = run_max(hl(:, 1), n)
      z(:, 3) = run_min(hl(:, size(hl, 2)), n)
      z(:, 2) = (z(:, 1) + z(:, 3)) / 2.0_dp
      if (present(include_lag)) then
         if (include_lag) then
            do j = 1, 3
               z(:, j) = shift(z(:, j), 1)
            end do
         end if
      end if
   end function

   pure function emv(hl, volume, n, vol_divisor, ma_type) result(z)
      real(dp), intent(in) :: hl(:, :) !! Two columns: high and low; rows are times.
      real(dp), intent(in) :: volume(:) !! Volumes aligned with the price rows.
      integer, intent(in) :: n !! Positive smoothing window for ease of movement.
      real(dp), intent(in), optional :: vol_divisor !! Volume scale divisor, default 10000.
      character(*), intent(in), optional :: ma_type !! Smoothing selector, default SMA.
      real(dp) :: z(size(hl, 1), 2), scale
      z = nan_value()
      if (size(hl, 2) /= 2 .or. size(volume) /= size(hl, 1)) return
      scale = 10000.0_dp
      if (present(vol_divisor)) scale = vol_divisor
      z(:, 1) = divide(momentum((hl(:, 1) + hl(:, 2)) / 2.0_dp, 1), &
         divide(divide(volume, scale), hl(:, 1) - hl(:, 2)))
      z(:, 2) = moving_average(z(:, 1), n, ma_type)
   end function

   pure function gmma(x, periods, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: periods(:) !! Concatenated short and long windows in desired column order.
      character(*), intent(in), optional :: ma_type !! Smoothing selector, default EMA.
      real(dp) :: z(size(x), size(periods))
      integer :: j
      do j = 1, size(periods)
         z(:, j) = moving_average(x, periods(j), ma_type, 'EMA')
      end do
   end function

   pure function kst(x, n, n_roc, n_sig, weights, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n(:) !! Smoothing window for each rate-of-change component.
      integer, intent(in) :: n_roc(:) !! Return lags, one per smoothing window.
      integer, intent(in) :: n_sig !! Positive final signal smoothing window.
      real(dp), intent(in), optional :: weights(:) !! Component weights; default 1:size(n).
      character(*), intent(in), optional :: ma_type !! Smoothing selector for all stages, default SMA.
      real(dp) :: z(size(x), 2), w(size(n))
      integer :: j
      z = nan_value()
      if (size(n) /= size(n_roc) .or. size(n) == 0) return
      w = [(real(j, dp), j = 1, size(n))]
      if (present(weights)) then
         if (size(weights) /= size(n)) return
         w = weights
      end if
      z(:, 1) = 0.0_dp
      do j = 1, size(n)
         z(:, 1) = z(:, 1) + w(j) * moving_average(roc(x, n_roc(j)), n(j), ma_type)
      end do
      z(:, 1) = 100.0_dp * z(:, 1)
      z(:, 2) = moving_average(z(:, 1), n_sig, ma_type)
   end function

   pure function macd(x, n_fast, n_slow, n_sig, percent, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n_fast !! Positive fast smoothing period.
      integer, intent(in) :: n_slow !! Positive slow smoothing period.
      integer, intent(in) :: n_sig !! Positive signal smoothing period.
      logical, intent(in), optional :: percent !! Percentage ratio (default) rather than difference.
      character(*), intent(in), optional :: ma_type !! Common smoothing selector, default EMA.
      real(dp) :: z(size(x), 2), fast(size(x)), slow(size(x))
      logical :: pct
      pct = .true.
      if (present(percent)) pct = percent
      fast = moving_average(x, n_fast, ma_type, 'EMA')
      slow = moving_average(x, n_slow, ma_type, 'EMA')
      z(:, 1) = fast - slow
      if (pct) z(:, 1) = 100.0_dp * (divide(fast, slow) - 1.0_dp)
      z(:, 2) = moving_average(z(:, 1), n_sig, ma_type, 'EMA')
   end function

   pure function mfi(hlc, volume, n) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      real(dp), intent(in) :: volume(:) !! Volumes aligned with all price rows.
      integer, intent(in) :: n !! Positive money flow summation period.
      real(dp) :: z(size(hlc, 1)), p(size(hlc, 1)), up(size(hlc, 1)), down(size(hlc, 1))
      integer :: i
      z = nan_value()
      if (size(volume) /= size(z)) return
      p = prices(hlc)
      up = nan_value()
      down = nan_value()
      do i = 2, size(z)
         if (ieee_is_nan(p(i)) .or. ieee_is_nan(p(i - 1))) cycle
         up(i) = 0.0_dp
         down(i) = 0.0_dp
         if (p(i) > p(i - 1)) up(i) = p(i) * volume(i)
         if (p(i) < p(i - 1)) down(i) = p(i) * volume(i)
      end do
      up = run_sum(up, n)
      down = run_sum(down, n)
      z = 100.0_dp - divide(100.0_dp, 1.0_dp + divide(up, down))
      where (down == 0.0_dp) z = 100.0_dp
      where (down == 0.0_dp .and. up == 0.0_dp) z = 50.0_dp
   end function

   pure function obv(x, volume) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices with leading NaNs permitted.
      real(dp), intent(in) :: volume(:) !! Volumes aligned with every price.
      real(dp) :: z(size(x)), total, changes(size(x))
      integer :: i
      z = nan_value()
      if (size(volume) /= size(x) .or. size(x) == 0) return
      changes = roc(x, 1)
      total = volume(1)
      z(1) = total
      do i = 2, size(x)
         if (ieee_is_nan(changes(i))) then
            total = nan_value()
         else if (abs(changes(i)) < sqrt(epsilon(1.0_dp))) then
            ! Upstream treats sufficiently small logarithmic changes as ties.
         else if (changes(i) > 0.0_dp) then
            total = total + volume(i)
         else
            total = total - volume(i)
         end if
         z(i) = total
      end do
   end function

   pure function pbands(x, n, sd, fast_n, centered, long_average, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological univariate prices.
      integer, intent(in) :: n !! Positive main price-band window.
      real(dp), intent(in), optional :: sd !! Standard-deviation multiplier, default 2.
      integer, intent(in), optional :: fast_n !! Fast moving-average window, default 2.
      logical, intent(in), optional :: centered !! Add smoothed standardized difference to center; default false.
      logical, intent(in), optional :: long_average !! Use twice n for center correction; default false.
      character(*), intent(in), optional :: ma_type !! Common smoothing selector, default SMA.
      real(dp) :: z(size(x), 3), avg(size(x)), spread(size(x)), diff(size(x)), mult
      integer :: nf, nc
      logical :: adjust
      nf = 2
      mult = 2.0_dp
      adjust = .false.
      if (present(fast_n)) nf = fast_n
      if (present(sd)) mult = sd
      if (present(centered)) adjust = centered
      avg = moving_average(x, n, ma_type)
      diff = avg - moving_average(x, nf, ma_type)
      spread = run_sd(diff, n, sample = .false.)
      z(:, 2) = avg
      if (adjust) then
         nc = n
         if (present(long_average)) then
            if (long_average) nc = 2 * n
         end if
         z(:, 2) = avg + moving_average(divide(diff, spread), nc, ma_type)
      end if
      z(:, 1) = z(:, 2) - mult * spread
      z(:, 3) = z(:, 2) + mult * spread
   end function

   pure function snr(hlc, n, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      integer, intent(in) :: n !! Positive signal/noise lag and ATR period.
      character(*), intent(in), optional :: ma_type !! ATR smoother, default Wilder.
      real(dp) :: z(size(hlc, 1)), a(size(hlc, 1), 4)
      z = nan_value()
      if (size(hlc, 2) /= 3) return
      a = atr(hlc, n, ma_type)
      z = divide(abs(momentum(hlc(:, 3), n)), a(:, 2))
   end function

   pure function tdi(x, n, multiple) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n !! Positive momentum lag and base sum period.
      integer, intent(in), optional :: multiple !! Positive longer-window multiplier, default 2.
      real(dp) :: z(size(x), 2), mom(size(x))
      integer :: m
      m = 2
      if (present(multiple)) m = multiple
      mom = momentum(x, n)
      where (ieee_is_nan(mom)) mom = 0.0_dp
      z(:, 2) = run_sum(mom, n)
      z(:, 1) = abs(z(:, 2)) - (run_sum(abs(mom), n * m) - run_sum(abs(mom), n))
   end function

   pure function trix(x, n, n_sig, percent, ma_type) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      integer, intent(in) :: n !! Positive period for each of the three smoothing stages.
      integer, intent(in) :: n_sig !! Positive final signal smoothing period.
      logical, intent(in), optional :: percent !! Percentage simple change (default), otherwise difference.
      character(*), intent(in), optional :: ma_type !! Common smoothing selector, default EMA.
      real(dp) :: z(size(x), 2), avg(size(x))
      logical :: pct
      pct = .true.
      if (present(percent)) pct = percent
      avg = moving_average(x, n, ma_type, 'EMA')
      avg = moving_average(avg, n, ma_type, 'EMA')
      avg = moving_average(avg, n, ma_type, 'EMA')
      z(:, 1) = momentum(avg, 1)
      if (pct) z(:, 1) = 100.0_dp * roc(avg, 1, .true.)
      z(:, 2) = moving_average(z(:, 1), n_sig, ma_type, 'EMA')
   end function

   pure function vhf(hlc, n) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      integer, intent(in) :: n !! Positive vertical-horizontal filter window.
      real(dp) :: z(size(hlc, 1)), h(size(hlc, 1)), l(size(hlc, 1)), c(size(hlc, 1))
      call high_low_close(hlc, h, l, c)
      z = divide(run_max(h, n) - run_min(l, n), run_sum(abs(momentum(c, 1)), n))
   end function

   pure function wpr(hlc, n) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      integer, intent(in) :: n !! Positive Williams percent-range window.
      real(dp) :: z(size(hlc, 1)), h(size(hlc, 1)), l(size(hlc, 1)), c(size(hlc, 1))
      call high_low_close(hlc, h, l, c)
      h = run_max(h, n)
      l = run_min(l, n)
      z = divide(h - c, h - l)
      where (h == l .and. h == c) z = 0.5_dp
   end function

   pure function chaikin_ad(hlc, volume) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      real(dp), intent(in) :: volume(:) !! Volumes aligned with all price rows.
      real(dp) :: z(size(hlc, 1))
      z = nan_value()
      if (size(volume) /= size(z)) return
      z = run_sum(clv(hlc) * volume, 1, cumulative = .true.)
   end function

   pure function chaikin_volatility(hl, n, ma_type) result(z)
      real(dp), intent(in) :: hl(:, :) !! Two columns: high and low; rows are times.
      integer, intent(in) :: n !! Positive smoothing window and change lag.
      character(*), intent(in), optional :: ma_type !! Range smoothing selector, default EMA.
      real(dp) :: z(size(hl, 1))
      z = nan_value()
      if (size(hl, 2) /= 2) return
      z = roc(moving_average(hl(:, 1) - hl(:, 2), n, ma_type, 'EMA'), n, .true.)
   end function

   pure function keltner_channels(hlc, n, atr_multiplier, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      integer, intent(in) :: n !! Positive mean and ATR period.
      real(dp), intent(in), optional :: atr_multiplier !! ATR band-width multiplier, default 2.
      character(*), intent(in), optional :: ma_type !! Center smoothing selector, default EMA.
      real(dp) :: z(size(hlc, 1), 3), a(size(hlc, 1), 4), mult
      mult = 2.0_dp
      if (present(atr_multiplier)) mult = atr_multiplier
      a = atr(hlc, n)
      z(:, 2) = moving_average(prices(hlc), n, ma_type, 'EMA')
      z(:, 1) = z(:, 2) - mult * a(:, 2)
      z(:, 3) = z(:, 2) + mult * a(:, 2)
   end function

   pure function stoch(hlc, n_fast_k, n_fast_d, n_slow_d, bounded, smooth, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      integer, intent(in) :: n_fast_k !! Positive high-low lookback window.
      integer, intent(in) :: n_fast_d !! Positive fast D smoothing window.
      integer, intent(in) :: n_slow_d !! Positive slow D smoothing window.
      logical, intent(in), optional :: bounded !! Include current high/low in range (default true).
      integer, intent(in), optional :: smooth !! Numerator/denominator smoothing window, default 1.
      character(*), intent(in), optional :: ma_type !! Common smoothing selector, default SMA.
      real(dp) :: z(size(hlc, 1), 3), h(size(hlc, 1)), l(size(hlc, 1)), c(size(hlc, 1))
      integer :: sm
      sm = 1
      if (present(smooth)) sm = smooth
      call high_low_close(hlc, h, l, c)
      if (present(bounded)) then
         if (.not. bounded .and. size(h) > 0) then
            h = [h(1), h(:size(h) - 1)]
            l = [l(1), l(:size(l) - 1)]
         end if
      end if
      h = run_max(h, n_fast_k)
      l = run_min(l, n_fast_k)
      z(:, 1) = divide(moving_average(c - l, sm, ma_type), moving_average(h - l, sm, ma_type))
      where (.not. ieee_is_finite(z(:, 1)) .and. .not. ieee_is_nan(z(:, 1))) z(:, 1) = 0.5_dp
      z(:, 2) = moving_average(z(:, 1), n_fast_d, ma_type)
      z(:, 3) = moving_average(z(:, 2), n_slow_d, ma_type)
   end function

   pure function smi(hlc, n, n_fast, n_slow, n_sig, bounded, ma_type) result(z)
      real(dp), intent(in) :: hlc(:, :) !! One close column or three high, low, close columns.
      integer, intent(in) :: n !! Positive high-low lookback window.
      integer, intent(in) :: n_fast !! Positive second-stage smoothing period.
      integer, intent(in) :: n_slow !! Positive first-stage smoothing period.
      integer, intent(in) :: n_sig !! Positive signal smoothing period.
      logical, intent(in), optional :: bounded !! Include current price in range (default true).
      character(*), intent(in), optional :: ma_type !! Common smoothing selector, default EMA.
      real(dp) :: z(size(hlc, 1), 2), h(size(hlc, 1)), l(size(hlc, 1)), c(size(hlc, 1))
      real(dp) :: hi(size(hlc, 1)), lo(size(hlc, 1)), num(size(hlc, 1)), den(size(hlc, 1))
      logical :: bound
      bound = .true.
      if (present(bounded)) bound = bounded
      call high_low_close(hlc, h, l, c)
      if (bound) then
         hi = run_max(h, n)
         lo = run_min(l, n)
      else
         if (size(h) == 0) return
         hi = run_max([h(1), h(:size(h) - 1)], n)
         ! TTR 0.24.4 deliberately uses runMax for low in the unbounded branch.
         lo = run_max([l(1), l(:size(l) - 1)], n)
      end if
      where (ieee_is_nan(hi)) hi = h
      where (ieee_is_nan(lo)) lo = l
      num = moving_average(c - (hi + lo) / 2.0_dp, n_slow, ma_type, 'EMA')
      den = moving_average(hi - lo, n_slow, ma_type, 'EMA')
      num = moving_average(num, n_fast, ma_type, 'EMA')
      den = moving_average(den, n_fast, ma_type, 'EMA')
      z(:, 1) = 200.0_dp * divide(num, den)
      z(:, 2) = moving_average(z(:, 1), n_sig, ma_type, 'EMA')
   end function

   pure function ultimate_oscillator(hlc, n, weights) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      integer, intent(in) :: n(3) !! Three positive buying-pressure and range windows.
      real(dp), intent(in), optional :: weights(3) !! Window weights, default [4,2,1].
      real(dp) :: z(size(hlc, 1)), a(size(hlc, 1), 3), bp(size(hlc, 1)), w(3)
      integer :: j
      z = nan_value()
      if (size(hlc, 2) /= 3) return
      w = [4.0_dp, 2.0_dp, 1.0_dp]
      if (present(weights)) w = weights
      if (sum(w) == 0.0_dp) return
      a = tr(hlc)
      bp = hlc(:, 3) - a(:, 3)
      z = 0.0_dp
      do j = 1, 3
         z = z + w(j) * divide(run_sum(bp, n(j)), run_sum(a(:, 1), n(j)))
      end do
      z = 100.0_dp * z / sum(w)
   end function

   pure function williams_ad(hlc) result(z)
      real(dp), intent(in) :: hlc(:, :) !! Three columns: high, low, close; rows are times.
      real(dp) :: z(size(hlc, 1)), a(size(hlc, 1), 3), change(size(hlc, 1))
      integer :: i
      z = nan_value()
      if (size(hlc, 2) /= 3) return
      a = tr(hlc)
      change = momentum(hlc(:, 3), 1)
      do i = 1, size(z)
         if (ieee_is_nan(change(i))) cycle
         if (change(i) > 0.0_dp) then
            z(i) = hlc(i, 3) - a(i, 3)
         else if (change(i) < 0.0_dp) then
            z(i) = hlc(i, 3) - a(i, 2)
         else
            z(i) = 0.0_dp
         end if
      end do
      z = run_sum(z, 1, cumulative = .true.)
   end function
end module ttr_indicators
