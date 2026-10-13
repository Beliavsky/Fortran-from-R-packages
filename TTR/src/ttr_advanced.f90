! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from TTR 0.24.4, Joshua M. Ulrich and contributors; see NOTICE.md.
module ttr_advanced
   use ttr_base, only: dp, nan_value, first_valid, divide, shift, roc, run_sum, &
      run_mean, run_sd, run_var, run_cov
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: sar, zigzag, adj_ratios, aroon, cti, volatility, roll_sfm

   interface sar
      module procedure sar_matrix, sar_vectors
   end interface
   interface zigzag
      module procedure zigzag_matrix, zigzag_vector
   end interface
   interface aroon
      module procedure aroon_matrix, aroon_vector
   end interface
   interface volatility
      module procedure volatility_matrix, volatility_vector
   end interface
contains
   pure function sar_matrix(hl, accel) result(z)
      real(dp), intent(in) :: hl(:, :) !! Chronological high/low matrix, with exactly two columns.
      real(dp), intent(in), optional :: accel(:) !! Initial/increment and maximum acceleration, default [0.02,0.2].
      real(dp) :: z(size(hl, 1))
      z = nan_value()
      if (size(hl, 2) /= 2) return
      z = sar_vectors(hl(:, 1), hl(:, 2), accel)
   end function

   pure function sar_vectors(high, low, accel) result(z)
      real(dp), intent(in) :: high(:) !! Chronological high prices, allowing leading NaNs only.
      real(dp), intent(in) :: low(:) !! Aligned low prices of the same length as high.
      real(dp), intent(in), optional :: accel(:) !! Two acceleration settings, default [0.02,0.2]; upstream equality cap retained.
      real(dp) :: z(size(high)), gaps(size(high)), a(2), ep, prior_ep, af, prior_af, lo, hi
      integer :: b, i, signal, prior_signal
      z = nan_value()
      if (size(low) /= size(high)) return
      b = max(first_valid(high), first_valid(low))
      if (b >= size(high)) return
      a = [0.02_dp, 0.2_dp]
      if (present(accel)) then
         if (size(accel) /= 2) return
         a = accel
      end if
      if (any(ieee_is_nan(a))) return
      gaps = run_sd(high - low, size(high) - b + 1)
      z(b) = low(b) - gaps(size(high))
      ep = high(b)
      af = a(1)
      signal = 1
      do i = b + 1, size(high)
         prior_signal = signal
         prior_ep = ep
         prior_af = af
         lo = min(low(i - 1), low(i))
         hi = max(high(i - 1), high(i))
         if (prior_signal == 1) then
            signal = -1
            if (low(i) > z(i - 1)) signal = 1
            ep = max(hi, prior_ep)
         else
            signal = 1
            if (high(i) < z(i - 1)) signal = -1
            ep = min(lo, prior_ep)
         end if
         if (signal == prior_signal) then
            z(i) = z(i - 1) + (prior_ep - z(i - 1)) * prior_af
            ! Upstream tests equality, rather than clamping a possible overshoot.
            af = a(2)
            if (prior_af /= a(2)) af = a(1) + prior_af
            if (signal == 1) then
               if (ep <= prior_ep) af = prior_af
               z(i) = min(z(i), lo)
            else
               if (ep >= prior_ep) af = prior_af
               z(i) = max(z(i), hi)
            end if
         else
            af = a(1)
            z(i) = ep
         end if
      end do
   end function

   pure function zigzag_vector(price, change, percent, retrace, last_extreme) result(z)
      real(dp), intent(in) :: price(:) !! Close prices, allowing leading NaNs only.
      real(dp), intent(in), optional :: change !! Minimum move, default 10 percentage points.
      logical, intent(in), optional :: percent !! Interpret change in percentage points (default), else price units.
      logical, intent(in), optional :: retrace !! Use a fraction of the preceding move; default false.
      logical, intent(in), optional :: last_extreme !! Choose the last equal extreme (default), else the first.
      real(dp) :: z(size(price)), hl(size(price), 2)
      hl(:, 1) = price
      hl(:, 2) = price
      z = zigzag_matrix(hl, change, percent, retrace, last_extreme)
   end function

   pure function zigzag_matrix(hl, change, percent, retrace, last_extreme) result(z)
      real(dp), intent(in) :: hl(:, :) !! Time rows and one close or two high/low columns; leading NaNs allowed.
      real(dp), intent(in), optional :: change !! Minimum move, default 10 percentage points.
      logical, intent(in), optional :: percent !! Interpret change in percentage points (default), else price units.
      logical, intent(in), optional :: retrace !! Use a fraction of the preceding move; default false.
      logical, intent(in), optional :: last_extreme !! Choose the last equal extreme (default), else the first.
      real(dp) :: z(size(hl, 1)), high(size(hl, 1)), low(size(hl, 1))
      real(dp) :: amount, reference, inflection, extreme_min, extreme_max, local_min, local_max
      integer :: b, i, signal, ref_index, inf_index, left, j
      logical :: use_percent, use_retrace, use_last
      z = nan_value()
      if (size(hl, 2) /= 1 .and. size(hl, 2) /= 2) return
      high = hl(:, 1)
      low = hl(:, size(hl, 2))
      b = max(first_valid(high), first_valid(low))
      if (b >= size(high)) return
      amount = 10.0_dp
      if (present(change)) amount = change
      if (ieee_is_nan(amount) .or. amount < 0.0_dp) return
      use_percent = .true.
      if (present(percent)) use_percent = percent
      use_retrace = .false.
      if (present(retrace)) use_retrace = retrace
      use_last = .true.
      if (present(last_extreme)) use_last = last_extreme
      if (use_percent) amount = amount / 100.0_dp
      reference = (high(b) + low(b)) / 2.0_dp
      ref_index = b
      inflection = (high(b + 1) + low(b + 1)) / 2.0_dp
      inf_index = b + 1
      signal = 0
      do i = b + 1, size(high)
         if (use_percent) then
            extreme_min = inflection * (1.0_dp - amount)
            extreme_max = inflection * (1.0_dp + amount)
         else
            extreme_min = inflection - amount
            extreme_max = inflection + amount
         end if
         local_max = max(inflection, high(i))
         local_min = min(inflection, low(i))
         if (signal == 0) then
            if (use_retrace) then
               signal = -1
               if (inflection >= reference) signal = 1
            else
               if (local_min <= extreme_min) signal = -1
               if (local_max >= extreme_max) signal = 1
            end if
         end if
         if (signal == -1) then
            if (low(i) == local_min) then
               if (use_last .or. low(i) /= low(i - 1)) then
                  inflection = low(i)
                  inf_index = i
               end if
            end if
            if (use_retrace) extreme_max = inflection + (reference - inflection) * amount
            if (high(i) >= extreme_max) then
               z(ref_index) = reference
               reference = inflection
               ref_index = inf_index
               inflection = high(i)
               inf_index = i
               signal = 1
               cycle
            end if
         end if
         if (signal == 1) then
            if (high(i) == local_max) then
               if (use_last .or. high(i) /= high(i - 1)) then
                  inflection = high(i)
                  inf_index = i
               end if
            end if
            if (use_retrace) extreme_min = inflection - (inflection - reference) * amount
            if (low(i) <= extreme_min) then
               z(ref_index) = reference
               reference = inflection
               ref_index = inf_index
               inflection = low(i)
               inf_index = i
               signal = -1
               cycle
            end if
         end if
      end do
      z(ref_index) = reference
      z(inf_index) = inflection
      ! zoo::na.approx interpolates interior gaps but leaves end gaps missing.
      left = 0
      do i = b, size(z)
         if (ieee_is_nan(z(i))) cycle
         if (left > 0) then
            do j = left + 1, i - 1
               z(j) = z(left) + (z(i) - z(left)) * real(j - left, dp) / real(i - left, dp)
            end do
         end if
         left = i
      end do
   end function

   pure function adj_ratios(splits, dividends, close) result(z)
      real(dp), intent(in), optional :: splits(:) !! Aligned split ratios, using NaN where no split occurs.
      real(dp), intent(in), optional :: dividends(:) !! Aligned cash dividends, NaN for no event; close required.
      real(dp), intent(in), optional :: close(:) !! Aligned closing prices, required for dividends; no date merging occurs.
      real(dp), allocatable :: z(:, :)
      integer :: m, i
      m = 0
      if (present(splits)) m = max(m, size(splits))
      if (present(dividends)) m = max(m, size(dividends))
      if (present(close)) m = max(m, size(close))
      allocate(z(m, 2))
      z = nan_value()
      if (m == 0) return
      if (present(splits)) then
         if (size(splits) /= m) return
      end if
      if (present(dividends)) then
         if (size(dividends) /= m .or. .not. present(close)) return
      end if
      if (present(close)) then
         if (size(close) /= m) return
      end if
      z(m, :) = 1.0_dp
      do i = m, 2, -1
         z(i - 1, :) = z(i, :)
         if (present(splits)) then
            if (.not. ieee_is_nan(splits(i))) z(i - 1, 1) = z(i, 1) * splits(i)
         end if
         if (present(dividends) .and. present(close)) then
            if (.not. ieee_is_nan(dividends(i))) z(i - 1, 2) = z(i, 2) * &
               (1.0_dp - divide(dividends(i), close(i - 1)))
         end if
      end do
   end function

   pure function aroon_vector(price, n) result(z)
      real(dp), intent(in) :: price(:) !! Close prices with optional leading NaNs.
      integer, intent(in), optional :: n !! Positive lookback, default 20; subsequent extrema cover n+1 values.
      real(dp) :: z(size(price), 3), hl(size(price), 1)
      hl(:, 1) = price
      z = aroon_matrix(hl, n)
   end function

   pure function aroon_matrix(hl, n) result(z)
      real(dp), intent(in) :: hl(:, :) !! Time rows and one close or two high/low columns.
      integer, intent(in), optional :: n !! Positive lookback, default 20; subsequent extrema cover n+1 values.
      real(dp) :: z(size(hl, 1), 3)
      integer :: width
      z = nan_value()
      if (size(hl, 2) /= 1 .and. size(hl, 2) /= 2) return
      width = 20
      if (present(n)) width = n
      z(:, 1) = aroon_max(hl(:, 1), width)
      z(:, 2) = aroon_max(-hl(:, size(hl, 2)), width)
      z(:, 3) = z(:, 1) - z(:, 2)
   end function

   pure function aroon_max(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological values, including upstream leading-NaN initialization behavior.
      integer, intent(in) :: n !! Positive lookback; first result uses n prices, later windows use n+1.
      real(dp) :: z(size(x)), maximum
      integer :: first, age, i, j
      z = nan_value()
      first = first_valid(x)
      if (n < 1 .or. n > size(x) - first + 1) return
      ! Deliberately initialize from x(1), matching the original C routine.
      maximum = x(1)
      age = 0
      do i = 1, first + n - 2
         if (x(i) >= maximum) then
            maximum = x(i)
            age = 0
         end if
         age = age + 1
      end do
      do i = first + n - 1, size(x)
         if (age > n) then
            maximum = x(i)
            age = 0
            do j = 1, n
               if (x(i - j) > maximum) then
                  maximum = x(i - j)
                  age = j
               end if
            end do
         else
            if (x(i) >= maximum) then
               maximum = x(i)
               age = 0
            end if
         end if
         z(i) = 100.0_dp * real(n - age, dp) / real(n, dp)
         age = age + 1
      end do
   end function

   pure function cti(price, n, slope) result(z)
      real(dp), intent(in) :: price(:) !! Chronological close prices; any NaN invalidates only windows containing it.
      integer, intent(in), optional :: n !! Positive Spearman-correlation window, default 20.
      real(dp), intent(in), optional :: slope !! Desired trend-line slope, default 1; zero yields NaNs.
      real(dp) :: z(size(price)), ranks(size(price)), trend(size(price)), direction, center, denom
      integer :: width, i, j, k, less, equal
      z = nan_value()
      width = 20
      if (present(n)) width = n
      direction = 1.0_dp
      if (present(slope)) direction = slope
      if (width < 2 .or. width > size(price)) return
      if (.not. ieee_is_finite(direction) .or. direction == 0.0_dp) return
      center = real(width + 1, dp) / 2.0_dp
      do j = 1, width
         trend(j) = sign(1.0_dp, direction) * (real(j, dp) - center)
      end do
      do i = width, size(price)
         if (any(ieee_is_nan(price(i - width + 1:i)))) cycle
         do j = 1, width
            less = 0
            equal = 0
            do k = i - width + 1, i
               if (price(k) < price(i - width + j)) less = less + 1
               if (price(k) == price(i - width + j)) equal = equal + 1
            end do
            ranks(j) = real(less, dp) + real(equal + 1, dp) / 2.0_dp - center
         end do
         denom = sqrt(sum(ranks(1:width) ** 2) * sum(trend(1:width) ** 2))
         z(i) = divide(sum(ranks(1:width) * trend(1:width)), denom)
      end do
   end function

   pure function volatility_vector(price, n, calc, annual_periods, mean0, alpha, k, discrete) result(z)
      real(dp), intent(in) :: price(:) !! Chronological close prices; only calc='close' accepts this form.
      integer, intent(in), optional :: n !! Number of price periods, default 10; close estimator uses n-1 returns.
      character(*), intent(in), optional :: calc !! Estimator name or unique prefix; default 'close'.
      real(dp), intent(in), optional :: annual_periods !! Annualization periods, default 260, nonnegative.
      logical, intent(in), optional :: mean0 !! Use a zero return mean for close estimator; default false.
      real(dp), intent(in), optional :: alpha !! Yang-Zhang weight parameter, default 1.34, ignored when k given.
      real(dp), intent(in), optional :: k !! Explicit Yang-Zhang weight; documented upstream override repaired.
      logical, intent(in), optional :: discrete !! Simple returns for close estimator; default log returns.
      real(dp) :: z(size(price)), ohlc(size(price), 1)
      ohlc(:, 1) = price
      z = volatility_matrix(ohlc, n, calc, annual_periods, mean0, alpha, k, discrete)
   end function

   pure function volatility_matrix(ohlc, n, calc, annual_periods, mean0, alpha, k, discrete) result(z)
      real(dp), intent(in) :: ohlc(:, :) !! Time rows with open/high/low/close columns, or one close column.
      integer, intent(in), optional :: n !! Number of periods, default 10; close estimator uses n-1 returns.
      character(*), intent(in), optional :: calc !! close, garman.klass, parkinson, rogers.satchell, gk.yz or yang.zhang.
      real(dp), intent(in), optional :: annual_periods !! Annualization periods, default 260, nonnegative.
      logical, intent(in), optional :: mean0 !! Use zero mean in close-to-close variance; default false.
      real(dp), intent(in), optional :: alpha !! Yang-Zhang parameter, default 1.34; ignored when k supplied.
      real(dp), intent(in), optional :: k !! Explicit Yang-Zhang weight; repairs the documented R override.
      logical, intent(in), optional :: discrete !! Simple returns for close estimator; default logarithmic.
      real(dp) :: z(size(ohlc, 1)), r(size(ohlc, 1)), co(size(ohlc, 1)), gap(size(ohlc, 1))
      real(dp) :: rs(size(ohlc, 1)), hl(size(ohlc, 1)), scale, shape, weight
      integer :: width, method, j, matches
      logical :: zero_mean
      character(16), parameter :: names(6) = [character(16) :: &
         'close', 'garman.klass', 'parkinson', 'rogers.satchell', 'gk.yz', 'yang.zhang']
      character(:), allocatable :: choice
      z = nan_value()
      width = 10
      if (present(n)) width = n
      scale = 260.0_dp
      if (present(annual_periods)) scale = annual_periods
      if (width < 1 .or. width > size(ohlc, 1) .or. scale < 0.0_dp) return
      choice = 'close'
      if (present(calc)) choice = trim(calc)
      method = 0
      matches = 0
      if (len(choice) == 0) return
      do j = 1, size(names)
         if (index(trim(names(j)), choice) == 1) then
            method = j
            matches = matches + 1
         end if
      end do
      if (matches /= 1) return
      if (method == 1) then
         if (size(ohlc, 2) == 1) then
            r = roc(ohlc(:, 1), 1, discrete)
         else if (size(ohlc, 2) >= 4) then
            r = roc(ohlc(:, 4), 1, discrete)
         else
            return
         end if
         zero_mean = .false.
         if (present(mean0)) zero_mean = mean0
         if (zero_mean) then
            z = sqrt(scale) * sqrt(divide(run_sum(r ** 2, width - 1), real(width - 2, dp)))
         else
            z = sqrt(scale) * run_sd(r, width - 1)
         end if
         return
      end if
      if (size(ohlc, 2) < 4) return
      hl = log(divide(ohlc(:, 2), ohlc(:, 3)))
      co = log(divide(ohlc(:, 4), ohlc(:, 1)))
      select case(method)
      case(2)
         r = 0.5_dp * hl ** 2 - (2.0_dp * log(2.0_dp) - 1.0_dp) * co ** 2
         z = sqrt(scale / real(width, dp) * run_sum(r, width))
      case(3)
         z = sqrt(scale / (4.0_dp * real(width, dp) * log(2.0_dp)) * run_sum(hl ** 2, width))
      case(4)
         r = log(divide(ohlc(:, 2), ohlc(:, 4))) * log(divide(ohlc(:, 2), ohlc(:, 1))) + &
            log(divide(ohlc(:, 3), ohlc(:, 4))) * log(divide(ohlc(:, 3), ohlc(:, 1)))
         z = sqrt(scale / real(width, dp) * run_sum(r, width))
      case(5)
         gap = log(divide(ohlc(:, 1), shift(ohlc(:, 4), 1)))
         r = gap ** 2 + 0.5_dp * hl ** 2 - (2.0_dp * log(2.0_dp) - 1.0_dp) * co ** 2
         z = sqrt(scale / real(width, dp) * run_sum(r, width))
      case(6)
         if (width < 2) return
         shape = 1.34_dp
         if (present(alpha)) shape = alpha
         weight = divide(shape - 1.0_dp, shape + real(width + 1, dp) / real(width - 1, dp))
         if (present(k)) weight = k
         gap = log(divide(ohlc(:, 1), shift(ohlc(:, 4), 1)))
         r = log(divide(ohlc(:, 2), ohlc(:, 4))) * log(divide(ohlc(:, 2), ohlc(:, 1))) + &
            log(divide(ohlc(:, 3), ohlc(:, 4))) * log(divide(ohlc(:, 3), ohlc(:, 1)))
         rs = sqrt(scale / real(width, dp) * run_sum(r, width))
         z = sqrt(scale * run_var(gap, width) + weight * scale * run_var(co, width) + &
            (1.0_dp - weight) * rs ** 2)
      end select
   end function

   pure function roll_sfm(asset, benchmark, n) result(z)
      real(dp), intent(in) :: asset(:) !! Chronological excess security returns, allowing leading NaNs only.
      real(dp), intent(in) :: benchmark(:) !! Same-length excess benchmark returns aligned with asset.
      integer, intent(in), optional :: n !! Rolling regression window, default 60; R-squared needs more than two values.
      real(dp) :: z(size(asset), 3), beta(size(asset)), residual(size(asset)), asset_var(size(asset))
      real(dp) :: a_sum(size(asset)), b_sum(size(asset)), width_real
      integer :: width
      z = nan_value()
      if (size(asset) /= size(benchmark)) return
      width = 60
      if (present(n)) width = n
      if (width < 1 .or. width > size(asset)) return
      width_real = real(width, dp)
      beta = divide(run_cov(asset, benchmark, width), run_var(benchmark, width))
      z(:, 1) = run_mean(asset, width) - beta * run_mean(benchmark, width)
      z(:, 2) = beta
      a_sum = run_sum(asset, width)
      b_sum = run_sum(benchmark, width)
      residual = divide(width_real * run_sum(asset ** 2, width) - a_sum ** 2 - &
         beta ** 2 * (width_real * run_sum(benchmark ** 2, width) - b_sum ** 2), &
         width_real * (width_real - 2.0_dp))
      asset_var = run_var(asset, width) * divide(width_real - 1.0_dp, width_real - 2.0_dp)
      z(:, 3) = 1.0_dp - divide(residual, asset_var)
   end function
end module ttr_advanced
