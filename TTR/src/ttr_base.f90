! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from TTR 0.24.4, Joshua M. Ulrich and contributors; see NOTICE.md.
module ttr_base
   use r_kinds, only: dp
   use r_descriptive, only: r_mean, r_variance, r_covariance, r_correlation, r_sd
   use r_quantiles, only: r_median
   use r_ordering, only: r_sort_values_in_place
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf
   implicit none
   private
   public :: dp, na_check, run_sum, run_mean, run_min, run_max, run_median, run_mad
   public :: run_cov, run_cor, run_var, run_sd, run_percent_rank, wilder_sum
   public :: sma, ema, dema, wma, vwma, vwap, evwma, zlema, hma, alma
   public :: momentum, roc, lags, growth
   ! Shared implementation utilities, intentionally public to the indicator module.
   public :: nan_value, first_valid, moving_average, shift, divide, prices, high_low_close
contains
   pure function nan_value() result(z)
      real(dp) :: z
      z = ieee_value(0.0_dp, ieee_quiet_nan)
   end function

   pure elemental function divide(a, b) result(z)
      real(dp), intent(in) :: a !! Numerator; missing values propagate.
      real(dp), intent(in) :: b !! Denominator; zero yields signed infinity or NaN.
      real(dp) :: z
      z = nan_value()
      if (ieee_is_nan(a) .or. ieee_is_nan(b)) return
      if (b == 0.0_dp) then
         if (a /= 0.0_dp) z = sign(ieee_value(0.0_dp, ieee_positive_inf), a) * sign(1.0_dp, b)
      else
         z = a / b
      end if
   end function

   pure subroutine na_check(x, n, n_missing, begin, valid)
      real(dp), intent(in) :: x(:) !! Chronological observations with optional leading NaNs.
      integer, intent(in) :: n !! Nonnegative initialization window offset.
      integer, intent(out) :: n_missing !! Count of missing observations, matching R NAs.
      integer, intent(out) :: begin !! One-based R beg index, n plus missing count.
      logical, intent(out) :: valid !! True when all missing observations precede observed data.
      n_missing = count(ieee_is_nan(x))
      begin = n + n_missing
      valid = n >= 0
      if (n_missing < size(x)) valid = valid .and. .not. any(ieee_is_nan(x(n_missing + 1:)))
   end subroutine

   pure integer function first_valid(x) result(k)
      real(dp), intent(in) :: x(:) !! Series whose missing observations must be leading only.
      integer :: m, b
      logical :: ok
      call na_check(x, 0, m, b, ok)
      k = size(x) + 1
      if (ok) k = m + 1
   end function

   pure function prices(x) result(z)
      real(dp), intent(in) :: x(:, :) !! One close column or three high, low, close columns.
      real(dp) :: z(size(x, 1))
      z = nan_value()
      if (size(x, 2) == 1) z = x(:, 1)
      if (size(x, 2) == 3) z = sum(x, dim = 2) / 3.0_dp
   end function

   pure subroutine high_low_close(x, h, l, c)
      real(dp), intent(in) :: x(:, :) !! One close column or three high, low, close columns.
      real(dp), intent(out) :: h(size(x, 1)) !! High prices, or close for univariate input.
      real(dp), intent(out) :: l(size(x, 1)) !! Low prices, or close for univariate input.
      real(dp), intent(out) :: c(size(x, 1)) !! Closing prices.
      h = nan_value()
      l = nan_value()
      c = nan_value()
      if (size(x, 2) == 1) then
         h = x(:, 1)
         l = h
         c = h
      else if (size(x, 2) == 3) then
         h = x(:, 1)
         l = x(:, 2)
         c = x(:, 3)
      end if
   end subroutine

   pure function shift(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological vector to lag or lead.
      integer, intent(in) :: n !! Positive for lag, negative for lead; vacated positions are NaN.
      real(dp) :: z(size(x))
      integer :: i, j
      z = nan_value()
      do i = 1, size(x)
         j = i - n
         if (j >= 1 .and. j <= size(x)) z(i) = x(j)
      end do
   end function

   pure function rolling(x, n, operation, cumulative, non_unique) result(z)
      real(dp), intent(in) :: x(:) !! Univariate observations, allowing leading NaNs only.
      integer, intent(in) :: n !! Minimum count; fixed window width unless cumulative.
      integer, intent(in) :: operation !! 1 sum, 2 mean, 3 minimum, 4 maximum, 5 median.
      logical, intent(in), optional :: cumulative !! Expanding from first observation; default false.
      integer, intent(in), optional :: non_unique !! Even median: -1 lower, 0 mean, 1 upper.
      real(dp) :: z(size(x)), total
      integer :: i, j, k, rule
      logical :: cumul
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      cumul = .false.
      if (present(cumulative)) cumul = cumulative
      rule = 0
      if (present(non_unique)) rule = non_unique
      if (operation == 1 .or. operation == 2) then
         ! Preserve C runsum addition/subtraction order: rank-based indicators
         ! can distinguish even tiny negative rounding differences from ties.
         total = 0.0_dp
         do i = k, k + n - 1
            total = total + x(i)
         end do
         z(k + n - 1) = total
         do i = k + n, size(x)
            total = total + x(i)
            if (.not. cumul) total = total - x(i - n)
            z(i) = total
         end do
         if (operation == 2) then
            do i = k + n - 1, size(x)
               j = n
               if (cumul) j = i - k + 1
               z(i) = z(i) / real(j, dp)
            end do
         end if
         return
      end if
      do i = k + n - 1, size(x)
         j = i - n + 1
         if (cumul) j = k
         select case(operation)
         case(1)
            z(i) = sum(x(j:i))
         case(2)
            z(i) = r_mean(x(j:i))
         case(3)
            z(i) = minval(x(j:i))
         case(4)
            z(i) = maxval(x(j:i))
         case(5)
            z(i) = median_rule(x(j:i), rule)
         end select
      end do
   end function

   pure function median_rule(x, rule) result(z)
      real(dp), intent(in) :: x(:) !! Nonmissing observations to summarize.
      integer, intent(in) :: rule !! Even sample convention: negative lower, zero mean, positive upper.
      real(dp) :: z, work(size(x))
      integer :: m
      z = r_median(x)
      m = size(x)
      if (m == 0 .or. mod(m, 2) /= 0 .or. rule == 0) return
      work = x
      call r_sort_values_in_place(work)
      z = work(m / 2)
      if (rule > 0) z = work(m / 2 + 1)
   end function

   pure function run_sum(x, n, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations, allowing leading NaNs only.
      integer, intent(in) :: n !! Positive window width or minimum expanding sample count.
      logical, intent(in), optional :: cumulative !! Use all observed history; default false.
      real(dp) :: z(size(x))
      z = rolling(x, n, 1, cumulative)
   end function

   pure function run_mean(x, n, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations, allowing leading NaNs only.
      integer, intent(in) :: n !! Positive window width or minimum expanding sample count.
      logical, intent(in), optional :: cumulative !! Use all observed history; default false.
      real(dp) :: z(size(x))
      z = rolling(x, n, 2, cumulative)
   end function

   pure function run_min(x, n, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations, allowing leading NaNs only.
      integer, intent(in) :: n !! Positive window width or minimum expanding sample count.
      logical, intent(in), optional :: cumulative !! Use all observed history; default false.
      real(dp) :: z(size(x))
      z = rolling(x, n, 3, cumulative)
   end function

   pure function run_max(x, n, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations, allowing leading NaNs only.
      integer, intent(in) :: n !! Positive window width or minimum expanding sample count.
      logical, intent(in), optional :: cumulative !! Use all observed history; default false.
      real(dp) :: z(size(x))
      z = rolling(x, n, 4, cumulative)
   end function

   pure function run_median(x, n, cumulative, non_unique) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations, allowing leading NaNs only.
      integer, intent(in) :: n !! Positive window width or minimum expanding sample count.
      logical, intent(in), optional :: cumulative !! Use all observed history; default false.
      integer, intent(in), optional :: non_unique !! Even median rule: -1 lower, 0 mean (default), 1 upper.
      real(dp) :: z(size(x))
      z = rolling(x, n, 5, cumulative, non_unique)
   end function

   pure function run_cov(x, y, n, sample, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! First chronological series.
      real(dp), intent(in) :: y(:) !! Second series aligned with x; equal length required.
      integer, intent(in) :: n !! Positive window width or minimum expanding sample count.
      logical, intent(in), optional :: sample !! Divide by count minus one (default), otherwise count.
      logical, intent(in), optional :: cumulative !! Expanding covariance; default false.
      real(dp) :: z(size(x))
      integer :: i, j, k, ddof
      logical :: cumul
      z = nan_value()
      if (size(x) /= size(y)) return
      k = max(first_valid(x), first_valid(y))
      if (n < 1 .or. n > size(x) - k + 1) return
      ddof = 1
      if (present(sample)) then
         if (.not. sample) ddof = 0
      end if
      cumul = .false.
      if (present(cumulative)) cumul = cumulative
      do i = k + n - 1, size(x)
         j = i - n + 1
         if (cumul) j = k
         z(i) = r_covariance(x(j:i), y(j:i), ddof = ddof)
      end do
   end function

   pure function run_var(x, n, sample, cumulative, y) result(z)
      real(dp), intent(in) :: x(:) !! Chronological series.
      integer, intent(in) :: n !! Positive window width or minimum expanding count.
      logical, intent(in), optional :: sample !! Sample variance by default; false for population.
      logical, intent(in), optional :: cumulative !! Expanding variance; default false.
      real(dp), intent(in), optional :: y(:) !! Optional second series, matching R runVar covariance mode.
      real(dp) :: z(size(x))
      if (present(y)) then
         z = run_cov(x, y, n, sample, cumulative)
      else
         z = run_cov(x, x, n, sample, cumulative)
      end if
   end function

   pure function run_sd(x, n, sample, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! Chronological series with leading NaNs permitted.
      integer, intent(in) :: n !! Positive window width or minimum expanding count.
      logical, intent(in), optional :: sample !! Sample standard deviation by default.
      logical, intent(in), optional :: cumulative !! Expanding statistic; default false.
      real(dp) :: z(size(x))
      z = sqrt(run_var(x, n, sample, cumulative))
   end function

   pure function run_cor(x, y, n, sample, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! First chronological series.
      real(dp), intent(in) :: y(:) !! Aligned second series of equal length.
      integer, intent(in) :: n !! Positive window width or minimum expanding count.
      logical, intent(in), optional :: sample !! Sample covariance normalization; default true.
      logical, intent(in), optional :: cumulative !! Expanding correlation; default false.
      real(dp) :: z(size(x))
      z = nan_value()
      if (size(x) /= size(y)) return
      z = divide(run_cov(x, y, n, sample, cumulative), &
         run_sd(x, n, sample, cumulative) * run_sd(y, n, sample, cumulative))
   end function

   pure function run_mad(x, n, center, median, constant, non_unique, cumulative) result(z)
      real(dp), intent(in) :: x(:) !! Chronological series with leading NaNs permitted.
      integer, intent(in) :: n !! Positive window width or minimum expanding count.
      real(dp), intent(in), optional :: center(:) !! Per-row centers; default rolling median.
      logical, intent(in), optional :: median !! Median absolute deviation (default) or mean deviation.
      real(dp), intent(in), optional :: constant !! Scale for median deviation only; default 1.4826.
      integer, intent(in), optional :: non_unique !! Even deviation median: -1 lower, 0 mean, 1 upper.
      logical, intent(in), optional :: cumulative !! Expanding statistic; default false.
      real(dp) :: z(size(x)), centers(size(x)), scale
      integer :: i, j, k, rule
      logical :: med, cumul
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      centers = run_median(x, n, cumulative)
      if (present(center)) then
         if (size(center) /= size(x)) return
         centers = center
      end if
      med = .true.
      if (present(median)) med = median
      cumul = .false.
      if (present(cumulative)) cumul = cumulative
      rule = 0
      if (present(non_unique)) rule = non_unique
      scale = 1.4826_dp
      if (present(constant)) scale = constant
      do i = k + n - 1, size(x)
         j = i - n + 1
         if (cumul) j = k
         if (med) then
            z(i) = scale * median_rule(abs(x(j:i) - centers(i)), rule)
         else
            z(i) = r_mean(abs(x(j:i) - centers(i)))
         end if
      end do
   end function

   pure function run_percent_rank(x, n, cumulative, exact_multiplier) result(z)
      real(dp), intent(in) :: x(:) !! Chronological series; leading NaNs are allowed.
      integer, intent(in) :: n !! Window width or cumulative warm-up count.
      logical, intent(in), optional :: cumulative !! Preserve TTR cumulative initialization; default false.
      real(dp), intent(in), optional :: exact_multiplier !! Tie weight in [0,1], default 0.5.
      real(dp) :: z(size(x)), mult, num, diff
      integer :: i, j, k, left
      logical :: cumul
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      cumul = .false.
      if (present(cumulative)) cumul = cumulative
      mult = 0.5_dp
      if (present(exact_multiplier)) mult = exact_multiplier
      if (mult < 0.0_dp .or. mult > 1.0_dp) return
      if (.not. cumul .and. n == 1) then
         z = mult
         return
      end if
      if (cumul .and. n == 1) k = 1
      do i = k + n - 1, size(x)
         left = i - n + 1
         if (cumul) left = 1
         num = mult
         do j = left, i - 1
            diff = x(j) - x(i)
            if (diff < 0.0_dp) then
               num = num + 1.0_dp
            else if (abs(diff) < 1.0e-8_dp) then
               num = num + mult
            end if
         end do
         z(i) = num / real(i - left + 1, dp)
      end do
      if (cumul) z(k + n - 1) = mult
   end function

   pure function wilder_sum(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations; leading NaNs are allowed.
      integer, intent(in) :: n !! Positive smoothing period; seed follows TTR wilderSum.
      real(dp) :: z(size(x)), decay
      integer :: k, b, i
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      if (n == 1 .and. k > 1) return
      b = k + n - 1
      decay = real(n - 1, dp) / real(n, dp)
      z(b) = x(b) + sum(x(k:b - 1)) * decay
      do i = b + 1, size(x)
         z(i) = x(i) + z(i - 1) * decay
      end do
   end function

   pure function sma(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations; leading NaNs are allowed.
      integer, intent(in) :: n !! Positive moving-average window.
      real(dp) :: z(size(x))
      z = run_mean(x, n)
   end function

   pure function ema(x, n, wilder, ratio) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations; leading NaNs are allowed.
      integer, intent(in) :: n !! Positive seed window; required even when ratio is provided.
      logical, intent(in), optional :: wilder !! Use decay 1/n rather than 2/(n+1); default false.
      real(dp), intent(in), optional :: ratio !! Explicit positive decay, overriding the default.
      real(dp) :: z(size(x)), a
      integer :: k, b, i
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      a = 2.0_dp / real(n + 1, dp)
      if (present(wilder)) then
         if (wilder) a = 1.0_dp / real(n, dp)
      end if
      if (present(ratio)) a = ratio
      if (a <= 0.0_dp) return
      b = k + n - 1
      z(b) = r_mean(x(k:b))
      do i = b + 1, size(x)
         z(i) = a * x(i) + (1.0_dp - a) * z(i - 1)
      end do
   end function

   pure function dema(x, n, v, wilder, ratio) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations.
      integer, intent(in) :: n !! Positive EMA seed window.
      real(dp), intent(in), optional :: v !! Generalized DEMA factor in [0,1], default 1.
      logical, intent(in), optional :: wilder !! Wilder EMA decay; default false.
      real(dp), intent(in), optional :: ratio !! Explicit positive EMA decay.
      real(dp) :: z(size(x)), a(size(x)), weight
      weight = 1.0_dp
      if (present(v)) weight = v
      z = nan_value()
      if (weight < 0.0_dp .or. weight > 1.0_dp) return
      a = ema(x, n, wilder, ratio)
      z = (1.0_dp + weight) * a - weight * ema(a, n, wilder, ratio)
   end function

   pure function wma(x, n, weights) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations.
      integer, intent(in) :: n !! Positive moving-window length.
      real(dp), intent(in), optional :: weights(:) !! Fixed n weights or series-length weights; default 1:n.
      real(dp) :: z(size(x)), w(max(0, n))
      integer :: k, i
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      w = [(real(i, dp), i = 1, n)]
      if (present(weights)) then
         if (size(weights) == n) then
            if (any(ieee_is_nan(weights))) return
            w = weights
         else if (size(weights) == size(x)) then
            z = divide(run_sum(x * weights, n), run_sum(weights, n))
            return
         else
            return
         end if
      end if
      do i = k + n - 1, size(x)
         z(i) = divide(sum(x(i - n + 1:i) * w), sum(w))
      end do
   end function

   pure function vwma(x, volume, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      real(dp), intent(in) :: volume(:) !! Aligned volume weights, matching length of prices.
      integer, intent(in) :: n !! Positive moving-window length.
      real(dp) :: z(size(x))
      z = wma(x, n, volume)
   end function

   pure function vwap(x, volume, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      real(dp), intent(in) :: volume(:) !! Aligned volume weights, matching length of prices.
      integer, intent(in) :: n !! Positive moving-window length.
      real(dp) :: z(size(x))
      z = vwma(x, volume, n)
   end function

   pure function evwma(x, volume, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      real(dp), intent(in) :: volume(:) !! One constant volume or series-length aligned volumes.
      integer, intent(in) :: n !! Positive rolling volume-sum window.
      real(dp) :: z(size(x)), v(size(x)), vs(size(x))
      integer :: k, b, i
      z = nan_value()
      if (size(volume) == 1) then
         v = volume(1)
      else if (size(volume) == size(x)) then
         v = volume
      else
         return
      end if
      k = max(first_valid(x), first_valid(v))
      if (n < 1 .or. n > size(x) - k + 1) return
      vs = run_sum(v, n)
      b = k + n - 1
      z(b) = x(b)
      do i = b + 1, size(x)
         z(i) = divide((vs(i) - v(i)) * z(i - 1) + v(i) * x(i), vs(i))
      end do
   end function

   pure function zlema(x, n, ratio) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices with optional leading NaNs.
      integer, intent(in) :: n !! Positive seed window, required with ratio too.
      real(dp), intent(in), optional :: ratio !! Positive decay; lag is 1/ratio as in upstream C.
      real(dp) :: z(size(x)), a, lag, weight, value
      integer :: k, b, i, loc
      z = nan_value()
      k = first_valid(x)
      if (n < 1 .or. n > size(x) - k + 1) return
      a = 2.0_dp / real(n + 1, dp)
      if (present(ratio)) a = ratio
      if (a <= 0.0_dp) return
      lag = 1.0_dp / a
      weight = modulo(lag, 1.0_dp)
      b = k + n - 1
      z(b) = r_mean(x(k:b))
      do i = b + 1, size(x)
         loc = int(real(i - 1, dp) - lag) + 1
         if (loc < k .or. loc + 1 > size(x)) then
            z = nan_value()
            return
         end if
         value = 2.0_dp * x(i) - ((1.0_dp - weight) * x(loc) + weight * x(loc + 1))
         z(i) = a * value + (1.0_dp - a) * z(i - 1)
      end do
   end function

   pure function hma(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations.
      integer, intent(in) :: n !! Hull window, at least two.
      real(dp) :: z(size(x))
      z = nan_value()
      if (n < 2) return
      z = wma(2.0_dp * wma(x, n / 2) - wma(x, n), int(sqrt(real(n, dp))))
   end function

   pure function alma(x, n, offset, sigma) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations.
      integer, intent(in) :: n !! Positive Gaussian weight window.
      real(dp), intent(in), optional :: offset !! Gaussian center fraction in [0,1], default 0.85.
      real(dp), intent(in), optional :: sigma !! Positive inverse-width parameter, default 6.
      real(dp) :: z(size(x)), w(max(0, n)), off, sig, width, center
      integer :: i
      z = nan_value()
      if (n < 1) return
      off = 0.85_dp
      sig = 6.0_dp
      if (present(offset)) off = offset
      if (present(sigma)) sig = sigma
      if (off < 0.0_dp .or. off > 1.0_dp .or. sig <= 0.0_dp) return
      center = real(floor(off * real(n - 1, dp)), dp)
      width = real(n, dp) / sig
      w = [(exp(-(real(i - 1, dp) - center) ** 2 / (2.0_dp * width ** 2)), i = 1, n)]
      z = wma(x, n, w)
   end function

   pure function moving_average(x, n, method, default) result(z)
      real(dp), intent(in) :: x(:) !! Chronological series to smooth.
      integer, intent(in) :: n !! Positive smoothing period.
      character(*), intent(in), optional :: method !! SMA, EMA, Wilder, WMA, DEMA, HMA, ALMA or ZLEMA.
      character(*), intent(in), optional :: default !! Method when omitted; otherwise SMA.
      real(dp) :: z(size(x))
      character(:), allocatable :: name
      name = 'SMA'
      if (present(default)) name = default
      if (present(method)) name = method
      select case(name)
      case('SMA', 'sma')
         z = sma(x, n)
      case('EMA', 'ema')
         z = ema(x, n)
      case('Wilder', 'wilder')
         z = ema(x, n, .true.)
      case('WMA', 'wma')
         z = wma(x, n)
      case('DEMA', 'dema')
         z = dema(x, n)
      case('HMA', 'hma')
         z = hma(x, n)
      case('ALMA', 'alma')
         z = alma(x, n)
      case('ZLEMA', 'zlema')
         z = zlema(x, n)
      case default
         z = nan_value()
      end select
   end function

   pure function momentum(x, n) result(z)
      real(dp), intent(in) :: x(:) !! Chronological observations.
      integer, intent(in) :: n !! Nonnegative lag; output is always padded to input length.
      real(dp) :: z(size(x))
      z = nan_value()
      if (n < 0) return
      z = x - shift(x, n)
   end function

   pure function roc(x, n, discrete) result(z)
      real(dp), intent(in) :: x(:) !! Chronological positive prices for continuous returns.
      integer, intent(in) :: n !! Nonnegative lag; output is always input length.
      logical, intent(in), optional :: discrete !! Simple return if true; log return by default.
      real(dp) :: z(size(x)), prior(size(x))
      logical :: simple
      z = nan_value()
      if (n < 0) return
      simple = .false.
      if (present(discrete)) simple = discrete
      prior = shift(x, n)
      if (simple) then
         z = divide(x, prior) - 1.0_dp
      else
         where (x > 0.0_dp .and. prior > 0.0_dp) z = log(x) - log(prior)
      end if
   end function

   pure function lags(x, n) result(z)
      real(dp), intent(in) :: x(:, :) !! Rows are times; columns are input variables.
      integer, intent(in) :: n !! Nonnegative number of lags; incomplete rows are omitted.
      real(dp) :: z(max(0, size(x, 1) - max(0, n)), size(x, 2) * (max(0, n) + 1))
      integer :: j, m, rows
      z = nan_value()
      if (n < 0 .or. n >= size(x, 1)) return
      m = size(x, 2)
      rows = size(x, 1) - n
      do j = 0, n
         z(:, j * m + 1:(j + 1) * m) = x(n + 1 - j:n + rows - j, :)
      end do
   end function

   pure function growth(x, signals, discrete, n, na_pad) result(z)
      real(dp), intent(in) :: x(:) !! Chronological prices.
      real(dp), intent(in), optional :: signals(:) !! Position multipliers aligned with returned ROC values.
      logical, intent(in), optional :: discrete !! Simple returns when true; log returns by default as R.
      integer, intent(in), optional :: n !! Return lag, default one.
      logical, intent(in), optional :: na_pad !! Keep leading NaNs (default); false returns shortened result.
      real(dp), allocatable :: z(:)
      real(dp) :: r(size(x)), total
      integer :: lag, start, i, j
      lag = 1
      if (present(n)) lag = n
      start = 1
      if (present(na_pad)) then
         if (.not. na_pad) start = max(0, lag) + 1
      end if
      allocate(z(max(0, size(x) - start + 1)))
      z = nan_value()
      if (present(signals)) then
         if (size(signals) /= size(z)) return
      end if
      r = roc(x, lag, discrete)
      total = 1.0_dp
      do i = start, size(x)
         j = i - start + 1
         if (present(signals)) then
            total = total * (1.0_dp + r(i) * signals(j))
         else
            total = total * (1.0_dp + r(i))
         end if
         z(j) = total
      end do
   end function
end module ttr_base
