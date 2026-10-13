! SPDX-License-Identifier: LGPL-2.0-or-later
! Derived from fTrading 3042.79 by Diethelm Wuertz and contributors.
module ftrading_base
   use r_kinds, only: dp
   use r_descriptive, only: r_mean, r_variance
   use r_quantiles, only: r_median
   use ttr_base, only: nan_value, divide, ttr_mean => run_mean, ttr_min => run_min, ttr_max => run_max
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   private
   public :: dp, max_draw_down, sharpe_ratio, sterling_ratio, roll_fun, roll_var, sma, ewma
   public :: ema_ta, trade_signals, trade_lengths, hit_rate, rolling_operation
   abstract interface
      pure function rolling_operation(window) result(value)
         import dp
         real(dp), intent(in) :: window(:) !! Complete chronological window to reduce to one scalar.
         real(dp) :: value
      end function
   end interface
   interface roll_fun
      module procedure roll_named, roll_callback
   end interface
contains
   pure function scaled_mean(x) result(y)
      real(dp), intent(in) :: x(:) !! Observations rescaled only when direct summation can overflow.
      real(dp) :: y, magnitude
      y = r_mean(x)
      if (size(x) == 0) return
      if (.not. all(ieee_is_finite(x))) return
      if (ieee_is_finite(y)) return
      magnitude = maxval(abs(x))
      if (magnitude > 0.0_dp) y = r_mean(x / magnitude) * magnitude
   end function

   pure function scaled_variance(x) result(y)
      real(dp), intent(in) :: x(:) !! Window whose sample variance is rescaled only after finite-input overflow.
      real(dp) :: y, magnitude, normalized, deviation
      y = r_variance(x)
      if (size(x) < 2) return
      if (.not. all(ieee_is_finite(x))) return
      if (.not. ieee_is_nan(y)) return
      magnitude = maxval(abs(x))
      if (magnitude == 0.0_dp) return
      normalized = r_variance(x / magnitude)
      if (normalized == 0.0_dp) then
         y = 0.0_dp
      else if (normalized > 0.0_dp) then
         deviation = sqrt(normalized) * magnitude
         y = deviation * deviation
      end if
   end function

   pure subroutine max_draw_down(x, drawdown, from, to)
      real(dp), intent(in) :: x(:) !! Finite univariate levels, usually cumulative log prices.
      real(dp), intent(out) :: drawdown !! Largest absolute peak-to-trough decline; NaN for invalid data.
      integer, allocatable, intent(out) :: from(:) !! Last peak index preceding each tied maximum decline.
      integer, allocatable, intent(out) :: to(:) !! All one-based trough indices attaining the maximum decline.
      real(dp) :: declines(size(x)), peak
      integer :: i, j, last_peak
      drawdown = nan_value()
      allocate(from(0), to(0))
      if (size(x) == 0) return
      if (.not. all(ieee_is_finite(x))) return
      peak = x(1)
      do i = 1, size(x)
         peak = max(peak, x(i))
         declines(i) = peak - x(i)
      end do
      drawdown = maxval(declines)
      deallocate(from, to)
      allocate(from(count(declines == drawdown)), to(count(declines == drawdown)))
      j = 0
      last_peak = 1
      do i = 1, size(x)
         if (declines(i) == 0.0_dp) last_peak = i
         if (declines(i) == drawdown) then
            j = j + 1
            from(j) = last_peak
            to(j) = i
         end if
      end do
   end subroutine

   pure function sharpe_ratio(x, risk_free, scale, rates) result(value)
      real(dp), intent(in) :: x(:) !! Finite chronological levels; differences are the returns.
      real(dp), intent(in), optional :: risk_free !! Scalar return per observation; default zero.
      real(dp), intent(in), optional :: scale !! Annualizing multiplier; default sqrt(250).
      real(dp), intent(in), optional :: rates(:) !! Time-varying risk-free returns, exactly size(x)-1 values.
      real(dp) :: value, factor, rate, y(max(0, size(x) - 1))
      value = nan_value()
      if (size(x) < 3) return
      if (.not. all(ieee_is_finite(x))) return
      factor = sqrt(250.0_dp)
      if (present(scale)) factor = scale
      rate = 0.0_dp
      if (present(risk_free)) rate = risk_free
      y = x(2:) - x(:size(x) - 1)
      if (present(rates)) then
         if (size(rates) /= size(y)) return
         value = factor * divide(scaled_mean(y - rates), sqrt(scaled_variance(y)))
      else
         value = factor * divide(scaled_mean(y - rate), sqrt(scaled_variance(y)))
      end if
   end function

   pure function sterling_ratio(x) result(value)
      real(dp), intent(in) :: x(:) !! Finite chronological levels; at least two observations required.
      real(dp) :: value, drawdown
      integer, allocatable :: from(:), to(:)
      value = nan_value()
      if (size(x) < 2) return
      call max_draw_down(x, drawdown, from, to)
      value = divide(x(size(x)) - x(1), drawdown)
   end function

   pure function roll_callback(x, n, fun, trim, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Chronological observations; NaNs can be removed globally before rolling.
      integer, intent(in) :: n !! Positive number of observations per complete window.
      procedure(rolling_operation), pointer, intent(in) :: fun !! Pure scalar reduction applied independently to each full window.
      logical, intent(in), optional :: trim !! Remove first n-1 unavailable entries; default true.
      logical, intent(in), optional :: na_rm !! Remove all NaNs before constructing windows; default false.
      real(dp), allocatable :: y(:), work(:)
      logical :: remove, trimmed
      integer :: i, offset, m
      if (.not. associated(fun)) then
         allocate(y(0))
         return
      end if
      remove = .false.
      if (present(na_rm)) remove = na_rm
      trimmed = .true.
      if (present(trim)) trimmed = trim
      work = x
      if (remove) work = pack(x, .not. ieee_is_nan(x))
      m = size(work)
      offset = 0
      if (.not. trimmed) offset = max(0, n - 1)
      if (n < 1 .or. n > m) then
         allocate(y(0))
         return
      end if
      allocate(y(m - n + 1 + offset))
      y = nan_value()
      do i = 1, m - n + 1
         y(i + offset) = fun(work(i:i + n - 1))
      end do
   end function

   pure function roll_named(x, n, operation, trim, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Chronological observations, with optional NaNs.
      integer, intent(in) :: n !! Positive complete-window width, no larger than the retained series.
      character(len=*), intent(in) :: operation !! Scalar reduction: mean, sum, min, max, var, sd, or median.
      logical, intent(in), optional :: trim !! Trim unavailable leading entries; default true.
      logical, intent(in), optional :: na_rm !! Remove NaNs from the entire series before rolling; default false.
      real(dp), allocatable :: y(:), work(:), full(:)
      integer :: i, offset, m
      logical :: remove, trimmed
      remove = .false.
      if (present(na_rm)) remove = na_rm
      trimmed = .true.
      if (present(trim)) trimmed = trim
      work = x
      if (remove) work = pack(x, .not. ieee_is_nan(x))
      m = size(work)
      if (n < 1 .or. n > m) then
         allocate(y(0))
         return
      end if
      offset = 0
      if (.not. trimmed) offset = n - 1
      allocate(y(m - n + 1 + offset))
      y = nan_value()
      ! Reuse compatible TTR kernels on each valid window. This preserves fTrading's
      ! local missing-data behavior rather than TTR's leading-NA-only restriction.
      do i = 1, m - n + 1
         if (any(ieee_is_nan(work(i:i + n - 1)))) cycle
         select case(operation)
         case('mean')
            full = ttr_mean(work(i:i + n - 1), n)
            y(i + offset) = full(n)
            if (all(ieee_is_finite(work(i:i + n - 1)))) then
               if (.not. ieee_is_finite(full(n))) y(i + offset) = scaled_mean(work(i:i + n - 1))
            end if
         case('min')
            full = ttr_min(work(i:i + n - 1), n)
            y(i + offset) = full(n)
         case('max')
            full = ttr_max(work(i:i + n - 1), n)
            y(i + offset) = full(n)
         case('median')
            y(i + offset) = r_median(work(i:i + n - 1))
         case('sum')
            y(i + offset) = sum(work(i:i + n - 1))
         case('var')
            y(i + offset) = scaled_variance(work(i:i + n - 1))
         case('sd')
            y(i + offset) = sqrt(scaled_variance(work(i:i + n - 1)))
         end select
      end do
   end function

   pure function roll_var(x, n, trim, unbiased, na_rm) result(y)
      real(dp), intent(in) :: x(:) !! Chronological observations with NaNs propagated window-wise.
      integer, intent(in), optional :: n !! Complete-window width; default nine.
      logical, intent(in), optional :: trim !! Remove n-1 leading unavailable values; default true.
      logical, intent(in), optional :: unbiased !! Divide by n-1 when true (default), otherwise n.
      logical, intent(in), optional :: na_rm !! Remove all NaNs before rolling; default false.
      real(dp), allocatable :: y(:)
      integer :: width
      width = 9
      if (present(n)) width = n
      y = roll_named(x, width, 'var', trim, na_rm)
      if (present(unbiased)) then
         if (.not. unbiased .and. width > 0) y = y * real(width - 1, dp) / real(width, dp)
      end if
   end function

   pure function sma(x, n) result(y)
      real(dp), intent(in) :: x(:) !! Chronological observations with window-local NaN propagation.
      integer, intent(in), optional :: n !! Positive complete-window width; default five.
      real(dp), allocatable :: y(:)
      integer :: width
      width = 5
      if (present(n)) width = n
      y = roll_named(x, width, 'mean')
   end function

   pure function ema_ta(x, lambda, startup) result(y)
      use ttr_base, only: ttr_ema => ema
      real(dp), intent(in) :: x(:) !! Chronological observations; initial seed is a future startup-window mean.
      real(dp), intent(in), optional :: lambda !! Decay in [0,1), or period >=1 converted to 2/(period+1); default 0.1.
      integer, intent(in), optional :: startup !! Seed width, zero means floor(2/decay); default zero.
      real(dp) :: y(size(x)), work(size(x)), decay
      integer :: seed, i
      y = nan_value()
      if (size(x) == 0) return
      decay = 0.1_dp
      if (present(lambda)) decay = lambda
      if (.not. ieee_is_finite(decay)) return
      if (decay < 0.0_dp) return
      if (decay == 0.0_dp) then
         y = scaled_mean(x)
         return
      end if
      if (decay >= 1.0_dp) decay = 2.0_dp / (decay + 1.0_dp)
      seed = 0
      if (present(startup)) seed = startup
      if (seed < 0) return
      if (seed == 0) then
         if (2.0_dp / decay > real(size(x), dp) + 1.0_dp) return
         seed = floor(2.0_dp / decay)
      end if
      if (seed < 1 .or. seed > size(x)) return
      work = x
      work(1) = scaled_mean(x(:seed))
      if (all(ieee_is_finite(work))) then
         y = ttr_ema(work, 1, ratio = decay)
      else
         y(1) = work(1)
         do i = 2, size(x)
            y(i) = decay * x(i) + (1.0_dp - decay) * y(i - 1)
         end do
      end if
   end function

   pure function ewma(x, lambda, startup) result(y)
      real(dp), intent(in) :: x(:) !! Chronological observations with fTrading's startup-mean seed.
      real(dp), intent(in) :: lambda !! Decay in [0,1), or period >=1 converted to 2/(period+1).
      integer, intent(in), optional :: startup !! Seed width; zero/default selects floor(2/decay).
      real(dp) :: y(size(x))
      y = ema_ta(x, lambda, startup)
   end function

   pure subroutine trade_signals(positions, signals, indices)
      real(dp), intent(in) :: positions(:) !! Chronological trading positions, typically -1 or +1.
      real(dp), allocatable, intent(out) :: signals(:) !! Half-changes whose absolute magnitude equals one.
      integer, allocatable, intent(out) :: indices(:) !! One-based positions of returned trade signals.
      real(dp) :: changes(max(0, size(positions) - 1))
      integer :: i, k
      changes = (positions(2:) - positions(:max(0, size(positions) - 1))) / 2.0_dp
      allocate(signals(count(abs(changes) == 1.0_dp)), indices(count(abs(changes) == 1.0_dp)))
      k = 0
      do i = 1, size(changes)
         if (abs(changes(i)) == 1.0_dp) then
            k = k + 1
            signals(k) = changes(i)
            indices(k) = i + 1
         end if
      end do
   end subroutine

   pure function trade_lengths(times) result(lengths)
      real(dp), intent(in) :: times(:) !! Signal timestamps in a consistent numeric unit, usually days.
      real(dp) :: lengths(max(0, size(times) - 1))
      lengths = times(2:) - times(:max(0, size(times) - 1))
   end function

   pure function hit_rate(returns, positions) result(rate)
      real(dp), intent(in) :: returns(:) !! Chronological returns; NaN-derived indicators are omitted.
      real(dp), intent(in) :: positions(:) !! Aligned signed positions; lengths must match returns.
      real(dp) :: rate, indicators(size(returns)), scaled
      integer :: i, base
      rate = nan_value()
      if (size(returns) /= size(positions)) return
      do i = 1, size(returns)
         indicators(i) = nan_value()
         if (ieee_is_nan(returns(i))) cycle
         if (returns(i) == 0.0_dp) then
            indicators(i) = (positions(i) * 0.0_dp + 1.0_dp) / 2.0_dp
         else
            indicators(i) = (positions(i) * sign(1.0_dp, returns(i)) + 1.0_dp) / 2.0_dp
         end if
         if (ieee_is_nan(positions(i))) indicators(i) = nan_value()
      end do
      rate = r_mean(indicators, na_rm = .true.)
      if (.not. ieee_is_finite(rate)) return
      scaled = abs(rate) * 1000.0_dp
      if (scaled > real(huge(base) - 1, dp)) return
      base = floor(scaled)
      if (scaled - real(base, dp) > 0.5_dp) then
         base = base + 1
      else if (scaled - real(base, dp) == 0.5_dp) then
         if (mod(base, 2) /= 0) base = base + 1
      end if
      rate = sign(real(base, dp) / 1000.0_dp, rate)
   end function
end module ftrading_base
