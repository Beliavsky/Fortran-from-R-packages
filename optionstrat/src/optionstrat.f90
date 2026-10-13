! SPDX-License-Identifier: GPL-3.0-only
! Derived from optionstrat 1.4.1 by John T. Buynak. See NOTICE.md.
module optionstrat
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_distributions, only : r_pnorm, r_qnorm, r_dnorm
   use r_transforms, only : r_log1p
   implicit none
   private
   public :: dp, option_result, dual_option_result, price_range, vertical_result, double_vertical_result
   public :: callpremium, putpremium, calldelta, putdelta, calltheta, puttheta
   public :: callrho, putrho, optiongamma, optionvega, calleval, puteval, opteval
   public :: callgreek, putgreek, lambda, iv_calc, implied_volatility, r_cont, tdiff
   public :: prob_above, prob_below, prob_btwn, prob_above_quantile, prob_below_quantile, prob_btwn_quantile
   public :: vertical, dv, calltheta_standard, puttheta_standard
   public :: prob_above_quantile_standard, prob_btwn_quantile_standard

   type :: option_result
      real(dp) :: premium, delta, gamma, vega, theta, rho
   end type option_result
   type :: dual_option_result
      type(option_result) :: call, put
   end type dual_option_result
   type :: price_range
      real(dp) :: probability, lower_change, upper_change, lower_price, upper_price
   end type price_range
   type :: vertical_result
      real(dp) :: spot, short_strike, long_strike, max_profit, max_loss, breakeven
      real(dp) :: prob_be, prob_max_profit, prob_max_loss, initial_dc
   end type vertical_result
   type :: double_vertical_result
      real(dp) :: spot, long_put, short_put, short_call, long_call
      real(dp) :: lower_be, higher_be, prob_be, prob_max_profit, prob_max_loss
   end type double_vertical_result
contains
   pure real(dp) function nan_value() result(value)
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure elemental logical function valid_option(s, x, sigma, t, r, d) result(valid)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in) :: d !! Annual continuously compounded dividend yield, default zero.
      valid = all(ieee_is_finite([s, x, sigma, t, r, d]))
      valid = valid .and. s > 0.0_dp .and. x > 0.0_dp .and. sigma >= 0.0_dp .and. t >= 0.0_dp
      if (.not. valid) return
      ! Reject discount-factor overflow before evaluating exponentials.
      valid = -r * t < log(huge(1.0_dp)) - log(max(x, 1.0_dp))
      valid = valid .and. (-d) * t < log(huge(1.0_dp)) - log(max(s, 1.0_dp))
   end function valid_option

   pure elemental function evaluate(s, x, sigma, t, r, d, is_call) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in) :: d !! Annual continuously compounded dividend yield, default zero.
      logical, intent(in) :: is_call !! True selects call analytics, false selects put analytics.
      type(option_result) :: value
      real(dp) :: te, sd, xd, z1, z2, n1, n2, density, root_t, width, m, q
      q = nan_value()
      value = option_result(q, q, q, q, q, q)
      if (.not. valid_option(s, x, sigma, t, r, d)) return
      te = t
      if (te == 0.0_dp) te = 1.0e-200_dp
      root_t = sqrt(te)
      sd = s * exp(-d * te)
      xd = x * exp(-r * te)
      m = (log(s) - log(x)) + (r - d) * te
      width = sigma * root_t
      if (.not. ieee_is_finite(width)) return
      if (width == 0.0_dp) then
         n1 = 0.5_dp
         if (m > 0.0_dp) n1 = 1.0_dp
         if (m < 0.0_dp) n1 = 0.0_dp
         n2 = n1
         density = 0.0_dp
         if (m == 0.0_dp) density = r_dnorm(0.0_dp)
         value%gamma = 0.0_dp
         if (m == 0.0_dp) value%gamma = q
      else
         z1 = m / width + 0.5_dp * width
         z2 = m / width - 0.5_dp * width
         n1 = r_pnorm(z1)
         n2 = r_pnorm(z2)
         density = r_dnorm(z1)
         value%gamma = exp(-d * te) * density / s / width
      end if
      value%vega = 0.01_dp * sd * root_t * density
      if (is_call) then
         value%premium = max(0.0_dp, sd * n1 - xd * n2)
         value%delta = exp(-d * te) * n1
         value%rho = 0.01_dp * xd * te * n2
         ! Preserve upstream's extra 1 / te on the rate/dividend terms.
         value%theta = (-0.5_dp * sd * width * density - r * xd * n2 + d * sd * n1) / 365.0_dp / te
      else
         if (width > 0.0_dp) then
            n1 = r_pnorm(-z1)
            n2 = r_pnorm(-z2)
         else
            n1 = 1.0_dp - n1
            n2 = 1.0_dp - n2
         end if
         value%premium = max(0.0_dp, xd * n2 - sd * n1)
         value%delta = -exp(-d * te) * n1
         value%rho = -0.01_dp * xd * te * n2
         value%theta = (-0.5_dp * sd * width * density + r * xd * n2 - d * sd * n1) / 365.0_dp / te
      end if
   end function evaluate

   pure elemental function calleval(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: value
      real(dp) :: yield
      yield = 0.0_dp
      if (present(d)) yield = d
      value = evaluate(s, x, sigma, t, r, yield, .true.)
   end function calleval

   pure elemental function puteval(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: value
      real(dp) :: yield
      yield = 0.0_dp
      if (present(d)) yield = d
      value = evaluate(s, x, sigma, t, r, yield, .false.)
   end function puteval

   pure elemental real(dp) function callpremium(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      value = analytics%premium
   end function callpremium

   pure elemental real(dp) function putpremium(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = puteval(s, x, sigma, t, r, d)
      value = analytics%premium
   end function putpremium

   pure elemental real(dp) function calldelta(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      value = analytics%delta
   end function calldelta

   pure elemental real(dp) function putdelta(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = puteval(s, x, sigma, t, r, d)
      value = analytics%delta
   end function putdelta

   pure elemental real(dp) function calltheta(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      value = analytics%theta
   end function calltheta

   pure elemental real(dp) function puttheta(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = puteval(s, x, sigma, t, r, d)
      value = analytics%theta
   end function puttheta

   pure elemental real(dp) function callrho(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      value = analytics%rho
   end function callrho

   pure elemental real(dp) function putrho(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = puteval(s, x, sigma, t, r, d)
      value = analytics%rho
   end function putrho

   pure elemental real(dp) function optiongamma(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      value = analytics%gamma
   end function optiongamma

   pure elemental real(dp) function optionvega(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      value = analytics%vega
   end function optionvega

   pure elemental function opteval(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(dual_option_result) :: value
      value%call = calleval(s, x, sigma, t, r, d)
      value%put = puteval(s, x, sigma, t, r, d)
   end function opteval

   pure elemental real(dp) function callgreek(greek, s, x, sigma, t, r, d) result(value)
      character(len=*), intent(in) :: greek !! Lowercase premium, delta, gamma, vega, theta, or rho.
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = calleval(s, x, sigma, t, r, d)
      select case (greek)
      case ('premium')
         value = analytics%premium
      case ('delta')
         value = analytics%delta
      case ('gamma')
         value = analytics%gamma
      case ('vega')
         value = analytics%vega
      case ('theta')
         value = analytics%theta
      case ('rho')
         value = analytics%rho
      case default
         value = nan_value()
      end select
   end function callgreek

   pure elemental real(dp) function putgreek(greek, s, x, sigma, t, r, d) result(value)
      character(len=*), intent(in) :: greek !! Lowercase premium, delta, gamma, vega, theta, or rho.
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      analytics = puteval(s, x, sigma, t, r, d)
      select case (greek)
      case ('premium')
         value = analytics%premium
      case ('delta')
         value = analytics%delta
      case ('gamma')
         value = analytics%gamma
      case ('vega')
         value = analytics%vega
      case ('theta')
         value = analytics%theta
      case ('rho')
         value = analytics%rho
      case default
         value = nan_value()
      end select
   end function putgreek

   pure elemental real(dp) function lambda(type, s, x, sigma, t, r, d) result(value)
      character(len=*), intent(in) :: type !! Option type, exactly call or put.
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal, e.g. 0.20.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      value = nan_value()
      select case (type)
      case ('call')
         analytics = calleval(s, x, sigma, t, r, d)
      case ('put')
         analytics = puteval(s, x, sigma, t, r, d)
      case default
         return
      end select
      if (analytics%premium > 0.0_dp) value = analytics%delta * (s / analytics%premium)
   end function lambda

   pure elemental real(dp) function iv_calc(type, price, s, x, t, r, d) result(value)
      character(len=*), intent(in) :: type !! Option type, exactly call or put.
      real(dp), intent(in) :: price !! Finite nonnegative observed option premium in currency units.
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      real(dp) :: yield, distance, best, candidate, model, te
      integer :: k
      value = nan_value()
      yield = 0.0_dp
      if (present(d)) yield = d
      if (.not. valid_option(s, x, 0.0_dp, t, r, yield)) return
      if (.not. ieee_is_finite(price) .or. price < 0.0_dp) return
      if (type /= 'call' .and. type /= 'put') return
      ! The R grid contains NaN at sigma = 0 when the forward log-moneyness is zero.
      ! Its min/match sequence then selects the first point, regardless of price.
      te = t
      if (te == 0.0_dp) te = 1.0e-200_dp
      if (log(s / x) + (r - yield) * te == 0.0_dp) then
         value = 0.0_dp
         return
      end if
      best = huge(1.0_dp)
      do k = 0, 5000
         candidate = real(k, dp) * 0.001_dp
         if (type == 'call') then
            model = callpremium(s, x, candidate, t, r, yield)
         else
            model = putpremium(s, x, candidate, t, r, yield)
         end if
         distance = abs(price - model)
         if (distance < best) then
            best = distance
            value = candidate
         end if
      end do
   end function iv_calc

   pure elemental real(dp) function implied_volatility(type, price, s, x, t, r, d) result(value)
      character(len=*), intent(in) :: type !! Option type, exactly call or put.
      real(dp), intent(in) :: price !! Finite nonnegative observed option premium in currency units.
      real(dp), intent(in) :: s !! Positive spot price in currency units.
      real(dp), intent(in) :: x !! Positive strike price in the same currency as spot.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the upstream 1e-200 convention.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      real(dp) :: yield, low, high, mid, model, lower_price, upper_price, scale
      integer :: iteration
      value = nan_value()
      yield = 0.0_dp
      if (present(d)) yield = d
      if (.not. valid_option(s, x, 0.0_dp, t, r, yield)) return
      if (.not. ieee_is_finite(price) .or. t <= 0.0_dp) return
      select case (type)
      case ('call')
         lower_price = callpremium(s, x, 0.0_dp, t, r, yield)
         upper_price = s * exp(-yield * t)
      case ('put')
         lower_price = putpremium(s, x, 0.0_dp, t, r, yield)
         upper_price = x * exp(-r * t)
      case default
         return
      end select
      if (price < lower_price .or. price >= upper_price) return
      if (price == lower_price) then
         value = 0.0_dp
         return
      end if
      low = 0.0_dp
      high = 1.0_dp
      ! A total standard deviation of 64 saturates both normal tails.
      scale = max(1.0_dp, 64.0_dp / sqrt(t))
      do
         if (type == 'call') then
            model = callpremium(s, x, high, t, r, yield)
         else
            model = putpremium(s, x, high, t, r, yield)
         end if
         if (model >= price) exit
         if (high >= scale) return
         high = min(2.0_dp * high, scale)
      end do
      do iteration = 1, 160
         mid = low + 0.5_dp * (high - low)
         if (mid == low .or. mid == high) exit
         if (type == 'call') then
            model = callpremium(s, x, mid, t, r, yield)
         else
            model = putpremium(s, x, mid, t, r, yield)
         end if
         if (model < price) then
            low = mid
         else
            high = mid
         end if
      end do
      value = low + 0.5_dp * (high - low)
   end function implied_volatility

   pure elemental real(dp) function r_cont(r, n) result(value)
      real(dp), intent(in) :: r !! Nominal annual rate as a decimal, with 1 + r / n positive.
      real(dp), intent(in) :: n !! Positive compounding frequency per year, not restricted to integers.
      value = nan_value()
      if (.not. all(ieee_is_finite([r, n]))) return
      if (n <= 0.0_dp) return
      if (r / n <= -1.0_dp) return
      value = n * r_log1p(r / n)
   end function r_cont

   pure elemental integer function civil_day(date) result(day)
      character(len=*), intent(in) :: date !! ISO Gregorian date YYYY-MM-DD in years 0001 through 9999.
      integer :: year, month, day_month, j, digit, days(12), ios
      day = -1
      if (len_trim(date) /= 10) return
      if (date(5:5) /= '-' .or. date(8:8) /= '-') return
      do j = 1, 10
         if (j == 5 .or. j == 8) cycle
         digit = iachar(date(j:j)) - iachar('0')
         if (digit < 0 .or. digit > 9) return
      end do
      read(date(1:4), '(i4)', iostat = ios) year
      if (ios /= 0) return
      read(date(6:7), '(i2)', iostat = ios) month
      if (ios /= 0) return
      read(date(9:10), '(i2)', iostat = ios) day_month
      if (ios /= 0) return
      if (year < 1 .or. month < 1 .or. month > 12) return
      days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      if (mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)) days(2) = 29
      if (day_month < 1 .or. day_month > days(month)) return
      year = year - 1
      day = 365 * year + year / 4 - year / 100 + year / 400 + sum(days(:month - 1)) + day_month
   end function civil_day

   pure elemental real(dp) function tdiff(date1, date2, period) result(value)
      character(len=*), intent(in) :: date1 !! Earlier ISO date YYYY-MM-DD, without a timezone or time of day.
      character(len=*), intent(in) :: date2 !! Later ISO date YYYY-MM-DD; reversed dates give a negative difference.
      character(len=*), intent(in) :: period !! Either days or years; a year always means 365 days.
      integer :: first, last
      value = nan_value()
      first = civil_day(date1)
      last = civil_day(date2)
      if (first < 0 .or. last < 0) return
      if (period /= 'days' .and. period /= 'years') return
      value = real(last - first, dp)
      if (period == 'years') value = value / 365.0_dp
   end function tdiff

   pure elemental real(dp) function horizon_scale(asd, dsd, dte, tradedays) result(value)
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      real(dp) :: annual, daily, days, trading
      value = nan_value()
      annual = 0.0_dp
      daily = 0.0_dp
      days = 0.0_dp
      trading = 262.0_dp
      if (present(asd)) annual = asd
      if (present(dsd)) daily = dsd
      if (present(dte)) days = dte
      if (present(tradedays)) trading = tradedays
      if (.not. all(ieee_is_finite([annual, daily, days, trading]))) return
      if (min(annual, daily, days) < 0.0_dp .or. trading <= 0.0_dp) return
      if (annual == 0.0_dp .and. daily == 0.0_dp .and. days == 0.0_dp) return
      value = annual
      if (daily == 0.0_dp .and. annual > 0.0_dp) then
         value = annual / sqrt(trading) * sqrt(days)
      else if (annual == 0.0_dp .and. daily > 0.0_dp .and. days > 0.0_dp) then
         value = daily * sqrt(days)
      end if
   end function horizon_scale

   pure elemental real(dp) function return_cdf(spot, level, location, scale, above) result(value)
      real(dp), intent(in) :: spot !! Strictly positive current asset price in currency units.
      real(dp), intent(in) :: level !! Finite comparison price in the same currency as spot.
      real(dp), intent(in) :: location !! Finite arithmetic-return normal location as a decimal.
      real(dp), intent(in) :: scale !! Nonnegative normal return standard deviation over the horizon.
      logical, intent(in) :: above !! True requests the strict upper tail rather than the lower CDF.
      real(dp) :: change
      value = nan_value()
      if (.not. all(ieee_is_finite([spot, level, location, scale]))) return
      if (spot <= 0.0_dp .or. scale < 0.0_dp) return
      change = (level - spot) / spot
      if (scale == 0.0_dp) then
         value = 0.0_dp
         if (change >= location) value = 1.0_dp
         if (above) value = 1.0_dp - value
      else
         value = r_pnorm(change, location, scale, lower_tail = .not. above)
      end if
   end function return_cdf

   pure elemental real(dp) function prob_below(spot, upper, asd, dsd, dte, mean, tradedays) result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: upper !! Upper price threshold in currency units.
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Normal arithmetic-return location as a decimal, default zero; not multiplied by horizon.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      real(dp) :: location, scale
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      value = return_cdf(spot, upper, location, scale, .false.)
   end function prob_below

   pure elemental real(dp) function prob_above(spot, lower, asd, dsd, dte, mean, tradedays) result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: lower !! Lower price threshold in currency units.
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Normal arithmetic-return location as a decimal, default zero; not multiplied by horizon.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      real(dp) :: location, scale
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      value = return_cdf(spot, lower, location, scale, .true.)
   end function prob_above

   pure elemental real(dp) function prob_btwn(spot, lower, upper, asd, dsd, dte, mean, tradedays) result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: lower !! Lower price threshold in currency units.
      real(dp), intent(in) :: upper !! Upper price threshold in currency units.
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Normal arithmetic-return location as a decimal, default zero; not multiplied by horizon.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      value = prob_below(spot, upper, asd, dsd, dte, mean, tradedays) &
         - prob_below(spot, lower, asd, dsd, dte, mean, tradedays)
   end function prob_btwn

   pure elemental real(dp) function prob_below_quantile(spot, p, asd, dsd, dte, mean, tradedays) result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: p !! Requested probability in the closed interval [0, 1].
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Normal arithmetic-return location as a decimal, default zero; not multiplied by horizon.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      real(dp) :: location, scale
      value = nan_value()
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      if (.not. all(ieee_is_finite([spot, p, location, scale]))) return
      if (spot <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) return
      value = spot * (1.0_dp + r_qnorm(p, location, scale))
   end function prob_below_quantile

   pure elemental real(dp) function prob_above_quantile(spot, p, asd, dsd, dte, mean, tradedays) result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: p !! Requested probability in the closed interval [0, 1].
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Normal arithmetic-return location as a decimal, default zero; not multiplied by horizon.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      real(dp) :: location, scale
      value = nan_value()
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      if (.not. all(ieee_is_finite([spot, p, location, scale]))) return
      if (spot <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) return
      value = spot * (1.0_dp - r_qnorm(p, location, scale))
   end function prob_above_quantile

   pure elemental function prob_btwn_quantile(spot, p, asd, dsd, dte, mean, tradedays) result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: p !! Requested probability in the closed interval [0, 1].
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero; both volatilities positive means use asd as horizon scale.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero; annual conversion uses tradedays.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Normal arithmetic-return location as a decimal, default zero; not multiplied by horizon.
      real(dp), intent(in), optional :: tradedays !! Positive trading days per year, default 262.
      type(price_range) :: value
      real(dp) :: location, scale, q
      q = nan_value()
      value = price_range(q, q, q, q, q)
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      if (.not. all(ieee_is_finite([spot, p, location, scale]))) return
      if (spot <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) return
      value%probability = p
      value%lower_change = r_qnorm(0.5_dp * (1.0_dp - p), location, scale)
      value%upper_change = -value%lower_change
      value%lower_price = spot * (1.0_dp + value%lower_change)
      value%upper_price = spot * (1.0_dp + value%upper_change)
   end function prob_btwn_quantile

   pure elemental function vertical(options, s, x1, x2, t, r, sigma, sigma2, vol, d) result(value)
      character(len=*), intent(in) :: options !! Option family, exactly call or put.
      real(dp), intent(in) :: s !! Positive underlying spot price in currency units.
      real(dp), intent(in) :: x1 !! Positive short strike for vertical, or long put strike for dv.
      real(dp), intent(in) :: x2 !! Positive long strike for vertical, or short put strike for dv.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero retains the upstream tiny positive tenor.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in) :: sigma !! Nonnegative annual decimal volatility of the first option leg.
      real(dp), intent(in), optional :: sigma2 !! Annual decimal volatility of the second leg, default sigma.
      real(dp), intent(in), optional :: vol !! Nonnegative annual volatility for probability calculations, default sigma.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield as a decimal, default zero.
      type(vertical_result) :: value
      real(dp) :: second_vol, probability_vol, yield, c1, c2, q, days, daily, te
      q = nan_value()
      value = vertical_result(q, q, q, q, q, q, q, q, q, q)
      second_vol = sigma
      probability_vol = sigma
      yield = 0.0_dp
      if (present(sigma2)) second_vol = sigma2
      if (present(vol)) probability_vol = vol
      if (present(d)) yield = d
      if (.not. valid_option(s, x1, sigma, t, r, yield)) return
      if (.not. valid_option(s, x2, second_vol, t, r, yield)) return
      if (.not. ieee_is_finite(probability_vol) .or. probability_vol < 0.0_dp) return
      if (options /= 'call' .and. options /= 'put') return
      if (x1 == x2) return
      te = t
      if (te == 0.0_dp) te = 1.0e-200_dp
      days = te * 365.0_dp
      daily = probability_vol / 16.0_dp
      if (options == 'call') then
         c1 = callpremium(s, x1, sigma, t, r, yield)
         c2 = callpremium(s, x2, second_vol, t, r, yield)
      else
         c1 = putpremium(s, x1, sigma, t, r, yield)
         c2 = putpremium(s, x2, second_vol, t, r, yield)
      end if
      value%spot = s
      value%short_strike = x1
      value%long_strike = x2
      value%initial_dc = c1 - c2
      if (options == 'call' .and. x1 > x2) then
         value%max_profit = x1 - x2 + c1 - c2
         value%max_loss = c1 - c2
         value%breakeven = x2 + abs(c1 - c2)
         value%prob_be = prob_above(s, value%breakeven, dsd = daily, dte = days)
         value%prob_max_profit = prob_above(s, x1, dsd = daily, dte = days)
         value%prob_max_loss = prob_below(s, x2, dsd = daily, dte = days)
      else if (options == 'call') then
         value%max_profit = c1 - c2
         value%max_loss = x1 - x2 + c1 - c2
         value%breakeven = x1 + c1 - c2
         value%prob_be = prob_below(s, value%breakeven, dsd = daily, dte = days)
         value%prob_max_profit = prob_below(s, x1, dsd = daily, dte = days)
         value%prob_max_loss = prob_above(s, x2, dsd = daily, dte = days)
      else if (x1 < x2) then
         value%max_profit = x2 - x1 + c1 - c2
         value%max_loss = c1 - c2
         value%breakeven = x2 - abs(c1 - c2)
         value%prob_be = prob_below(s, value%breakeven, dsd = daily, dte = days)
         value%prob_max_profit = prob_below(s, x1, dsd = daily, dte = days)
         value%prob_max_loss = prob_above(s, x2, dsd = daily, dte = days)
      else
         value%max_profit = c1 - c2
         value%max_loss = x2 - x1 + c1 - c2
         value%breakeven = x1 - c1 + c2
         value%prob_be = prob_above(s, value%breakeven, dsd = daily, dte = days)
         value%prob_max_profit = prob_above(s, x1, dsd = daily, dte = days)
         value%prob_max_loss = prob_below(s, x2, dsd = daily, dte = days)
      end if
   end function vertical

   pure elemental function dv(s, x1, x2, x3, x4, t, r, sigma, sigma2, sigma3, sigma4, vol, d) result(value)
      real(dp), intent(in) :: s !! Positive underlying spot price in currency units.
      real(dp), intent(in) :: x1 !! Positive short strike for vertical, or long put strike for dv.
      real(dp), intent(in) :: x2 !! Positive long strike for vertical, or short put strike for dv.
      real(dp), intent(in) :: x3 !! Positive short call strike; dv requires x1 < x2 <= x3 < x4.
      real(dp), intent(in) :: x4 !! Positive long call strike; dv requires x1 < x2 <= x3 < x4.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero retains the upstream tiny positive tenor.
      real(dp), intent(in) :: r !! Finite annual continuously compounded rate as a decimal.
      real(dp), intent(in) :: sigma !! Nonnegative annual decimal volatility of the first option leg.
      real(dp), intent(in), optional :: sigma2 !! Annual decimal volatility of the second leg, default sigma.
      real(dp), intent(in), optional :: sigma3 !! Annual decimal volatility of the third leg, default sigma.
      real(dp), intent(in), optional :: sigma4 !! Annual decimal volatility of the fourth leg, default sigma.
      real(dp), intent(in), optional :: vol !! Nonnegative annual volatility for probability calculations, default sigma.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield as a decimal, default zero.
      type(double_vertical_result) :: value
      real(dp) :: v2, v3, v4, probability_vol, yield, o1, o2, o3, o4, q, days, daily, credit
      q = nan_value()
      value = double_vertical_result(q, q, q, q, q, q, q, q, q, q)
      v2 = sigma
      v3 = sigma
      v4 = sigma
      probability_vol = sigma
      yield = 0.0_dp
      if (present(sigma2)) v2 = sigma2
      if (present(sigma3)) v3 = sigma3
      if (present(sigma4)) v4 = sigma4
      if (present(vol)) probability_vol = vol
      if (present(d)) yield = d
      if (.not. valid_option(s, x1, sigma, t, r, yield)) return
      if (.not. valid_option(s, x2, v2, t, r, yield)) return
      if (.not. valid_option(s, x3, v3, t, r, yield)) return
      if (.not. valid_option(s, x4, v4, t, r, yield)) return
      if (.not. ieee_is_finite(probability_vol) .or. probability_vol < 0.0_dp) return
      if (x1 >= x2 .or. x2 > x3 .or. x3 >= x4) return
      days = max(t, 1.0e-200_dp) * 365.0_dp
      daily = probability_vol / 16.0_dp
      o1 = putpremium(s, x1, sigma, t, r, yield)
      o2 = putpremium(s, x2, v2, t, r, yield)
      o3 = callpremium(s, x3, v3, t, r, yield)
      o4 = callpremium(s, x4, v4, t, r, yield)
      credit = (o2 + o3) - (o1 + o4)
      value%spot = s
      value%long_put = x1
      value%short_put = x2
      value%short_call = x3
      value%long_call = x4
      value%lower_be = x2 - credit
      value%higher_be = x3 + credit
      value%prob_be = prob_btwn(s, value%lower_be, value%higher_be, dsd = daily, dte = days)
      value%prob_max_profit = prob_btwn(s, x2, x3, dsd = daily, dte = days)
      value%prob_max_loss = prob_below(s, x1, dsd = daily, dte = days) &
         + prob_above(s, x4, dsd = daily, dte = days)
   end function dv


   pure elemental real(dp) function calltheta_standard(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive underlying spot price in currency units.
      real(dp), intent(in) :: x !! Positive option strike in currency units.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the tiny-tenor convention.
      real(dp), intent(in) :: r !! Annual continuously compounded decimal risk-free rate.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      real(dp) :: yield, diffusion
      yield = 0.0_dp
      if (present(d)) yield = d
      analytics = calleval(s, x, sigma, t, r, yield)
      diffusion = 0.0_dp
      if (sigma > 0.0_dp) diffusion = -0.5_dp * s * (sigma * (sigma * (s * analytics%gamma)))
      value = (diffusion + yield * s * analytics%delta &
         - r * (s * analytics%delta - analytics%premium)) / 365.0_dp
   end function calltheta_standard

   pure elemental real(dp) function puttheta_standard(s, x, sigma, t, r, d) result(value)
      real(dp), intent(in) :: s !! Positive underlying spot price in currency units.
      real(dp), intent(in) :: x !! Positive option strike in currency units.
      real(dp), intent(in) :: sigma !! Nonnegative annual volatility as a decimal.
      real(dp), intent(in) :: t !! Nonnegative tenor in years; zero uses the tiny-tenor convention.
      real(dp), intent(in) :: r !! Annual continuously compounded decimal risk-free rate.
      real(dp), intent(in), optional :: d !! Annual continuously compounded dividend yield, default zero.
      type(option_result) :: analytics
      real(dp) :: yield, diffusion
      yield = 0.0_dp
      if (present(d)) yield = d
      analytics = puteval(s, x, sigma, t, r, yield)
      diffusion = 0.0_dp
      if (sigma > 0.0_dp) diffusion = -0.5_dp * s * (sigma * (sigma * (s * analytics%gamma)))
      value = (diffusion + yield * s * analytics%delta &
         - r * (s * analytics%delta - analytics%premium)) / 365.0_dp
   end function puttheta_standard

   pure elemental real(dp) function prob_above_quantile_standard(spot, p, asd, dsd, dte, mean, tradedays) &
      result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: p !! Requested probability in [0, 1].
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Arithmetic-return location as a decimal, default zero.
      real(dp), intent(in), optional :: tradedays !! Positive number of trading days per year, default 262.
      real(dp) :: location, scale, q
      q = nan_value()
      value = q
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      if (.not. all(ieee_is_finite([spot, p, location, scale]))) return
      if (spot <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) return
      value = spot * (1.0_dp + r_qnorm(p, location, scale, lower_tail = .false.))
   end function prob_above_quantile_standard

   pure elemental function prob_btwn_quantile_standard(spot, p, asd, dsd, dte, mean, tradedays) &
      result(value)
      real(dp), intent(in) :: spot !! Positive current asset price in currency units.
      real(dp), intent(in) :: p !! Requested probability in [0, 1].
      real(dp), intent(in), optional :: asd !! Annual decimal volatility, default zero.
      real(dp), intent(in), optional :: dsd !! Daily decimal volatility, default zero.
      real(dp), intent(in), optional :: dte !! Nonnegative horizon in days, default zero.
      real(dp), intent(in), optional :: mean !! Arithmetic-return location as a decimal, default zero.
      real(dp), intent(in), optional :: tradedays !! Positive number of trading days per year, default 262.
      type(price_range) :: value
      real(dp) :: location, scale, q
      q = nan_value()
      value = price_range(q, q, q, q, q)
      location = 0.0_dp
      if (present(mean)) location = mean
      scale = horizon_scale(asd, dsd, dte, tradedays)
      if (.not. all(ieee_is_finite([spot, p, location, scale]))) return
      if (spot <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) return
      value%probability = p
      value%lower_change = r_qnorm(0.5_dp * (1.0_dp - p), location, scale)
      value%upper_change = r_qnorm(0.5_dp * (1.0_dp - p), location, scale, lower_tail = .false.)
      value%lower_price = spot * (1.0_dp + value%lower_change)
      value%upper_price = spot * (1.0_dp + value%upper_change)
   end function prob_btwn_quantile_standard
end module optionstrat
