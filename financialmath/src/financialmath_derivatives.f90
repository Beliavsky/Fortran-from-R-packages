! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module financialmath_derivatives
   use financialmath_kinds, only : dp
   use financialmath_types, only : option_order1_t, payoff_table_t, forward_result_t
   use financialmath_math, only : normal_cdf, normal_pdf
   implicit none
   private
   public :: bls_order1, black_scholes_call, black_scholes_put
   public :: option_call, option_put, covered_call, covered_put, protective_put
   public :: straddle, straddle_bls, strangle, strangle_bls
   public :: bull_call, bull_call_bls, bear_call, bear_call_bls
   public :: collar, collar_bls, butterfly_spread, butterfly_spread_bls
   public :: forward_contract, forward_prepaid, forward_price, prepaid_forward_price

contains

   pure real(dp) function black_scholes_call(spot, strike, rate, time, volatility, dividend_yield) result(value)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      real(dp), intent(in), optional :: dividend_yield
      real(dp) :: q, d1, d2
      q = 0.0_dp
      if (present(dividend_yield)) q = dividend_yield
      if (time <= 0.0_dp) then
         value = max(spot-strike, 0.0_dp)
      else if (volatility <= 0.0_dp) then
         value = max(spot*exp(-q*time)-strike*exp(-rate*time), 0.0_dp)
      else
         d1 = (log(spot/strike)+(rate-q+0.5_dp*volatility*volatility)*time)/(volatility*sqrt(time))
         d2 = d1-volatility*sqrt(time)
         value = spot*exp(-q*time)*normal_cdf(d1)-strike*exp(-rate*time)*normal_cdf(d2)
      end if
   end function black_scholes_call

   pure real(dp) function black_scholes_put(spot, strike, rate, time, volatility, dividend_yield) result(value)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      real(dp), intent(in), optional :: dividend_yield
      real(dp) :: q
      q = 0.0_dp
      if (present(dividend_yield)) q = dividend_yield
      value = black_scholes_call(spot, strike, rate, time, volatility, q) - &
         spot*exp(-q*time)+strike*exp(-rate*time)
   end function black_scholes_put

   function bls_order1(spot, strike, rate, time, volatility, dividend_yield) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility, dividend_yield
      type(option_order1_t) :: out
      real(dp) :: d1, d2, root_t

      out%status%ok = .false.
      if (spot <= 0.0_dp .or. strike <= 0.0_dp .or. time <= 0.0_dp .or. volatility <= 0.0_dp) then
         out%status%message = 'spot, strike, time, and volatility must be positive'
         return
      end if
      root_t = sqrt(time)
      d1 = (log(spot/strike)+(rate-dividend_yield+0.5_dp*volatility*volatility)*time)/(volatility*root_t)
      d2 = d1-volatility*root_t
      out%call_price = spot*exp(-dividend_yield*time)*normal_cdf(d1) - &
         strike*exp(-rate*time)*normal_cdf(d2)
      out%put_price = -spot*exp(-dividend_yield*time)*normal_cdf(-d1) + &
         strike*exp(-rate*time)*normal_cdf(-d2)
      out%call_delta = exp(-dividend_yield*time)*normal_cdf(d1)
      out%put_delta = exp(-dividend_yield*time)*(normal_cdf(d1)-1.0_dp)
      out%call_theta = -volatility*spot*exp(-dividend_yield*time)*normal_pdf(d1)/(2.0_dp*root_t) + &
         dividend_yield*spot*normal_cdf(-d1)*exp(-dividend_yield*time) - &
         rate*strike*exp(-rate*time)*normal_cdf(d2)
      out%put_theta = -volatility*spot*exp(-dividend_yield*time)*normal_pdf(d1)/(2.0_dp*root_t) - &
         dividend_yield*spot*normal_cdf(d1)*exp(-dividend_yield*time) + &
         rate*strike*exp(-rate*time)*normal_cdf(-d2)
      out%vega = spot*root_t*exp(-dividend_yield*time)*normal_pdf(d1)
      out%status%ok = .true.
   end function bls_order1

   pure subroutine set_grid(grid, reference, upper)
      real(dp), allocatable, intent(out) :: grid(:)
      real(dp), intent(in) :: reference, upper
      integer :: j
      allocate(grid(11))
      do j = 1, 11
         grid(j) = upper*real(j-1, dp)/10.0_dp
      end do
      if (reference >= 0.0_dp .and. reference <= upper) grid(6) = reference
   end subroutine set_grid

   function option_call(spot, strike, rate, time, volatility, position, price) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      character(len=*), intent(in) :: position
      real(dp), intent(in), optional :: price
      type(payoff_table_t) :: out
      real(dp) :: premium, sign
      integer :: j

      premium = black_scholes_call(spot, strike, rate, time, volatility)
      if (present(price)) premium = price
      if (trim(position) == 'long') then
         sign = 1.0_dp
      else if (trim(position) == 'short') then
         sign = -1.0_dp
      else
         out%status%ok = .false.
         out%status%message = 'position must be long or short'
         return
      end if
      call set_grid(out%stock, spot, 2.0_dp*strike)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(1))
      do j = 1, size(out%stock)
         out%payoff(j) = sign*max(out%stock(j)-strike, 0.0_dp)
         out%profit(j) = out%payoff(j)-sign*premium*exp(rate*time)
      end do
      out%premiums = [premium]
   end function option_call

   function option_put(spot, strike, rate, time, volatility, position, price) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      character(len=*), intent(in) :: position
      real(dp), intent(in), optional :: price
      type(payoff_table_t) :: out
      real(dp) :: premium, sign
      integer :: j

      premium = black_scholes_put(spot, strike, rate, time, volatility)
      if (present(price)) premium = price
      if (trim(position) == 'long') then
         sign = 1.0_dp
      else if (trim(position) == 'short') then
         sign = -1.0_dp
      else
         out%status%ok = .false.
         out%status%message = 'position must be long or short'
         return
      end if
      call set_grid(out%stock, strike, 2.0_dp*strike)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(1))
      do j = 1, size(out%stock)
         out%payoff(j) = sign*max(strike-out%stock(j), 0.0_dp)
         out%profit(j) = out%payoff(j)-sign*premium*exp(rate*time)
      end do
      out%premiums = [premium]
   end function option_put

   function covered_call(spot, strike, rate, time, volatility, price) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      real(dp), intent(in), optional :: price
      type(payoff_table_t) :: out
      real(dp) :: premium
      integer :: j
      premium = black_scholes_call(spot, strike, rate, time, volatility)
      if (present(price)) premium = price
      call set_grid(out%stock, strike, 2.0_dp*strike)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(1))
      do j = 1, size(out%stock)
         out%payoff(j) = min(out%stock(j), strike)
         out%profit(j) = out%payoff(j)+premium*exp(rate*time)-spot
      end do
      out%premiums = [premium]
   end function covered_call

   function covered_put(spot, strike, rate, time, volatility, price) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      real(dp), intent(in), optional :: price
      type(payoff_table_t) :: out
      real(dp) :: premium
      integer :: j
      premium = black_scholes_put(spot, strike, rate, time, volatility)
      if (present(price)) premium = price
      call set_grid(out%stock, strike, 2.0_dp*strike)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(1))
      do j = 1, size(out%stock)
         if (out%stock(j) <= strike) then
            out%payoff(j) = spot-strike
         else
            out%payoff(j) = spot-out%stock(j)
         end if
         out%profit(j) = out%payoff(j)+premium*exp(rate*time)
      end do
      out%premiums = [premium]
   end function covered_put

   function protective_put(spot, strike, rate, time, volatility, price) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      real(dp), intent(in), optional :: price
      type(payoff_table_t) :: out
      real(dp) :: premium
      integer :: j
      premium = black_scholes_put(spot, strike, rate, time, volatility)
      if (present(price)) premium = price
      call set_grid(out%stock, strike, 2.0_dp*strike)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(1))
      do j = 1, size(out%stock)
         out%payoff(j) = max(out%stock(j), strike)-spot
         out%profit(j) = out%payoff(j)-premium*exp(rate*time)
      end do
      out%premiums = [premium]
   end function protective_put

   function straddle(spot, strike, rate, time, call_price, put_price, position) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, call_price, put_price
      character(len=*), intent(in) :: position
      type(payoff_table_t) :: out
      real(dp) :: sign
      integer :: j
      if (trim(position) == 'long') then
         sign = 1.0_dp
      else if (trim(position) == 'short') then
         sign = -1.0_dp
      else
         out%status%ok = .false.
         out%status%message = 'position must be long or short'
         return
      end if
      call set_grid(out%stock, spot, 2.0_dp*strike)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(3))
      do j = 1, size(out%stock)
         out%payoff(j) = sign*abs(out%stock(j)-strike)
         out%profit(j) = out%payoff(j)-sign*(call_price+put_price)*exp(rate*time)
      end do
      out%premiums = [call_price, put_price, sign*(call_price+put_price)]
   end function straddle

   function straddle_bls(spot, strike, rate, time, volatility, position) result(out)
      real(dp), intent(in) :: spot, strike, rate, time, volatility
      character(len=*), intent(in) :: position
      type(payoff_table_t) :: out
      out = straddle(spot, strike, rate, time, &
         black_scholes_call(spot, strike, rate, time, volatility), &
         black_scholes_put(spot, strike, rate, time, volatility), position)
   end function straddle_bls

   function strangle(spot, strike_put, strike_call, rate, time, put_price, call_price) result(out)
      real(dp), intent(in) :: spot, strike_put, strike_call, rate, time, put_price, call_price
      type(payoff_table_t) :: out
      integer :: j
      call set_grid(out%stock, spot, strike_call+strike_put)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(3))
      do j = 1, size(out%stock)
         out%payoff(j) = max(strike_put-out%stock(j), 0.0_dp) + max(out%stock(j)-strike_call, 0.0_dp)
         out%profit(j) = out%payoff(j)-(put_price+call_price)*exp(rate*time)
      end do
      out%premiums = [put_price, call_price, put_price+call_price]
   end function strangle

   function strangle_bls(spot, strike_put, strike_call, rate, time, volatility) result(out)
      real(dp), intent(in) :: spot, strike_put, strike_call, rate, time, volatility
      type(payoff_table_t) :: out
      out = strangle(spot, strike_put, strike_call, rate, time, &
         black_scholes_put(spot, strike_put, rate, time, volatility), &
         black_scholes_call(spot, strike_call, rate, time, volatility))
   end function strangle_bls

   function bull_call(spot, strike_low, strike_high, rate, time, low_price, high_price) result(out)
      real(dp), intent(in) :: spot, strike_low, strike_high, rate, time, low_price, high_price
      type(payoff_table_t) :: out
      integer :: j
      call set_grid(out%stock, spot, strike_high+strike_low)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(3))
      do j = 1, size(out%stock)
         out%payoff(j) = max(out%stock(j)-strike_low, 0.0_dp)-max(out%stock(j)-strike_high, 0.0_dp)
         out%profit(j) = out%payoff(j)+(high_price-low_price)*exp(rate*time)
      end do
      out%premiums = [low_price, high_price, low_price-high_price]
   end function bull_call

   function bull_call_bls(spot, strike_low, strike_high, rate, time, volatility) result(out)
      real(dp), intent(in) :: spot, strike_low, strike_high, rate, time, volatility
      type(payoff_table_t) :: out
      out = bull_call(spot, strike_low, strike_high, rate, time, &
         black_scholes_call(spot, strike_low, rate, time, volatility), &
         black_scholes_call(spot, strike_high, rate, time, volatility))
   end function bull_call_bls

   function bear_call(spot, strike_low, strike_high, rate, time, low_price, high_price) result(out)
      real(dp), intent(in) :: spot, strike_low, strike_high, rate, time, low_price, high_price
      type(payoff_table_t) :: out
      integer :: j
      call set_grid(out%stock, spot, strike_high+strike_low)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(3))
      do j = 1, size(out%stock)
         out%payoff(j) = -max(out%stock(j)-strike_low, 0.0_dp)+max(out%stock(j)-strike_high, 0.0_dp)
         out%profit(j) = out%payoff(j)+(low_price-high_price)*exp(rate*time)
      end do
      out%premiums = [low_price, high_price, high_price-low_price]
   end function bear_call

   function bear_call_bls(spot, strike_low, strike_high, rate, time, volatility) result(out)
      real(dp), intent(in) :: spot, strike_low, strike_high, rate, time, volatility
      type(payoff_table_t) :: out
      out = bear_call(spot, strike_low, strike_high, rate, time, &
         black_scholes_call(spot, strike_low, rate, time, volatility), &
         black_scholes_call(spot, strike_high, rate, time, volatility))
   end function bear_call_bls

   function collar(spot, strike_put, strike_call, rate, time, put_price, call_price) result(out)
      real(dp), intent(in) :: spot, strike_put, strike_call, rate, time, put_price, call_price
      type(payoff_table_t) :: out
      integer :: j
      call set_grid(out%stock, spot, strike_call+strike_put)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(3))
      do j = 1, size(out%stock)
         out%payoff(j) = max(strike_put-out%stock(j), 0.0_dp)-max(out%stock(j)-strike_call, 0.0_dp)
         out%profit(j) = out%payoff(j)+(call_price-put_price)*exp(rate*time)
      end do
      out%premiums = [put_price, call_price, put_price-call_price]
   end function collar

   function collar_bls(spot, strike_put, strike_call, rate, time, volatility) result(out)
      real(dp), intent(in) :: spot, strike_put, strike_call, rate, time, volatility
      type(payoff_table_t) :: out
      out = collar(spot, strike_put, strike_call, rate, time, &
         black_scholes_put(spot, strike_put, rate, time, volatility), &
         black_scholes_call(spot, strike_call, rate, time, volatility))
   end function collar_bls

   function butterfly_spread(spot, strike_low, strike_mid, strike_high, rate, time, &
      low_price, mid_price, high_price) result(out)
      real(dp), intent(in) :: spot, strike_low, strike_mid, strike_high, rate, time
      real(dp), intent(in) :: low_price, mid_price, high_price
      type(payoff_table_t) :: out
      integer :: j
      call set_grid(out%stock, spot, strike_high+strike_low)
      allocate(out%payoff(size(out%stock)), out%profit(size(out%stock)), out%premiums(6))
      do j = 1, size(out%stock)
         out%payoff(j) = max(out%stock(j)-strike_low, 0.0_dp) - &
            2.0_dp*max(out%stock(j)-strike_mid, 0.0_dp) + max(out%stock(j)-strike_high, 0.0_dp)
         out%profit(j) = out%payoff(j)+(2.0_dp*mid_price-low_price-high_price)*exp(rate*time)
      end do
      out%premiums = [low_price, mid_price, high_price, low_price+high_price, &
         2.0_dp*mid_price, low_price+high_price-2.0_dp*mid_price]
   end function butterfly_spread

   function butterfly_spread_bls(spot, strike_low, strike_mid, strike_high, rate, time, volatility) result(out)
      real(dp), intent(in) :: spot, strike_low, strike_mid, strike_high, rate, time, volatility
      type(payoff_table_t) :: out
      out = butterfly_spread(spot, strike_low, strike_mid, strike_high, rate, time, &
         black_scholes_call(spot, strike_low, rate, time, volatility), &
         black_scholes_call(spot, strike_mid, rate, time, volatility), &
         black_scholes_call(spot, strike_high, rate, time, volatility))
   end function butterfly_spread_bls

   function forward_price(spot, time, force_rate, dividend_structure, dividend, dividend_frequency, &
      dividend_yield, growth_rate) result(price)
      real(dp), intent(in) :: spot, time, force_rate, dividend, dividend_frequency, dividend_yield, growth_rate
      character(len=*), intent(in) :: dividend_structure
      real(dp) :: price, j, n
      select case (trim(dividend_structure))
      case ('none')
         price = spot*exp(force_rate*time)
      case ('continuous')
         price = spot*exp((force_rate-dividend_yield)*time)
      case ('discrete')
         j = exp(force_rate/dividend_frequency)-1.0_dp
         n = time*dividend_frequency
         if (abs(growth_rate-j) < 1.0e-12_dp) then
            price = spot*exp(force_rate*time)-dividend*n/(1.0_dp+j)*exp(force_rate*time)
         else if (growth_rate > -0.999999_dp) then
            price = spot*exp(force_rate*time)-dividend*(1.0_dp-((1.0_dp+growth_rate)/(1.0_dp+j))**n)/ &
               (j-growth_rate)*exp(force_rate*time)
         else
            price = spot*exp(force_rate*time)-dividend*((1.0_dp+j)**n-1.0_dp)/j
         end if
      case default
         price = huge(1.0_dp)
      end select
   end function forward_price

   function prepaid_forward_price(spot, time, force_rate, dividend_structure, dividend, dividend_frequency, &
      dividend_yield, growth_rate) result(price)
      real(dp), intent(in) :: spot, time, force_rate, dividend, dividend_frequency, dividend_yield, growth_rate
      character(len=*), intent(in) :: dividend_structure
      real(dp) :: price
      price = forward_price(spot, time, force_rate, dividend_structure, dividend, dividend_frequency, &
         dividend_yield, growth_rate)*exp(-force_rate*time)
   end function prepaid_forward_price

   function forward_contract(spot, time, force_rate, position, dividend_structure, dividend, &
      dividend_frequency, dividend_yield, growth_rate) result(out)
      real(dp), intent(in) :: spot, time, force_rate, dividend, dividend_frequency, dividend_yield, growth_rate
      character(len=*), intent(in) :: position, dividend_structure
      type(forward_result_t) :: out
      real(dp) :: sign
      integer :: j
      if (trim(position) == 'long') then
         sign = 1.0_dp
      else if (trim(position) == 'short') then
         sign = -1.0_dp
      else
         out%status%ok = .false.
         out%status%message = 'position must be long or short'
         return
      end if
      out%delivery_price = forward_price(spot, time, force_rate, dividend_structure, dividend, &
         dividend_frequency, dividend_yield, growth_rate)
      out%prepaid_price = out%delivery_price*exp(-force_rate*time)
      call set_grid(out%stock, out%delivery_price, 2.0_dp*out%delivery_price)
      allocate(out%payoff(size(out%stock)))
      do j = 1, size(out%stock)
         out%payoff(j) = sign*(out%stock(j)-out%delivery_price)
      end do
   end function forward_contract

   function forward_prepaid(spot, time, force_rate, position, dividend_structure, dividend, &
      dividend_frequency, dividend_yield, growth_rate) result(out)
      real(dp), intent(in) :: spot, time, force_rate, dividend, dividend_frequency, dividend_yield, growth_rate
      character(len=*), intent(in) :: position, dividend_structure
      type(forward_result_t) :: out
      real(dp) :: sign
      integer :: j
      if (trim(position) == 'long') then
         sign = 1.0_dp
      else if (trim(position) == 'short') then
         sign = -1.0_dp
      else
         out%status%ok = .false.
         out%status%message = 'position must be long or short'
         return
      end if
      out%delivery_price = forward_price(spot, time, force_rate, dividend_structure, dividend, &
         dividend_frequency, dividend_yield, growth_rate)
      out%prepaid_price = out%delivery_price*exp(-force_rate*time)
      call set_grid(out%stock, out%prepaid_price, 2.0_dp*out%prepaid_price)
      allocate(out%payoff(size(out%stock)))
      do j = 1, size(out%stock)
         out%payoff(j) = sign*(out%stock(j)-out%prepaid_price)
      end do
   end function forward_prepaid

end module financialmath_derivatives
