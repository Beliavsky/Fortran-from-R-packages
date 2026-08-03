! SPDX-License-Identifier: GPL-2.0-or-later
module ragtop_black_scholes
   use ragtop_kinds, only : dp
   use ragtop_constants, only : call_option, put_option, ragtop_ok, ragtop_invalid_argument
   use ragtop_types, only : option_value, dividend_schedule, market_spec
   use ragtop_math, only : normal_cdf, normal_pdf
   use ragtop_term_structures, only : effective_rate, effective_volatility, effective_hazard
   use ragtop_cashflows, only : time_adjusted_dividends
   implicit none
   private
   public :: black_scholes, black_scholes_on_term_structures

contains

   function black_scholes(callput, spot, strike, rate, maturity, volatility, default_intensity, &
                          dividend_rate, borrow_cost, dividends) result(ans)
      integer, intent(in) :: callput
      real(dp), intent(in) :: spot, strike, rate, maturity, volatility
      real(dp), intent(in), optional :: default_intensity, dividend_rate, borrow_cost
      type(dividend_schedule), intent(in), optional :: dividends
      type(option_value) :: ans
      real(dp) :: h, qdiv, borrow, adjusted_spot, sd, q, d1, d2, v, delta, vega
      real(dp) :: surv, default_value
      real(dp), allocatable :: div_amount(:)
      h = 0.0_dp
      qdiv = 0.0_dp
      borrow = 0.0_dp
      if (present(default_intensity)) h = default_intensity
      if (present(dividend_rate)) qdiv = dividend_rate
      if (present(borrow_cost)) borrow = borrow_cost
      ans%status = ragtop_ok
      if ((callput /= call_option .and. callput /= put_option) .or. spot < 0.0_dp .or. &
          strike < 0.0_dp .or. maturity < 0.0_dp .or. volatility < 0.0_dp) then
         ans%status = ragtop_invalid_argument
         return
      end if
      adjusted_spot = spot
      if (present(dividends)) then
         if (allocated(dividends%time)) then
            allocate(div_amount(1))
            call time_adjusted_dividends(dividends,0.0_dp,rate,0.0_dp,[spot],spot,div_amount)
            adjusted_spot = max(0.0_dp,spot-div_amount(1))
         end if
      end if
      if (maturity <= epsilon(1.0_dp)) then
         ans%price = max(0.0_dp,real(callput,dp)*(adjusted_spot-strike))
         if (real(callput,dp)*(adjusted_spot-strike) > 0.0_dp) ans%delta = real(callput,dp)
         ans%vega = 0.0_dp
         return
      end if
      if (strike <= epsilon(1.0_dp)) then
         ans%price = adjusted_spot
         ans%delta = 1.0_dp
         ans%vega = 0.0_dp
         return
      end if
      sd = volatility*sqrt(maturity)
      q = qdiv+borrow-h
      surv = exp(-h*maturity)
      default_value = max(0.0_dp,-real(callput,dp)*strike*exp(-rate*maturity))
      if (sd <= sqrt(epsilon(1.0_dp))) then
         v = max(0.0_dp,real(callput,dp)*(adjusted_spot*exp(-q*maturity)-strike*exp(-rate*maturity)))
         ans%price = surv*v+(1.0_dp-surv)*default_value
         if (v > 0.0_dp) ans%delta = surv*real(callput,dp)*exp(-q*maturity)
         ans%vega = 0.0_dp
         return
      end if
      d1 = (log(adjusted_spot/strike)+(rate-q)*maturity+0.5_dp*sd*sd)/sd
      d2 = d1-sd
      v = real(callput,dp)*(adjusted_spot*exp(-q*maturity)*normal_cdf(real(callput,dp)*d1)- &
                            strike*exp(-rate*maturity)*normal_cdf(real(callput,dp)*d2))
      delta = exp(-q*maturity)*real(callput,dp)*normal_cdf(real(callput,dp)*d1)
      vega = adjusted_spot*exp(-q*maturity)*normal_pdf(d1)*sqrt(maturity)
      ans%price = surv*v+(1.0_dp-surv)*default_value
      ans%delta = surv*delta
      ans%vega = surv*vega
   end function black_scholes

   function black_scholes_on_term_structures(callput, spot, strike, maturity, market, dividends) result(ans)
      integer, intent(in) :: callput
      real(dp), intent(in) :: spot, strike, maturity
      type(market_spec), intent(in) :: market
      type(dividend_schedule), intent(in), optional :: dividends
      type(option_value) :: ans
      real(dp) :: r, v, h
      r = effective_rate(market,maturity)
      v = effective_volatility(market,maturity)
      h = effective_hazard(market,spot)
      if (present(dividends)) then
         ans = black_scholes(callput,spot,strike,r,maturity,v,h,market%dividend_rate,market%borrow_cost,dividends)
      else
         ans = black_scholes(callput,spot,strike,r,maturity,v,h,market%dividend_rate,market%borrow_cost)
      end if
   end function black_scholes_on_term_structures

end module ragtop_black_scholes
