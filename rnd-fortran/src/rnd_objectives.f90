! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_objectives
   use rnd_kinds, only : dp
   use rnd_types, only : option_prices, am_price_result
   use rnd_special, only : normal_cdf, beta_function, is_finite_real
   use rnd_pricing, only : price_gb_option, price_mln_option, price_ew_option, price_am_option
   implicit none
   private
   public :: bsm_objective, gb_objective, mln_objective, ew_objective, mln_am_objective

contains

   real(dp) function bsm_objective(theta, s0, r, te, dividend_yield, market_calls, &
         call_strikes, call_weights, market_puts, put_strikes, put_weights, lambda) result(value)
      real(dp), intent(in) :: theta(:), s0, r, te, dividend_yield
      real(dp), intent(in) :: market_calls(:), call_strikes(:), call_weights(:)
      real(dp), intent(in) :: market_puts(:), put_strikes(:), put_weights(:), lambda
      real(dp) :: mu, zeta, discount, expected_value
      real(dp), allocatable :: d1(:), d2(:), c1(:), c2(:), theoretical_calls(:), theoretical_puts(:)
      mu = theta(1)
      zeta = theta(2)
      if (zeta <= 0.0_dp) then
         value = 1.0e7_dp
         return
      end if
      discount = exp(-r*te)
      expected_value = exp(mu+0.5_dp*zeta*zeta)
      allocate(d1(size(call_strikes)),d2(size(call_strikes)),theoretical_calls(size(call_strikes)))
      d1 = (log(call_strikes)-mu-zeta*zeta)/zeta
      d2 = (log(call_strikes)-mu)/zeta
      theoretical_calls = discount*(expected_value*(1.0_dp-normal_cdf(d1)) &
         - call_strikes*(1.0_dp-normal_cdf(d2)))
      allocate(c1(size(put_strikes)),c2(size(put_strikes)),theoretical_puts(size(put_strikes)))
      c1 = (log(put_strikes)-mu-zeta*zeta)/zeta
      c2 = (log(put_strikes)-mu)/zeta
      theoretical_puts = discount*(put_strikes*normal_cdf(c2)-normal_cdf(c1)*expected_value)
      value = sum(call_weights*(theoretical_calls-market_calls)**2) &
         + sum(put_weights*(theoretical_puts-market_puts)**2) &
         + lambda*(s0*exp(-dividend_yield*te)-expected_value*discount)**2
   end function bsm_objective

   real(dp) function gb_objective(theta, r, te, dividend_yield, s0, market_calls, &
         call_strikes, call_weights, market_puts, put_strikes, put_weights, lambda) result(value)
      real(dp), intent(in) :: theta(:), r, te, dividend_yield, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:), call_weights(:)
      real(dp), intent(in) :: market_puts(:), put_strikes(:), put_weights(:), lambda
      real(dp) :: a, b, v, w, expected_value, discount
      type(option_prices) :: calls, puts
      a = theta(1)
      b = theta(2)
      v = theta(3)
      w = theta(4)
      if (a <= 0.0_dp .or. a > 20.0_dp .or. b <= 0.0_dp .or. v <= 0.0_dp &
            .or. w <= 0.0_dp .or. w <= 4.0_dp/a) then
         value = 1.0e7_dp
         return
      end if
      expected_value = b*beta_function(v+1.0_dp/a,w-1.0_dp/a)/beta_function(v,w)
      if (.not. is_finite_real(expected_value)) then
         value = 1.0e7_dp
         return
      end if
      discount = exp(-r*te)
      calls = price_gb_option(r,te,s0,call_strikes,dividend_yield,a,b,v,w)
      puts = price_gb_option(r,te,s0,put_strikes,dividend_yield,a,b,v,w)
      value = sum(call_weights*(calls%call-market_calls)**2) &
         + sum(put_weights*(puts%put-market_puts)**2) &
         + lambda*(s0*exp(-dividend_yield*te)-expected_value*discount)**2
   end function gb_objective

   real(dp) function mln_objective(theta, r, dividend_yield, te, s0, market_calls, &
         call_strikes, call_weights, market_puts, put_strikes, put_weights, lambda) result(value)
      real(dp), intent(in) :: theta(:), r, dividend_yield, te, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:), call_weights(:)
      real(dp), intent(in) :: market_puts(:), put_strikes(:), put_weights(:), lambda
      real(dp) :: alpha1, meanlog1, meanlog2, sdlog1, sdlog2, expected_value
      type(option_prices) :: calls, puts
      alpha1 = theta(1)
      meanlog1 = theta(2)
      meanlog2 = theta(3)
      sdlog1 = theta(4)
      sdlog2 = theta(5)
      if (alpha1 < 0.0_dp .or. alpha1 > 1.0_dp .or. sdlog1 <= 0.0_dp .or. sdlog2 <= 0.0_dp) then
         value = 1.0e7_dp
         return
      end if
      expected_value = alpha1*exp(meanlog1+0.5_dp*sdlog1*sdlog1) &
         +(1.0_dp-alpha1)*exp(meanlog2+0.5_dp*sdlog2*sdlog2)
      calls = price_mln_option(r,te,dividend_yield,call_strikes,alpha1,meanlog1,meanlog2,sdlog1,sdlog2)
      puts = price_mln_option(r,te,dividend_yield,put_strikes,alpha1,meanlog1,meanlog2,sdlog1,sdlog2)
      value = sum(call_weights*(calls%call-market_calls)**2) &
         + sum(put_weights*(puts%put-market_puts)**2) &
         + lambda*(s0-exp(dividend_yield*te)*expected_value*exp(-r*te))**2
   end function mln_objective

   real(dp) function ew_objective(theta, r, dividend_yield, te, s0, market_calls, &
         call_strikes, call_weights, lambda) result(value)
      real(dp), intent(in) :: theta(:), r, dividend_yield, te, s0
      real(dp), intent(in) :: market_calls(:), call_strikes(:), call_weights(:), lambda
      real(dp) :: sigma, skew, kurt, expected_value, discount, m
      type(option_prices) :: calls
      logical, allocatable :: valid(:)
      sigma = theta(1)
      skew = theta(2)
      kurt = theta(3)
      if (sigma <= 0.0_dp .or. kurt < 0.0_dp) then
         value = 1.0e7_dp
         return
      end if
      discount = exp(-r*te)
      m = log(s0)+(r-dividend_yield-0.5_dp*sigma*sigma)*te
      expected_value = exp(m+0.5_dp*sigma*sigma*te)
      calls = price_ew_option(r,te,s0,call_strikes,sigma,dividend_yield,skew,kurt)
      allocate(valid(size(call_strikes)))
      valid = is_finite_real(calls%call)
      if (.not. any(valid)) then
         value = 1.0e7_dp
      else
         value = sum(pack(call_weights*(calls%call-market_calls)**2,valid)) &
            + lambda*(s0*exp(-dividend_yield*te)-expected_value*discount)**2
      end if
   end function ew_objective

   real(dp) function mln_am_objective(theta, s0, r, te, market_calls, call_weights, &
         market_puts, put_weights, strikes, lambda) result(value)
      real(dp), intent(in) :: theta(:), s0, r, te
      real(dp), intent(in) :: market_calls(:), call_weights(:), market_puts(:), put_weights(:)
      real(dp), intent(in) :: strikes(:), lambda
      real(dp) :: w1, w2, u1, u2, u3, sigma1, sigma2, sigma3, p1, p2, expected_f0
      real(dp) :: theoretical_call, theoretical_put
      real(dp) :: call_sum, put_sum, call_weight, put_weight
      type(am_price_result) :: price
      integer :: i
      w1 = theta(1)
      w2 = theta(2)
      u1 = theta(3)
      u2 = theta(4)
      u3 = theta(5)
      sigma1 = theta(6)
      sigma2 = theta(7)
      sigma3 = theta(8)
      p1 = theta(9)
      p2 = theta(10)
      if (w1 < 0.0_dp .or. w1 > 1.0_dp .or. w2 < 0.0_dp .or. w2 > 1.0_dp &
            .or. min(sigma1,sigma2,sigma3) <= 0.0_dp .or. p1 < 0.0_dp .or. p2 < 0.0_dp &
            .or. p1+p2 > 1.0_dp) then
         value = 1.0e7_dp
         return
      end if
      expected_f0 = p1*exp(u1+0.5_dp*sigma1*sigma1) &
         + p2*exp(u2+0.5_dp*sigma2*sigma2) &
         + (1.0_dp-p1-p2)*exp(u3+0.5_dp*sigma3*sigma3)
      call_sum = 0.0_dp
      put_sum = 0.0_dp
      do i = 1, size(strikes)
         if (strikes(i) < expected_f0) then
            call_weight = w1
            put_weight = w2
         else
            call_weight = w2
            put_weight = w1
         end if
         price = price_am_option(strikes(i),r,te,call_weight,u1,u2,u3,sigma1,sigma2,sigma3,p1,p2)
         theoretical_call = price%call
         price = price_am_option(strikes(i),r,te,put_weight,u1,u2,u3,sigma1,sigma2,sigma3,p1,p2)
         theoretical_put = price%put
         call_sum = call_sum+call_weights(i)*(theoretical_call-market_calls(i))**2
         put_sum = put_sum+put_weights(i)*(theoretical_put-market_puts(i))**2
      end do
      value = call_sum+put_sum+lambda*(s0-expected_f0*exp(-r*te))**2
   end function mln_am_objective

end module rnd_objectives
