! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_pricing
   use rnd_kinds, only : dp
   use rnd_types, only : option_prices, am_price_result
   use rnd_special, only : normal_cdf
   use rnd_densities, only : pgb, approximate_max
   implicit none
   private
   public :: price_bsm_option, price_gb_option, price_mln_option
   public :: price_ew_option, price_shimko_option, price_am_option

contains

   function price_bsm_option(s0, strike, r, te, sigma, dividend_yield) result(out)
      real(dp), intent(in) :: s0, strike(:), r, te, sigma, dividend_yield
      type(option_prices) :: out
      real(dp), allocatable :: d1(:), d2(:)
      integer :: n
      n = size(strike)
      allocate(out%call(n), out%put(n), d1(n), d2(n))
      if (sigma <= 0.0_dp .or. te <= 0.0_dp) then
         out%call = max(0.0_dp, s0*exp(-dividend_yield*te)-strike*exp(-r*te))
         out%put = max(0.0_dp, strike*exp(-r*te)-s0*exp(-dividend_yield*te))
         return
      end if
      d1 = (log(s0/strike)+(r-dividend_yield+0.5_dp*sigma*sigma)*te)/(sigma*sqrt(te))
      d2 = d1-sigma*sqrt(te)
      out%call = s0*exp(-dividend_yield*te)*normal_cdf(d1) &
         - strike*exp(-r*te)*normal_cdf(d2)
      out%put = strike*exp(-r*te)*normal_cdf(-d2) &
         - s0*exp(-dividend_yield*te)*normal_cdf(-d1)
   end function price_bsm_option

   function price_gb_option(r, te, s0, strike, dividend_yield, a, b, v, w) result(out)
      real(dp), intent(in) :: r, te, s0, strike(:), dividend_yield, a, b, v, w
      type(option_prices) :: out
      real(dp), allocatable :: prob1(:), prob2(:)
      integer :: n
      n = size(strike)
      allocate(out%call(n), out%put(n), prob1(n), prob2(n))
      prob1 = pgb(strike, a, b, v+1.0_dp/a, w-1.0_dp/a)
      prob2 = pgb(strike, a, b, v, w)
      out%call = s0*exp(-dividend_yield*te)*(1.0_dp-prob1) &
         - strike*exp(-r*te)*(1.0_dp-prob2)
      out%put = out%call-s0*exp(-dividend_yield*te)+strike*exp(-r*te)
   end function price_gb_option

   function price_mln_option(r, te, dividend_yield, strike, alpha1, meanlog1, meanlog2, &
         sdlog1, sdlog2) result(out)
      real(dp), intent(in) :: r, te, dividend_yield, strike(:)
      real(dp), intent(in) :: alpha1, meanlog1, meanlog2, sdlog1, sdlog2
      type(option_prices) :: out
      real(dp), allocatable :: u1(:), u2(:), c1(:), c2(:)
      real(dp) :: alpha2, expected1, expected2, discount, s0
      integer :: n
      n = size(strike)
      allocate(out%call(n), out%put(n), u1(n), u2(n), c1(n), c2(n))
      alpha2 = 1.0_dp-alpha1
      expected1 = exp(meanlog1+0.5_dp*sdlog1*sdlog1)
      expected2 = exp(meanlog2+0.5_dp*sdlog2*sdlog2)
      discount = exp(-r*te)
      s0 = exp((dividend_yield-r)*te)*(alpha1*expected1+alpha2*expected2)
      u1 = (log(strike)-meanlog1)/sdlog1
      c1 = discount*(expected1*(1.0_dp-normal_cdf(u1-sdlog1)) &
         - strike*(1.0_dp-normal_cdf(u1)))
      u2 = (log(strike)-meanlog2)/sdlog2
      c2 = discount*(expected2*(1.0_dp-normal_cdf(u2-sdlog2)) &
         - strike*(1.0_dp-normal_cdf(u2)))
      out%call = alpha1*c1+alpha2*c2
      out%put = out%call-s0*exp(-dividend_yield*te)+strike*discount
   end function price_mln_option

   function price_ew_option(r, te, s0, strike, sigma, dividend_yield, skew, kurt) result(out)
      real(dp), intent(in) :: r, te, s0, strike(:), sigma, dividend_yield, skew, kurt
      type(option_prices) :: out
      type(option_prices) :: bsm
      real(dp), allocatable :: density(:), first_derivative(:), second_derivative(:)
      real(dp) :: discount, v, m, skew_lognorm, kurt_lognorm, cumul_lognorm
      integer :: n
      n = size(strike)
      allocate(out%call(n), out%put(n), density(n), first_derivative(n), second_derivative(n))
      bsm = price_bsm_option(s0, strike, r, te, sigma, dividend_yield)
      discount = exp(-r*te)
      v = sqrt(exp(sigma*sigma*te)-1.0_dp)
      m = log(s0)+(r-dividend_yield-0.5_dp*sigma*sigma)*te
      skew_lognorm = 3.0_dp*v+v**3
      kurt_lognorm = 16.0_dp*v**2+15.0_dp*v**4+6.0_dp*v**6+v**8
      cumul_lognorm = (s0*exp((r-dividend_yield)*te)*v)**2
      density = exp(-0.5_dp*((log(strike)-m)/(sigma*sqrt(te)))**2) &
         /(strike*sigma*sqrt(te)*sqrt(2.0_dp*acos(-1.0_dp)))
      first_derivative = -(1.0_dp+(log(strike)-m)/(te*sigma*sigma))*density/strike
      second_derivative = -(2.0_dp+(log(strike)-m)/(te*sigma*sigma))*first_derivative/strike &
         - density/(strike*strike*sigma*sigma)
      out%call = bsm%call-discount*(skew-skew_lognorm)*cumul_lognorm**1.5_dp*first_derivative/6.0_dp &
         + discount*(kurt-kurt_lognorm)*cumul_lognorm**2*second_derivative/24.0_dp
      out%put = out%call+strike*exp(-r*te)-s0*exp(-dividend_yield*te)
   end function price_ew_option

   function price_shimko_option(r, te, s0, strike, dividend_yield, a0, a1, a2) result(out)
      real(dp), intent(in) :: r, te, s0, strike(:), dividend_yield, a0, a1, a2
      type(option_prices) :: out
      type(option_prices) :: one
      real(dp) :: local_sigma
      integer :: i, n
      n = size(strike)
      allocate(out%call(n), out%put(n))
      do i = 1, n
         local_sigma = a0+a1*strike(i)+a2*strike(i)*strike(i)
         one = price_bsm_option(s0, strike(i:i), r, te, local_sigma, dividend_yield)
         out%call(i) = one%call(1)
         out%put(i) = one%put(1)
      end do
   end function price_shimko_option

   function price_am_option(strike, r, te, weight, u1, u2, u3, sigma1, sigma2, sigma3, &
         p1, p2) result(out)
      real(dp), intent(in) :: strike, r, te, weight, u1, u2, u3
      real(dp), intent(in) :: sigma1, sigma2, sigma3, p1, p2
      type(am_price_result) :: out
      real(dp) :: mu(3), sigma(3), proportion(3), above_terms(3)
      real(dp), parameter :: eps = 100.0_dp*epsilon(1.0_dp)
      mu = [u1, u2, u3]
      sigma = [sigma1, sigma2, sigma3]
      proportion = [p1, p2, 1.0_dp-p1-p2]
      out%expected_f0 = sum(proportion*exp(mu+0.5_dp*sigma*sigma))
      above_terms = 1.0_dp-normal_cdf((log(strike)-mu)/sigma)
      out%prob_above = sum(proportion*above_terms)
      out%prob_below = 1.0_dp-out%prob_above
      if (out%prob_above > eps) then
         out%conditional_above = sum(proportion*exp(mu+0.5_dp*sigma*sigma) &
            *(1.0_dp-normal_cdf((log(strike)-mu-sigma*sigma)/sigma)))/out%prob_above
      else
         out%conditional_above = strike
      end if
      if (out%prob_below > eps) then
         out%conditional_below = (out%expected_f0-out%conditional_above*out%prob_above) &
            /out%prob_below
      else
         out%conditional_below = strike
      end if
      out%call = weight*(out%conditional_above-strike)*out%prob_above &
         +(1.0_dp-weight)*approximate_max(out%expected_f0-strike, &
         exp(-r*te)*(out%conditional_above-strike)*out%prob_above)
      out%call = max(0.0_dp, out%call)
      out%put = weight*(strike-out%conditional_below)*out%prob_below &
         +(1.0_dp-weight)*approximate_max(strike-out%expected_f0, &
         exp(-r*te)*(strike-out%conditional_below)*out%prob_below)
      out%put = max(0.0_dp, out%put)
   end function price_am_option

end module rnd_pricing
