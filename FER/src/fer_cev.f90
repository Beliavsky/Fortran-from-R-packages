! SPDX-License-Identifier: GPL-2.0-or-later
module fer_cev
   use fer_kinds, only : dp
   use fer_special, only : noncentral_chisq_cdf, regularized_gamma_q
   implicit none
   private
   public :: cev_price, cev_mass_zero
contains
   real(dp) function cev_price(strike, forward, texp, sigma, beta, cp, df) result(price)
      real(dp), intent(in) :: strike, forward, texp, sigma, beta
      integer, intent(in) :: cp
      real(dp), intent(in), optional :: df
      real(dp) :: disc, betac, scale, strike_cov, forward_cov, deg, term1, term2
      disc = 1.0_dp
      if (present(df)) disc = df
      betac = 1.0_dp-beta
      scale = (betac*sigma)**2*texp
      strike_cov = strike**(2.0_dp*betac)/scale
      forward_cov = forward**(2.0_dp*betac)/scale
      deg = 1.0_dp/betac
      term1 = noncentral_chisq_cdf(strike_cov,deg+2.0_dp,forward_cov,cp<0)
      term2 = noncentral_chisq_cdf(forward_cov,deg,strike_cov,cp>0)
      price = real(cp,dp)*disc*(forward*term1-strike*term2)
   end function cev_price

   real(dp) function cev_mass_zero(forward, texp, sigma, beta) result(mass)
      real(dp), intent(in) :: forward, texp, sigma, beta
      real(dp) :: betac, scale, x, shape
      betac = 1.0_dp-beta
      scale = (betac*sigma)**2*texp
      x = 0.5_dp*forward**(2.0_dp*betac)/scale
      shape = 0.5_dp/betac
      mass = regularized_gamma_q(shape,x)
   end function cev_mass_zero
end module fer_cev
