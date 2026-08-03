! SPDX-License-Identifier: GPL-2.0-or-later
module fer_vanilla
   use fer_kinds, only : dp
   use fer_special, only : normal_pdf, normal_cdf
   implicit none
   private
   public :: bachelier_price, bachelier_impvol
   public :: black_scholes_price, black_scholes_impvol
contains
   elemental real(dp) function bachelier_price(strike, forward, texp, sigma, cp, df) result(price)
      real(dp), intent(in) :: strike, forward, texp, sigma
      integer, intent(in) :: cp
      real(dp), intent(in), optional :: df
      real(dp) :: disc, stdev, dn
      disc = 1.0_dp
      if (present(df)) disc = df
      stdev = max(abs(sigma)*sqrt(max(texp,0.0_dp)), 1.0e-32_dp)
      dn = real(cp,dp)*(forward-strike)/stdev
      price = disc*(real(cp,dp)*(forward-strike)*normal_cdf(dn) + stdev*normal_pdf(dn))
   end function bachelier_price

   real(dp) function bachelier_impvol(price, strike, forward, texp, cp, df) result(vol)
      real(dp), intent(in) :: price, strike, forward, texp
      integer, intent(in) :: cp
      real(dp), intent(in), optional :: df
      real(dp) :: disc, intrinsic, lo, hi, mid, pmid
      integer :: iter
      disc = 1.0_dp
      if (present(df)) disc = df
      intrinsic = disc*max(real(cp,dp)*(forward-strike),0.0_dp)
      if (price <= intrinsic + 32.0_dp*epsilon(price)) then
         vol = 0.0_dp
         return
      end if
      lo = 0.0_dp
      hi = max(1.0_dp, abs(forward-strike)/sqrt(max(texp,epsilon(1.0_dp))))
      do while (bachelier_price(strike,forward,texp,hi,cp,disc) < price)
         hi = 2.0_dp*hi
         if (hi > huge(1.0_dp)/4.0_dp) exit
      end do
      do iter = 1, 120
         mid = 0.5_dp*(lo+hi)
         pmid = bachelier_price(strike,forward,texp,mid,cp,disc)
         if (pmid < price) then
            lo = mid
         else
            hi = mid
         end if
      end do
      vol = 0.5_dp*(lo+hi)
   end function bachelier_impvol

   elemental real(dp) function black_scholes_price(strike, forward, texp, sigma, cp, df) result(price)
      real(dp), intent(in) :: strike, forward, texp, sigma
      integer, intent(in) :: cp
      real(dp), intent(in), optional :: df
      real(dp) :: disc, stdev, d1, d2
      disc = 1.0_dp
      if (present(df)) disc = df
      if (strike <= 0.0_dp .or. forward <= 0.0_dp) then
         price = disc*max(real(cp,dp)*(forward-strike),0.0_dp)
         return
      end if
      stdev = max(abs(sigma)*sqrt(max(texp,0.0_dp)), epsilon(1.0_dp))
      d1 = log(forward/strike)/stdev + 0.5_dp*stdev
      d2 = d1-stdev
      price = disc*real(cp,dp)*(forward*normal_cdf(real(cp,dp)*d1)-strike*normal_cdf(real(cp,dp)*d2))
   end function black_scholes_price

   real(dp) function black_scholes_impvol(price, strike, forward, texp, cp, df) result(vol)
      real(dp), intent(in) :: price, strike, forward, texp
      integer, intent(in) :: cp
      real(dp), intent(in), optional :: df
      real(dp) :: disc, intrinsic, maximum, lo, hi, mid, pmid
      integer :: iter
      disc = 1.0_dp
      if (present(df)) disc = df
      intrinsic = disc*max(real(cp,dp)*(forward-strike),0.0_dp)
      maximum = disc*merge(forward,strike,cp>0)
      if (price <= intrinsic + 32.0_dp*epsilon(price)) then
         vol = 0.0_dp
         return
      end if
      if (price >= maximum) then
         vol = huge(1.0_dp)
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      do while (black_scholes_price(strike,forward,texp,hi,cp,disc) < price)
         hi = 2.0_dp*hi
         if (hi >= 1.0e4_dp) exit
      end do
      do iter = 1, 120
         mid = 0.5_dp*(lo+hi)
         pmid = black_scholes_price(strike,forward,texp,mid,cp,disc)
         if (pmid < price) then
            lo = mid
         else
            hi = mid
         end if
      end do
      vol = 0.5_dp*(lo+hi)
   end function black_scholes_impvol
end module fer_vanilla
