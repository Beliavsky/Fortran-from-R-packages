! SPDX-License-Identifier: GPL-2.0-or-later
module fer_spread
   use fer_kinds, only : dp
   use fer_special, only : normal_cdf
   use fer_vanilla, only : black_scholes_price, bachelier_price
   implicit none
   private
   public :: switch_margrabe, spread_kirk, spread_bjerksund_2014, spread_bachelier
contains
   elemental real(dp) function switch_margrabe(forward1,forward2,texp,sigma1,sigma2,corr,cp,df) result(value)
      real(dp), intent(in) :: forward1,forward2,texp,sigma1,sigma2,corr,df
      integer, intent(in) :: cp
      real(dp) :: vol
      vol = sqrt(max(0.0_dp,sigma1*sigma1-2.0_dp*corr*sigma1*sigma2+sigma2*sigma2))
      value = black_scholes_price(forward2,forward1,texp,vol,cp,df)
   end function switch_margrabe

   elemental real(dp) function spread_kirk(strike,forward1,forward2,texp,sigma1,sigma2,corr,cp,df) result(value)
      real(dp), intent(in) :: strike,forward1,forward2,texp,sigma1,sigma2,corr,df
      integer, intent(in) :: cp
      real(dp) :: kp, km, s1, s2
      kp = max(strike,0.0_dp)
      km = min(strike,0.0_dp)
      s1 = sigma1*forward1/(forward1-km)
      s2 = sigma2*forward2/(forward2+kp)
      value = switch_margrabe(forward1-km,forward2+kp,texp,s1,s2,corr,cp,df)
   end function spread_kirk

   elemental real(dp) function spread_bjerksund_2014(strike,forward1,forward2,texp,sigma1,sigma2,corr,cp,df) result(value)
      real(dp), intent(in) :: strike,forward1,forward2,texp,sigma1,sigma2,corr,df
      integer, intent(in) :: cp
      real(dp) :: std11,std12,std22,a,b,std,d0,d1,d2,d3
      std11=sigma1*sigma1*texp; std12=sigma1*sigma2*texp; std22=sigma2*sigma2*texp
      a=forward2+strike; b=forward2/a
      std=sqrt(std11-2.0_dp*b*corr*std12+b*b*std22)
      d0=log(forward1/a)
      d1=(d0+0.5_dp*std11-b*(corr*std12-0.5_dp*b*std22))/std
      d2=(d0-0.5_dp*std11+corr*std12+b*(0.5_dp*b-1.0_dp)*std22)/std
      d3=(d0-0.5_dp*std11+0.5_dp*b*b*std22)/std
      value=df*real(cp,dp)*(forward1*normal_cdf(real(cp,dp)*d1)- &
         forward2*normal_cdf(real(cp,dp)*d2)-strike*normal_cdf(real(cp,dp)*d3))
   end function spread_bjerksund_2014

   elemental real(dp) function spread_bachelier(strike,forward1,forward2,texp,sigma1,sigma2,corr,cp,df) result(value)
      real(dp), intent(in) :: strike,forward1,forward2,texp,sigma1,sigma2,corr,df
      integer, intent(in) :: cp
      real(dp) :: sigma
      sigma=sqrt(max(0.0_dp,sigma1*sigma1-2.0_dp*corr*sigma1*sigma2+sigma2*sigma2))
      value=bachelier_price(strike,forward1-forward2,texp,sigma,cp,df)
   end function spread_bachelier
end module fer_spread
