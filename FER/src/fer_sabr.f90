! SPDX-License-Identifier: GPL-2.0-or-later
module fer_sabr
   use fer_kinds, only : dp
   use fer_special, only : normal_cdf
   use fer_vanilla, only : black_scholes_price
   implicit none
   private
   public :: sabr_hagan_2002, sabr_hagan_price, nsvh1_choi_2019
contains
   elemental real(dp) function sabr_hagan_2002(strike, forward, texp, sigma, vov, rho, beta) result(vol_bs)
      real(dp), intent(in) :: strike, forward, texp, sigma, vov, rho, beta
      real(dp) :: betac, betac2, rho2, pfs, lf, lf2, pre1, pre2, z, y, xz
      betac = 1.0_dp-beta
      betac2 = betac*betac
      rho2 = rho*rho
      pfs = (forward*strike)**(0.5_dp*betac)
      lf = log(forward/strike)
      lf2 = lf*lf
      pre1 = pfs*(1.0_dp+betac2*lf2/24.0_dp*(1.0_dp+betac2*lf2/80.0_dp))
      pre2 = 1.0_dp+texp*((2.0_dp-3.0_dp*rho2)*vov*vov/24.0_dp + &
         sigma*(vov*rho*beta/(4.0_dp*pfs)+betac2*sigma/(24.0_dp*pfs*pfs)))
      if (abs(vov) <= epsilon(1.0_dp)) then
         xz = 1.0_dp
      else
         z = pfs*lf*vov/sigma
         y = sqrt(1.0_dp+z*(z-2.0_dp*rho))
         if (z >= 1.0e-5_dp) then
            xz = log((y+z-rho)/(1.0_dp-rho))/z
         else if (z <= -1.0e-5_dp) then
            xz = log((1.0_dp+rho)/(y-z+rho))/z
         else
            xz = 1.0_dp+0.5_dp*z*(rho+z*((rho2-1.0_dp/3.0_dp)+ &
               (5.0_dp*rho2-3.0_dp)*rho*z/4.0_dp))
         end if
      end if
      vol_bs = sigma*pre2/(pre1*xz)
   end function sabr_hagan_2002

   elemental real(dp) function sabr_hagan_price(strike,forward,texp,sigma,vov,rho,beta,cp,df) result(price)
      real(dp), intent(in) :: strike,forward,texp,sigma,vov,rho,beta,df
      integer, intent(in) :: cp
      price = black_scholes_price(strike,forward,texp,sabr_hagan_2002(strike,forward,texp,sigma,vov,rho,beta),cp,df)
   end function sabr_hagan_price

   elemental real(dp) function nsvh1_choi_2019(strike,forward,texp,sigma,vov,rho,cp,df) result(price)
      real(dp), intent(in) :: strike,forward,texp,sigma,vov,rho,df
      integer, intent(in) :: cp
      real(dp) :: rhoc, vsqt, vvar, d, p
      if (abs(vov) <= sqrt(epsilon(1.0_dp))) then
         price = df*max(real(cp,dp)*(forward-strike),0.0_dp)
         return
      end if
      rhoc = sqrt(1.0_dp-rho*rho)
      vsqt = vov*sqrt(texp)
      vvar = exp(0.5_dp*texp*vov*vov)
      d = (asinh(((forward-strike)*vov/sigma-vvar*rho)/rhoc)+atanh(rho))/vsqt
      p = (forward-strike-sigma*rho*vvar/vov)*normal_cdf(real(cp,dp)*d)
      p = p+0.5_dp*sigma*vvar/vov*((1.0_dp+rho)*normal_cdf(real(cp,dp)*(d+vsqt))- &
         (1.0_dp-rho)*normal_cdf(real(cp,dp)*(d-vsqt)))
      price = df*real(cp,dp)*p
   end function nsvh1_choi_2019
end module fer_sabr
