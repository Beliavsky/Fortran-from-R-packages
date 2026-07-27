! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_european
   use optionpricing_kinds, only : dp
   use optionpricing_math, only : normal_cdf, normal_pdf
   use optionpricing_types, only : european_result
   implicit none
   private
   public :: bs_european_call, bs_european_put, bs_ec, bs_ep
contains
   pure function bs_european_call(t, k, r, sigma, s0) result(res)
      real(dp), intent(in) :: t, k, r, sigma, s0
      type(european_result) :: res
      real(dp) :: d1, d2, root_t, disc
      if (t <= 0.0_dp) then
         res%price=max(s0-k,0.0_dp)
         res%delta=merge(1.0_dp,0.0_dp,s0>k)
         return
      end if
      if (sigma <= 0.0_dp) then
         disc=exp(-r*t)
         res%price=max(s0-k*disc,0.0_dp)
         res%delta=merge(1.0_dp,0.0_dp,s0>k*disc)
         return
      end if
      root_t=sqrt(t)
      d1=(log(s0/k)+(r+0.5_dp*sigma*sigma)*t)/(sigma*root_t)
      d2=d1-sigma*root_t
      disc=exp(-r*t)
      res%price=s0*normal_cdf(d1)-k*disc*normal_cdf(d2)
      res%delta=normal_cdf(d1)
      res%gamma=normal_pdf(d1)/(s0*sigma*root_t)
      ! Upstream R used Phi(d1) instead of phi(d1).
      res%upstream_gamma=res%delta/(s0*sigma*root_t)
   end function bs_european_call

   pure function bs_european_put(t, k, r, sigma, s0) result(res)
      real(dp), intent(in) :: t, k, r, sigma, s0
      type(european_result) :: res
      real(dp) :: d1, d2, root_t, disc
      if (t <= 0.0_dp) then
         res%price=max(k-s0,0.0_dp)
         res%delta=merge(-1.0_dp,0.0_dp,s0<k)
         return
      end if
      if (sigma <= 0.0_dp) then
         disc=exp(-r*t)
         res%price=max(k*disc-s0,0.0_dp)
         res%delta=merge(-1.0_dp,0.0_dp,s0<k*disc)
         return
      end if
      root_t=sqrt(t)
      d1=(log(s0/k)+(r+0.5_dp*sigma*sigma)*t)/(sigma*root_t)
      d2=d1-sigma*root_t
      disc=exp(-r*t)
      res%price=-s0*normal_cdf(-d1)+k*disc*normal_cdf(-d2)
      res%delta=-normal_cdf(-d1)
      res%gamma=normal_pdf(d1)/(s0*sigma*root_t)
      res%upstream_gamma=normal_cdf(d1)/(s0*sigma*root_t)
   end function bs_european_put

   pure function bs_ec(t,k,r,sigma,s0) result(res)
      real(dp), intent(in) :: t,k,r,sigma,s0
      type(european_result) :: res
      res=bs_european_call(t,k,r,sigma,s0)
   end function bs_ec

   pure function bs_ep(t,k,r,sigma,s0) result(res)
      real(dp), intent(in) :: t,k,r,sigma,s0
      type(european_result) :: res
      res=bs_european_put(t,k,r,sigma,s0)
   end function bs_ep
end module optionpricing_european
