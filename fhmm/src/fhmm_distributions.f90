! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_distributions
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use fhmm_kinds, only: dp, pi, log_two_pi, tiny_prob
   use fhmm_types, only: dist_normal, dist_lognormal, dist_student_t, dist_gamma, dist_poisson
   use fhmm_math, only: normal_cdf, normal_quantile, regularized_gamma_p, regularized_beta
   use fhmm_math, only: random_normal, random_gamma, random_poisson
   implicit none
   private
   public :: distribution_pdf, distribution_logpdf, distribution_cdf
   public :: distribution_quantile, distribution_random, distribution_mean

contains

   pure real(dp) function distribution_logpdf(family,x,mu,sigma,df) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: x,mu,sigma,df
      real(dp) :: z,shape,scale
      integer :: k
      select case(family)
      case(dist_normal)
         if (sigma <= 0.0_dp) then
            value=log(tiny_prob)
         else
            z=(x-mu)/sigma
            value=-0.5_dp*log_two_pi-log(sigma)-0.5_dp*z*z
         end if
      case(dist_lognormal)
         if (x <= 0.0_dp .or. sigma <= 0.0_dp) then
            value=log(tiny_prob)
         else
            z=(log(x)-mu)/sigma
            value=-0.5_dp*log_two_pi-log(sigma)-log(x)-0.5_dp*z*z
         end if
      case(dist_student_t)
         if (sigma <= 0.0_dp .or. df <= 0.0_dp) then
            value=log(tiny_prob)
         else
            z=(x-mu)/sigma
            value=log_gamma(0.5_dp*(df+1.0_dp))-log_gamma(0.5_dp*df) &
               -0.5_dp*log(df*pi)-log(sigma) &
               -0.5_dp*(df+1.0_dp)*log(1.0_dp+z*z/df)
         end if
      case(dist_gamma)
         if (mu <= 0.0_dp .or. sigma <= 0.0_dp .or. x <= 0.0_dp) then
            value=log(tiny_prob)
         else
            shape=mu*mu/(sigma*sigma)
            scale=sigma*sigma/mu
            value=(shape-1.0_dp)*log(x)-x/scale-log_gamma(shape)-shape*log(scale)
         end if
      case(dist_poisson)
         if (mu <= 0.0_dp .or. x < 0.0_dp) then
            value=log(tiny_prob)
         else
            k=nint(x)
            if (abs(x-real(k,dp)) > 1.0e-8_dp) then
               value=log(tiny_prob)
            else
               value=real(k,dp)*log(mu)-mu-log_gamma(real(k+1,dp))
            end if
         end if
      case default
         value=log(tiny_prob)
      end select
      value=max(value,log(tiny_prob))
   end function distribution_logpdf

   pure real(dp) function distribution_pdf(family,x,mu,sigma,df) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: x,mu,sigma,df
      value=exp(distribution_logpdf(family,x,mu,sigma,df))
   end function distribution_pdf

   pure real(dp) function distribution_cdf(family,x,mu,sigma,df) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: x,mu,sigma,df
      real(dp) :: z,shape,scale,y,pbeta
      integer :: k,j
      select case(family)
      case(dist_normal)
         if (sigma <= 0.0_dp) then
            value=merge(1.0_dp,0.0_dp,x>=mu)
         else
            value=normal_cdf((x-mu)/sigma)
         end if
      case(dist_lognormal)
         if (x <= 0.0_dp) then
            value=0.0_dp
         else if (sigma <= 0.0_dp) then
            value=merge(1.0_dp,0.0_dp,x>=exp(mu))
         else
            value=normal_cdf((log(x)-mu)/sigma)
         end if
      case(dist_student_t)
         if (sigma <= 0.0_dp .or. df <= 0.0_dp) then
            value=0.0_dp
         else
            z=(x-mu)/sigma
            y=df/(df+z*z)
            pbeta=regularized_beta(y,0.5_dp*df,0.5_dp)
            value=merge(1.0_dp-0.5_dp*pbeta,0.5_dp*pbeta,z>=0.0_dp)
         end if
      case(dist_gamma)
         if (x <= 0.0_dp) then
            value=0.0_dp
         else if (mu <= 0.0_dp .or. sigma <= 0.0_dp) then
            value=0.0_dp
         else
            shape=mu*mu/(sigma*sigma)
            scale=sigma*sigma/mu
            value=regularized_gamma_p(shape,x/scale)
         end if
      case(dist_poisson)
         if (x < 0.0_dp .or. mu <= 0.0_dp) then
            value=0.0_dp
         else
            k=floor(x)
            value=0.0_dp
            do j=0,k
               value=value+exp(real(j,dp)*log(mu)-mu-log_gamma(real(j+1,dp)))
            end do
         end if
      case default
         value=0.0_dp
      end select
      value=min(max(value,0.0_dp),1.0_dp)
   end function distribution_cdf

   pure real(dp) function distribution_quantile(family,p,mu,sigma,df) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: p,mu,sigma,df
      real(dp) :: lo,hi,mid,c
      integer :: iter,k
      if (p <= 0.0_dp) then
         select case(family)
         case(dist_gamma,dist_lognormal,dist_poisson); value=0.0_dp
         case default; value=-huge(1.0_dp)
         end select
         return
      else if (p >= 1.0_dp) then
         value=huge(1.0_dp); return
      end if
      select case(family)
      case(dist_normal)
         value=mu+sigma*normal_quantile(p)
      case(dist_lognormal)
         value=exp(mu+sigma*normal_quantile(p))
      case(dist_poisson)
         k=0
         do while (distribution_cdf(family,real(k,dp),mu,sigma,df) < p .and. k < 1000000)
            k=k+1
         end do
         value=real(k,dp)
      case(dist_student_t)
         lo=mu-100.0_dp*max(sigma,1.0_dp); hi=mu+100.0_dp*max(sigma,1.0_dp)
         do iter=1,200
            mid=0.5_dp*(lo+hi); c=distribution_cdf(family,mid,mu,sigma,df)
            if(c<p) then; lo=mid; else; hi=mid; end if
         end do
         value=0.5_dp*(lo+hi)
      case(dist_gamma)
         lo=0.0_dp; hi=max(mu+10.0_dp*sigma,1.0_dp)
         do while(distribution_cdf(family,hi,mu,sigma,df)<p)
            hi=2.0_dp*hi
         end do
         do iter=1,200
            mid=0.5_dp*(lo+hi); c=distribution_cdf(family,mid,mu,sigma,df)
            if(c<p) then; lo=mid; else; hi=mid; end if
         end do
         value=0.5_dp*(lo+hi)
      case default
         value=0.0_dp
      end select
   end function distribution_quantile

   real(dp) function distribution_random(family,mu,sigma,df) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: mu,sigma,df
      select case(family)
      case(dist_normal)
         value=mu+sigma*random_normal()
      case(dist_lognormal)
         value=exp(mu+sigma*random_normal())
      case(dist_student_t)
         value=mu+sigma*random_normal()/sqrt(random_gamma(0.5_dp*df,2.0_dp)/df)
      case(dist_gamma)
         value=random_gamma(mu*mu/(sigma*sigma),sigma*sigma/mu)
      case(dist_poisson)
         value=real(random_poisson(mu),dp)
      case default
         value=0.0_dp
      end select
   end function distribution_random

   pure real(dp) function distribution_mean(family,mu,sigma,df) result(value)
      integer, intent(in) :: family
      real(dp), intent(in) :: mu,sigma,df
      if (family == dist_student_t .and. df <= 1.0_dp) then
         value=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      select case(family)
      case(dist_lognormal)
         value=exp(mu+0.5_dp*sigma*sigma)
      case default
         value=mu
      end select
   end function distribution_mean

end module fhmm_distributions
