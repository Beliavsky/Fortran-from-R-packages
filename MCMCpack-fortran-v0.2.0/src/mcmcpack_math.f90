! SPDX-License-Identifier: GPL-3.0-only
module mcmcpack_math
   use mcmcpack_kinds, only : dp,pi
   implicit none
   private
   public :: logistic, log1pexp, logsumexp, normal_pdf, normal_cdf, normal_quantile
   public :: log_beta_fn, log_choose, log_factorial, log1mexp
contains
   elemental real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         p=1.0_dp/(1.0_dp+exp(-x))
      else
         p=exp(x)/(1.0_dp+exp(x))
      end if
   end function logistic

   elemental real(dp) function log1pexp(x) result(v)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then; v=x+log(1.0_dp+exp(-x)); else; v=log(1.0_dp+exp(x)); end if
   end function log1pexp

   pure real(dp) function logsumexp(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x)==0) then; v=-huge(1.0_dp); return; end if
      m=maxval(x); v=m+log(sum(exp(x-m)))
   end function logsumexp

   elemental real(dp) function normal_pdf(x) result(v)
      real(dp), intent(in) :: x
      v=exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(v)
      real(dp), intent(in) :: x
      v=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function normal_quantile(p) result(x)
      ! Acklam rational approximation, followed by one Newton step.
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376d1,a2=2.209460984245205d2
      real(dp), parameter :: a3=-2.759285104469687d2,a4=1.383577518672690d2
      real(dp), parameter :: a5=-3.066479806614716d1,a6=2.506628277459239d0
      real(dp), parameter :: b1=-5.447609879822406d1,b2=1.615858368580409d2
      real(dp), parameter :: b3=-1.556989798598866d2,b4=6.680131188771972d1,b5=-1.328068155288572d1
      real(dp), parameter :: c1=-7.784894002430293d-3,c2=-3.223964580411365d-1
      real(dp), parameter :: c3=-2.400758277161838d0,c4=-2.549732539343734d0
      real(dp), parameter :: c5=4.374664141464968d0,c6=2.938163982698783d0
      real(dp), parameter :: d1=7.784695709041462d-3,d2=3.224671290700398d-1
      real(dp), parameter :: d3=2.445134137142996d0,d4=3.754408661907416d0
      real(dp), parameter :: pl=0.02425_dp, ph=1.0_dp-pl
      real(dp) :: q,r,e
      if (p <= 0.0_dp) then; x=-huge(1.0_dp); return; end if
      if (p >= 1.0_dp) then; x=huge(1.0_dp); return; end if
      if (p < pl) then
         q=sqrt(-2.0_dp*log(p))
         x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p > ph) then
         q=sqrt(-2.0_dp*log(1.0_dp-p))
         x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q=p-0.5_dp; r=q*q
         x=((((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q)/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
      e=normal_cdf(x)-p
      x=x-e/max(normal_pdf(x),tiny(1.0_dp))
   end function normal_quantile

   elemental real(dp) function log_beta_fn(a,b) result(v)
      real(dp),intent(in)::a,b
      v=log_gamma(a)+log_gamma(b)-log_gamma(a+b)
   end function log_beta_fn

   elemental real(dp) function log_factorial(n) result(v)
      integer,intent(in)::n
      if(n<0) then; v=huge(1.0_dp); else; v=log_gamma(real(n+1,dp)); end if
   end function log_factorial

   elemental real(dp) function log_choose(n,k) result(v)
      integer,intent(in)::n,k
      if(k<0 .or. k>n) then; v=-huge(1.0_dp); else
         v=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
      end if
   end function log_choose

   elemental real(dp) function log1mexp(x) result(v)
      real(dp),intent(in)::x
      if (x >= 0.0_dp) then; v=-huge(1.0_dp)
      else if (x < log(0.5_dp)) then; v=log(1.0_dp-exp(x))
      else; v=log(1.0_dp-exp(x)); end if
   end function log1mexp
end module mcmcpack_math
