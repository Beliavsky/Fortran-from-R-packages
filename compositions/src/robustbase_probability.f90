! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_probability
   use robustbase_kinds, only: dp
   implicit none
   private
   public :: normal_cdf, normal_quantile, chi_square_quantile, chi_square_cdf, regularized_gamma_p
contains
   elemental function normal_cdf(x) result(p)
      real(dp),intent(in)::x;real(dp)::p
      p=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function
   elemental function normal_quantile(p) result(x)
      real(dp),intent(in)::p;real(dp)::x,q,r
      real(dp),parameter::a(6)=[-3.969683028665376e1_dp,2.209460984245205e2_dp, &
         -2.759285104469687e2_dp,1.383577518672690e2_dp,-3.066479806614716e1_dp, &
         2.506628277459239_dp]
      real(dp),parameter::b(5)=[-5.447609879822406e1_dp,1.615858368580409e2_dp, &
         -1.556989798598866e2_dp,6.680131188771972e1_dp,-1.328068155288572e1_dp]
      real(dp),parameter::c(6)=[-7.784894002430293e-3_dp,-3.223964580411365e-1_dp, &
         -2.400758277161838_dp,-2.549732539343734_dp,4.374664141464968_dp, &
         2.938163982698783_dp]
      real(dp),parameter::d(4)=[7.784695709041462e-3_dp,3.224671290700398e-1_dp,2.445134137142996_dp,3.754408661907416_dp]
      if(p<=0.0_dp) then;x=-huge(1.0_dp);return;else if(p>=1.0_dp) then;x=huge(1.0_dp);return;end if
      if(p<0.02425_dp) then
         q=sqrt(-2.0_dp*log(p));x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if(p>0.97575_dp) then
         q=sqrt(-2.0_dp*log(1.0_dp-p))
         x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else
         q=p-0.5_dp;r=q*q;x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/(((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      end if
   end function

   function regularized_gamma_p(a,x) result(p)
      real(dp),intent(in)::a,x
      real(dp)::p,ap,sumv,del,b,c,d,h,an
      integer::n
      real(dp),parameter::eps=2.0e-14_dp,fpmin=1.0e-300_dp
      if(a<=0.0_dp .or. x<0.0_dp)error stop "regularized_gamma_p: invalid arguments"
      if(x<=tiny(1.0_dp))then
         p=0.0_dp
         return
      end if
      if(x<a+1.0_dp)then
         ap=a;sumv=1.0_dp/a;del=sumv
         do n=1,10000
            ap=ap+1.0_dp;del=del*x/ap;sumv=sumv+del
            if(abs(del)<=abs(sumv)*eps)exit
         end do
         p=sumv*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1.0_dp-a;c=1.0_dp/fpmin;d=1.0_dp/max(b,fpmin);h=d
         do n=1,10000
            an=-real(n,dp)*(real(n,dp)-a)
            b=b+2.0_dp;d=an*d+b;if(abs(d)<fpmin)d=fpmin
            c=b+an/c;if(abs(c)<fpmin)c=fpmin
            d=1.0_dp/d;del=d*c;h=h*del
            if(abs(del-1.0_dp)<=eps)exit
         end do
         p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      p=max(0.0_dp,min(1.0_dp,p))
   end function regularized_gamma_p

   function chi_square_cdf(x,df) result(p)
      real(dp),intent(in)::x,df
      real(dp)::p
      if(df<=0.0_dp)error stop "chi_square_cdf: df must be positive"
      if(x<=0.0_dp)then
         p=0.0_dp
      else
         p=regularized_gamma_p(0.5_dp*df,0.5_dp*x)
      end if
   end function chi_square_cdf

   function chi_square_quantile(p,df) result(q)
      real(dp),intent(in)::p,df;real(dp)::q,z
      z=normal_quantile(p)
      q=df*(1.0_dp-2.0_dp/(9.0_dp*df)+z*sqrt(2.0_dp/(9.0_dp*df)))**3
      q=max(q,0.0_dp)
   end function
end module robustbase_probability
