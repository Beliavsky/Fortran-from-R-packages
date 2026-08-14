module randtoolbox_math
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   public :: normal_quantile, gamma_q, chi_square_survival
contains
   pure elemental real(real64) function normal_quantile(p) result(x)
      real(real64), intent(in) :: p
      real(real64), parameter :: a(6)=[-3.969683028665376e1_real64,2.209460984245205e2_real64, &
         -2.759285104469687e2_real64,1.383577518672690e2_real64,-3.066479806614716e1_real64,2.506628277459239_real64]
      real(real64), parameter :: b(5)=[-5.447609879822406e1_real64,1.615858368580409e2_real64, &
         -1.556989798598866e2_real64,6.680131188771972e1_real64,-1.328068155288572e1_real64]
      real(real64), parameter :: c(6)=[-7.784894002430293e-3_real64,-3.223964580411365e-1_real64, &
         -2.400758277161838_real64,-2.549732539343734_real64,4.374664141464968_real64,2.938163982698783_real64]
      real(real64), parameter :: d(4)=[7.784695709041462e-3_real64,3.224671290700398e-1_real64, &
         2.445134137142996_real64,3.754408661907416_real64]
      real(real64), parameter :: plow=0.02425_real64, phigh=1.0_real64-plow
      real(real64) :: q,r,e,u
      if(p<=0.0_real64) then
         x=ieee_value(x,ieee_negative_inf); return
      else if(p>=1.0_real64) then
         x=ieee_value(x,ieee_positive_inf); return
      end if
      if(p<plow) then
         q=sqrt(-2.0_real64*log(p))
         x=(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_real64)
      else if(p<=phigh) then
         q=p-0.5_real64; r=q*q
         x=(((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
           (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_real64)
      else
         q=sqrt(-2.0_real64*log(1.0_real64-p))
         x=-(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_real64)
      end if
      ! One Halley refinement gives near-double precision.
      e=0.5_real64*erfc(-x/sqrt(2.0_real64))-p
      u=e*sqrt(2.0_real64*acos(-1.0_real64))*exp(0.5_real64*x*x)
      x=x-u/(1.0_real64+0.5_real64*x*u)
   end function normal_quantile

   pure real(real64) function gamma_q(a,x) result(q)
      real(real64), intent(in) :: a,x
      integer, parameter :: itmax=10000
      real(real64), parameter :: eps=8.0_real64*epsilon(1.0_real64), fpmin=tiny(1.0_real64)/eps
      integer :: n
      real(real64) :: ap,del,sum,b,c,d,h,an,gln,p
      if(a<=0.0_real64 .or. x<0.0_real64) then
         q=0.0_real64; return
      end if
      if(x<=0.0_real64) then
         q=1.0_real64; return
      end if
      gln=log_gamma(a)
      if(x<a+1.0_real64) then
         ap=a; sum=1.0_real64/a; del=sum
         do n=1,itmax
            ap=ap+1.0_real64; del=del*x/ap; sum=sum+del
            if(abs(del)<=abs(sum)*eps) exit
         end do
         p=sum*exp(-x+a*log(x)-gln)
         q=max(0.0_real64,min(1.0_real64,1.0_real64-p))
      else
         b=x+1.0_real64-a; c=1.0_real64/fpmin; d=1.0_real64/b; h=d
         do n=1,itmax
            an=-real(n,real64)*(real(n,real64)-a)
            b=b+2.0_real64; d=an*d+b
            if(abs(d)<fpmin)d=fpmin
            c=b+an/c; if(abs(c)<fpmin)c=fpmin
            d=1.0_real64/d; del=d*c; h=h*del
            if(abs(del-1.0_real64)<=eps) exit
         end do
         q=h*exp(-x+a*log(x)-gln)
         q=max(0.0_real64,min(1.0_real64,q))
      end if
   end function gamma_q

   pure real(real64) function chi_square_survival(x,df) result(p)
      real(real64), intent(in) :: x
      integer, intent(in) :: df
      if(df<=0) then
         p=0.0_real64
      else
         p=gamma_q(0.5_real64*real(df,real64),0.5_real64*x)
      end if
   end function chi_square_survival
end module randtoolbox_math
