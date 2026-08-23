module psdistr_math
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use psdistr_kinds, only : dp
   implicit none
   private
   public :: normal_pdf, normal_cdf, normal_quantile, signed_root
   public :: regularized_beta, beta_quantile

contains

   elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp), parameter :: inv_sqrt_2pi = 0.3989422804014326779399460599343819_dp
      y = inv_sqrt_2pi * exp(-0.5_dp*x*x)
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      y = 0.5_dp * erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function signed_root(x, n) result(y)
      real(dp), intent(in) :: x, n
      if (n <= 0.0_dp) then
         y = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (abs(x) <= tiny(1.0_dp)) then
         y = 0.0_dp
      else
         y = sign(abs(x)**(1.0_dp/n), x)
      end if
   end function signed_root

   elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: q, r, e
      real(dp), parameter :: a1=-3.969683028665376e+01_dp, a2= 2.209460984245205e+02_dp
      real(dp), parameter :: a3=-2.759285104469687e+02_dp, a4= 1.383577518672690e+02_dp
      real(dp), parameter :: a5=-3.066479806614716e+01_dp, a6= 2.506628277459239e+00_dp
      real(dp), parameter :: b1=-5.447609879822406e+01_dp, b2= 1.615858368580409e+02_dp
      real(dp), parameter :: b3=-1.556989798598866e+02_dp, b4= 6.680131188771972e+01_dp
      real(dp), parameter :: b5=-1.328068155288572e+01_dp
      real(dp), parameter :: c1=-7.784894002430293e-03_dp, c2=-3.223964580411365e-01_dp
      real(dp), parameter :: c3=-2.400758277161838e+00_dp, c4=-2.549732539343734e+00_dp
      real(dp), parameter :: c5= 4.374664141464968e+00_dp, c6= 2.938163982698783e+00_dp
      real(dp), parameter :: d1= 7.784695709041462e-03_dp, d2= 3.224671290700398e-01_dp
      real(dp), parameter :: d3= 2.445134137142996e+00_dp, d4= 3.754408661907416e+00_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      if (p < plow) then
         q=sqrt(-2.0_dp*log(p))
         x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p <= phigh) then
         q=p-0.5_dp; r=q*q
         x=(((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q=sqrt(-2.0_dp*log(1.0_dp-p))
         x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
      ! One Halley correction is enough for full double precision in the central range.
      if (abs(x) < 8.0_dp) then
         e = normal_cdf(x) - p
         x = x - e / max(normal_pdf(x), tiny(1.0_dp))
      end if
   end function normal_quantile

   pure real(dp) function beta_cf(a,b,x) result(cf)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=300
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      integer :: m, m2
      real(dp) :: aa,c,d,del,h,qab,qam,qap
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp
      d=1.0_dp-qab*x/qap
      if(abs(d)<fpmin) d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
         m2=2*m
         aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; h=h*d*c
         aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; del=d*c; h=h*del
         if(abs(del-1.0_dp)<=eps) exit
      end do
      cf=h
   end function beta_cf

   pure real(dp) function regularized_beta(x,a,b) result(y)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         y=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      if(x<=0.0_dp) then; y=0.0_dp; return; end if
      if(x>=1.0_dp) then; y=1.0_dp; return; end if
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if(x < (a+1.0_dp)/(a+b+2.0_dp)) then
         y=bt*beta_cf(a,b,x)/a
      else
         y=1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
      end if
      y=min(1.0_dp,max(0.0_dp,y))
   end function regularized_beta

   pure real(dp) function beta_quantile(p,a,b) result(x)
      real(dp), intent(in) :: p,a,b
      real(dp) :: lo,hi,mid
      integer :: i
      if(p<=0.0_dp) then; x=0.0_dp; return; end if
      if(p>=1.0_dp) then; x=1.0_dp; return; end if
      lo=0.0_dp; hi=1.0_dp
      do i=1,120
         mid=0.5_dp*(lo+hi)
         if(regularized_beta(mid,a,b)<p) then
            lo=mid
         else
            hi=mid
         end if
      end do
      x=0.5_dp*(lo+hi)
   end function beta_quantile
end module psdistr_math
