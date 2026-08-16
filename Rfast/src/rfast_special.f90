module rfast_special
   use, intrinsic :: iso_fortran_env, only : real64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private
   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt2 = sqrt(2.0_dp)

   public :: log_beta, log_choose, choose_real
   public :: digamma_r, trigamma_r
   public :: normal_pdf, normal_cdf, normal_quantile
   public :: reg_gamma_p, reg_gamma_q, reg_beta
   public :: chisq_cdf, student_t_cdf, f_cdf
   public :: log1p_r, expm1_r, nan_r

contains

   pure elemental real(dp) function nan_r() result(v)
      v = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_r

   pure elemental real(dp) function log1p_r(x) result(v)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-8_dp) then
         v = x - 0.5_dp*x*x + x*x*x/3.0_dp - 0.25_dp*x**4
      else
         v = log(1.0_dp + x)
      end if
   end function log1p_r

   pure elemental real(dp) function expm1_r(x) result(v)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-8_dp) then
         v = x + 0.5_dp*x*x + x*x*x/6.0_dp + x**4/24.0_dp
      else
         v = exp(x) - 1.0_dp
      end if
   end function expm1_r

   pure elemental real(dp) function log_beta(a, b) result(v)
      real(dp), intent(in) :: a, b
      v = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
   end function log_beta

   pure elemental real(dp) function log_choose(n, k) result(v)
      real(dp), intent(in) :: n, k
      if (k < 0.0_dp .or. k > n) then
         v = -huge(1.0_dp)
      else
         v = log_gamma(n + 1.0_dp) - log_gamma(k + 1.0_dp) - log_gamma(n - k + 1.0_dp)
      end if
   end function log_choose

   pure elemental real(dp) function choose_real(n, k) result(v)
      real(dp), intent(in) :: n, k
      real(dp) :: lc
      lc = log_choose(n, k)
      if (lc <= -0.5_dp*huge(1.0_dp)) then
         v = 0.0_dp
      else
         v = exp(lc)
      end if
   end function choose_real

   pure recursive elemental real(dp) function digamma_r(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, inv, inv2
      if (x <= 0.0_dp .and. abs(x - anint(x)) < 1.0e-14_dp) then
         y = huge(1.0_dp)
         return
      end if
      if (x < 0.0_dp) then
         y = digamma_r(1.0_dp - x) - pi / tan(pi*x)
         return
      end if
      z = x
      y = 0.0_dp
      do while (z < 8.0_dp)
         y = y - 1.0_dp/z
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z
      inv2 = inv*inv
      y = y + log(z) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
          inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp - inv2*5.0_dp/660.0_dp))))
   end function digamma_r

   pure recursive elemental real(dp) function trigamma_r(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, inv, inv2
      if (x <= 0.0_dp .and. abs(x - anint(x)) < 1.0e-14_dp) then
         y = huge(1.0_dp)
         return
      end if
      if (x < 0.0_dp) then
         y = pi*pi/(sin(pi*x)**2) - trigamma_r(1.0_dp - x)
         return
      end if
      z = x
      y = 0.0_dp
      do while (z < 8.0_dp)
         y = y + 1.0_dp/(z*z)
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z
      inv2 = inv*inv
      y = y + inv + 0.5_dp*inv2 + inv*inv2*(1.0_dp/6.0_dp - inv2*(1.0_dp/30.0_dp - &
          inv2*(1.0_dp/42.0_dp - inv2*(1.0_dp/30.0_dp - inv2*5.0_dp/66.0_dp))))
   end function trigamma_r

   pure elemental real(dp) function normal_pdf(x) result(v)
      real(dp), intent(in) :: x
      v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental real(dp) function normal_cdf(x) result(v)
      real(dp), intent(in) :: x
      v = 0.5_dp*erfc(-x/sqrt2)
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376e1_dp, a2=2.209460984245205e2_dp
      real(dp), parameter :: a3=-2.759285104469687e2_dp, a4=1.383577518672690e2_dp
      real(dp), parameter :: a5=-3.066479806614716e1_dp, a6=2.506628277459239_dp
      real(dp), parameter :: b1=-5.447609879822406e1_dp, b2=1.615858368580409e2_dp
      real(dp), parameter :: b3=-1.556989798598866e2_dp, b4=6.680131188771972e1_dp
      real(dp), parameter :: b5=-1.328068155288572e1_dp
      real(dp), parameter :: c1=-7.784894002430293e-3_dp, c2=-3.223964580411365e-1_dp
      real(dp), parameter :: c3=-2.400758277161838_dp, c4=-2.549732539343734_dp
      real(dp), parameter :: c5=4.374664141464968_dp, c6=2.938163982698783_dp
      real(dp), parameter :: d1=7.784695709041462e-3_dp, d2=3.224671290700398e-1_dp
      real(dp), parameter :: d3=2.445134137142996_dp, d4=3.754408661907416_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
             ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p > phigh) then
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
              ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q = p - 0.5_dp
         r = q*q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
             (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
      if (abs(x) < huge(1.0_dp)/2.0_dp) then
         x = x - (normal_cdf(x)-p)/normal_pdf(x)
      end if
   end function normal_quantile

   pure real(dp) function reg_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax=10000
      real(dp), parameter :: eps=2.0e-15_dp, fpmin=tiny(1.0_dp)/eps
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = nan_r()
         return
      end if
      if (x == 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n=1,itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n=1,itmax
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         p = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
      end if
      p = max(0.0_dp,min(1.0_dp,p))
   end function reg_gamma_p

   pure real(dp) function reg_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      q = 1.0_dp - reg_gamma_p(a,x)
      q = max(0.0_dp,min(1.0_dp,q))
   end function reg_gamma_q

   pure real(dp) function betacf(a,b,x) result(cf)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=10000
      real(dp), parameter :: eps=2.0e-15_dp, fpmin=tiny(1.0_dp)/eps
      integer :: m, m2
      real(dp) :: aa,c,d,del,h,qab,qam,qap
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp; d=1.0_dp-qab*x/qap
      if(abs(d)<fpmin)d=fpmin
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
         if(abs(del-1.0_dp)<=eps)exit
      end do
      cf=h
   end function betacf

   pure real(dp) function reg_beta(x,a,b) result(v)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v=nan_r(); return
      end if
      if(x<=0.0_dp)then; v=0.0_dp; return; end if
      if(x>=1.0_dp)then; v=1.0_dp; return; end if
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log1p_r(-x))
      if(x<(a+1.0_dp)/(a+b+2.0_dp))then
         v=bt*betacf(a,b,x)/a
      else
         v=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
      end if
      v=max(0.0_dp,min(1.0_dp,v))
   end function reg_beta

   pure elemental real(dp) function chisq_cdf(x, df) result(p)
      real(dp), intent(in) :: x, df
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         p = reg_gamma_p(0.5_dp*df,0.5_dp*x)
      end if
   end function chisq_cdf

   pure elemental real(dp) function student_t_cdf(t, df) result(p)
      real(dp), intent(in) :: t, df
      real(dp) :: x, ib
      if (df <= 0.0_dp) then
         p=nan_r(); return
      end if
      if (t == 0.0_dp) then
         p=0.5_dp; return
      end if
      x=df/(df+t*t)
      ib=reg_beta(x,0.5_dp*df,0.5_dp)
      if(t>0.0_dp)then
         p=1.0_dp-0.5_dp*ib
      else
         p=0.5_dp*ib
      end if
   end function student_t_cdf

   pure elemental real(dp) function f_cdf(x, df1, df2) result(p)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: z
      if (x <= 0.0_dp) then
         p=0.0_dp
      else
         z=(df1*x)/(df1*x+df2)
         p=reg_beta(z,0.5_dp*df1,0.5_dp*df2)
      end if
   end function f_cdf

end module rfast_special
