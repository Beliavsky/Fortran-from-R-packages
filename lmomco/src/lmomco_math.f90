module lmomco_math
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use lmomco_kinds, only : dp, pi
   implicit none
   private
   public :: normal_cdf, normal_pdf, normal_quantile
   public :: gamma_p, gamma_q, gamma_quantile
   public :: beta_inc, student_t_cdf, student_t_quantile
   public :: bessel_i_nu, integrate_simpson, clamp01

contains

   pure elemental real(dp) function clamp01(x) result(y)
      real(dp), intent(in) :: x
      y = max(0.0_dp, min(1.0_dp, x))
   end function clamp01

   pure elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: q, r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, -2.759285104469687e2_dp, &
          1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, -1.556989798598866e2_dp, &
          6.680131188771972e1_dp, -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, -2.400758277161838_dp, &
         -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, 2.445134137142996_dp, &
          3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
   end function normal_quantile

   pure real(dp) function gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 2.0e-15_dp, fpmin = 1.0e-300_dp
      integer :: n
      real(dp) :: sumv, del, ap, b, c, d, h, an, gln
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         p = 0.0_dp
         return
      end if
      gln = log_gamma(a)
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x + a*log(x) - gln)
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/max(b, fpmin)
         h = d
         do n = 1, itmax
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         p = 1.0_dp - exp(-x + a*log(x) - gln)*h
      end if
      p = clamp01(p)
   end function gamma_p

   pure real(dp) function gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      q = 1.0_dp - gamma_p(a, x)
   end function gamma_q

   real(dp) function gamma_quantile(p, shape) result(x)
      real(dp), intent(in) :: p, shape
      real(dp) :: lo, hi, mid
      integer :: i
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      end if
      if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      lo = 0.0_dp
      hi = max(1.0_dp, shape)
      do while (gamma_p(shape, hi) < p .and. hi < huge(1.0_dp)/4.0_dp)
         hi = hi*2.0_dp
      end do
      do i = 1, 160
         mid = 0.5_dp*(lo+hi)
         if (gamma_p(shape, mid) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
   end function gamma_quantile

   pure real(dp) function betacf(a,b,x) result(cf)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=300
      real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
      integer :: m, m2
      real(dp) :: qab,qap,qam,c,d,h,aa,del
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp
      d=1.0_dp-qab*x/qap
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

   pure real(dp) function beta_inc(x,a,b) result(p)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else if (x >= 1.0_dp) then
         p = 1.0_dp
      else
         bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
         if(x < (a+1.0_dp)/(a+b+2.0_dp))then
            p=bt*betacf(a,b,x)/a
         else
            p=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
         end if
         p=clamp01(p)
      end if
   end function beta_inc

   pure real(dp) function student_t_cdf(t, nu) result(p)
      real(dp), intent(in) :: t, nu
      real(dp) :: z, ib
      if (nu <= 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (abs(t) <= tiny(1.0_dp)) then
         p = 0.5_dp
         return
      end if
      z = nu/(nu+t*t)
      ib = beta_inc(z, 0.5_dp*nu, 0.5_dp)
      if (t > 0.0_dp) then
         p = 1.0_dp - 0.5_dp*ib
      else
         p = 0.5_dp*ib
      end if
   end function student_t_cdf

   real(dp) function student_t_quantile(p, nu) result(x)
      real(dp), intent(in) :: p, nu
      real(dp) :: lo, hi, mid
      integer :: i
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp); return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp); return
      end if
      lo = -1.0_dp; hi = 1.0_dp
      do while(student_t_cdf(lo,nu) > p)
         lo=2.0_dp*lo
      end do
      do while(student_t_cdf(hi,nu) < p)
         hi=2.0_dp*hi
      end do
      do i=1,180
         mid=0.5_dp*(lo+hi)
         if(student_t_cdf(mid,nu) < p)then
            lo=mid
         else
            hi=mid
         end if
      end do
      x=0.5_dp*(lo+hi)
   end function student_t_quantile

   pure real(dp) function bessel_i_nu(x, nu) result(v)
      real(dp), intent(in) :: x, nu
      integer :: k
      real(dp) :: term, sumv, halfx, lg
      if (x < 0.0_dp) then
         v = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         if (abs(nu) < 1.0e-14_dp) then
            v = 1.0_dp
         else if (nu > 0.0_dp) then
            v = 0.0_dp
         else
            v = huge(1.0_dp)
         end if
         return
      end if
      if (x > 80.0_dp) then
         v = exp(x)/sqrt(2.0_dp*pi*x) * &
             (1.0_dp - (4.0_dp*nu*nu-1.0_dp)/(8.0_dp*x) + &
             (4.0_dp*nu*nu-1.0_dp)*(4.0_dp*nu*nu-9.0_dp)/(2.0_dp*(8.0_dp*x)**2))
         return
      end if
      halfx = 0.5_dp*x
      lg = log_gamma(nu+1.0_dp)
      term = exp(nu*log(halfx)-lg)
      sumv = term
      do k=1,10000
         term = term*(halfx*halfx)/(real(k,dp)*(nu+real(k,dp)))
         sumv = sumv + term
         if(abs(term) <= 2.0e-15_dp*max(1.0_dp,abs(sumv))) exit
      end do
      v = sumv
   end function bessel_i_nu

   real(dp) function integrate_simpson(fun, a, b, n) result(v)
      interface
         function fun(x) result(y)
            import dp
            real(dp), intent(in) :: x
            real(dp) :: y
         end function fun
      end interface
      real(dp), intent(in) :: a,b
      integer, intent(in), optional :: n
      integer :: nn, i
      real(dp) :: h, s, x
      if (b <= a) then
         v=0.0_dp; return
      end if
      nn=800
      if(present(n)) nn=max(20,n)
      if(mod(nn,2)/=0) nn=nn+1
      h=(b-a)/real(nn,dp)
      s=fun(a)+fun(b)
      do i=1,nn-1
         x=a+h*real(i,dp)
         if(mod(i,2)==0)then
            s=s+2.0_dp*fun(x)
         else
            s=s+4.0_dp*fun(x)
         end if
      end do
      v=s*h/3.0_dp
   end function integrate_simpson

end module lmomco_math
