! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp, pi, sqrt2, log2pi
   implicit none
   private
   public :: log1p_v, expm1_v, log1pexp, log1mexp, logistic, logit
   public :: normal_pdf, normal_cdf, normal_quantile
   public :: digamma, trigamma, tetragamma, regularized_gamma_p, regularized_gamma_q
   public :: regularized_beta, beta_fn, log_beta_fn, gamma_quantile
   public :: zeta_fn, hurwitz_zeta, lambert_w0, expint_ei, expint_e1, lerch_phi
   public :: harmonic_generalized, kendall_tau

contains

   elemental real(dp) function nan_dp() result(y)
      y = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   elemental real(dp) function inf_dp() result(y)
      y = ieee_value(0.0_dp, ieee_positive_inf)
   end function inf_dp

   elemental real(dp) function log1p_v(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: u
      if (x == -1.0_dp) then
         y = -inf_dp()
      else if (x < -1.0_dp) then
         y = nan_dp()
      else if (abs(x) < 1.0e-8_dp) then
         y = x*(1.0_dp - x*(0.5_dp - x*(1.0_dp/3.0_dp - x*0.25_dp)))
      else
         u = 1.0_dp + x
         y = log(u)
      end if
   end function log1p_v

   elemental real(dp) function expm1_v(x) result(y)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-6_dp) then
         y = x*(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x*(1.0_dp/24.0_dp + x/120.0_dp))))
      else
         y = exp(x) - 1.0_dp
      end if
   end function expm1_v

   elemental real(dp) function log1pexp(x) result(y)
      real(dp), intent(in) :: x
      if (x > 35.0_dp) then
         y = x + exp(-x)
      else if (x < -35.0_dp) then
         y = exp(x)
      else if (x > 0.0_dp) then
         y = x + log1p_v(exp(-x))
      else
         y = log1p_v(exp(x))
      end if
   end function log1pexp

   elemental real(dp) function log1mexp(x) result(y)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         y = nan_dp()
      else if (x > -log(2.0_dp)) then
         y = log(-expm1_v(x))
      else
         y = log1p_v(-exp(x))
      end if
   end function log1mexp

   elemental real(dp) function logistic(x) result(p)
      real(dp), intent(in) :: x
      if (x >= 0.0_dp) then
         p = 1.0_dp/(1.0_dp + exp(-x))
      else
         p = exp(x)/(1.0_dp + exp(x))
      end if
   end function logistic

   elemental real(dp) function logit(p) result(x)
      real(dp), intent(in) :: p
      if (p <= 0.0_dp .or. p >= 1.0_dp) then
         if (p == 0.0_dp) then
            x = -inf_dp()
         else if (p == 1.0_dp) then
            x = inf_dp()
         else
            x = nan_dp()
         end if
      else
         x = log(p) - log1p_v(-p)
      end if
   end function logit

   elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x - 0.5_dp*log2pi)
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt2)
   end function normal_cdf

   elemental real(dp) function normal_quantile(p) result(x)
      ! Peter J. Acklam's rational approximation, followed by one Halley step.
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
      real(dp) :: q, r, e, u
      if (p <= 0.0_dp) then
         x = merge(-inf_dp(), nan_dp(), p == 0.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = merge(inf_dp(), nan_dp(), p == 1.0_dp)
         return
      end if
      if (p < plow) then
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
      e = normal_cdf(x) - p
      u = e / max(normal_pdf(x), tiny(1.0_dp))
      x = x - u/(1.0_dp + 0.5_dp*x*u)
   end function normal_quantile

   recursive elemental real(dp) function digamma(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, r, inv, inv2
      if (x <= 0.0_dp .and. abs(x-nint(x)) < 8.0_dp*epsilon(x)) then
         y = nan_dp(); return
      end if
      if (x < 0.0_dp) then
         y = digamma(1.0_dp-x) - pi/tan(pi*x)
         return
      end if
      z = x; r = 0.0_dp
      do while (z < 8.0_dp)
         r = r - 1.0_dp/z
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z; inv2 = inv*inv
      y = r + log(z) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - &
          inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp - inv2*5.0_dp/660.0_dp))))
   end function digamma

   recursive elemental real(dp) function trigamma(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, r, inv, inv2
      if (x <= 0.0_dp .and. abs(x-nint(x)) < 8.0_dp*epsilon(x)) then
         y = nan_dp(); return
      end if
      if (x < 0.0_dp) then
         y = pi*pi/(sin(pi*x)**2) - trigamma(1.0_dp-x)
         return
      end if
      z = x; r = 0.0_dp
      do while (z < 8.0_dp)
         r = r + 1.0_dp/(z*z)
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z; inv2 = inv*inv
      y = r + inv + 0.5_dp*inv2 + inv*inv2*(1.0_dp/6.0_dp - inv2*(1.0_dp/30.0_dp - &
          inv2*(1.0_dp/42.0_dp - inv2*(1.0_dp/30.0_dp - inv2*5.0_dp/66.0_dp))))
   end function trigamma

   elemental real(dp) function tetragamma(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: h
      if (x <= 0.0_dp) then
         y = nan_dp(); return
      end if
      h = max(1.0e-5_dp, sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x)))
      y = (trigamma(x+h)-trigamma(x-h))/(2.0_dp*h)
   end function tetragamma

   pure real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 4.0e-15_dp, fpmin = tiny(1.0_dp)/eps
      real(dp) :: ap, del, sumv, b, c, d, h, an, q
      integer :: n
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = nan_dp(); return
      else if (x == 0.0_dp) then
         p = 0.0_dp; return
      end if
      if (x < a + 1.0_dp) then
         ap = a; sumv = 1.0_dp/a; del = sumv
         do n = 1, itmax
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
         do n = 1, itmax
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d+b; if (abs(d)<fpmin) d=fpmin
            c = b+an/c; if (abs(c)<fpmin) c=fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         q = exp(-x+a*log(x)-log_gamma(a))*h
         p = max(0.0_dp,min(1.0_dp,1.0_dp-q))
      end if
   end function regularized_gamma_p

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      q = 1.0_dp - regularized_gamma_p(a,x)
   end function regularized_gamma_q

   pure real(dp) function beta_cf(a,b,x) result(h)
      real(dp), intent(in) :: a,b,x
      integer, parameter :: maxit=10000
      real(dp), parameter :: eps=4.0e-15_dp, fpmin=tiny(1.0_dp)/eps
      integer :: m,m2
      real(dp) :: aa,c,d,del,qab,qam,qap
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp; d=1.0_dp-qab*x/qap; if(abs(d)<fpmin)d=fpmin
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
   end function beta_cf

   pure real(dp) function regularized_beta(x,a,b) result(p)
      real(dp), intent(in) :: x,a,b
      real(dp) :: bt
      if (a<=0.0_dp .or. b<=0.0_dp .or. x<0.0_dp .or. x>1.0_dp) then
         p=nan_dp(); return
      else if (x==0.0_dp) then
         p=0.0_dp; return
      else if (x==1.0_dp) then
         p=1.0_dp; return
      end if
      bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log1p_v(-x))
      if(x<(a+1.0_dp)/(a+b+2.0_dp))then
         p=bt*beta_cf(a,b,x)/a
      else
         p=1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
      end if
      p=max(0.0_dp,min(1.0_dp,p))
   end function regularized_beta

   elemental real(dp) function log_beta_fn(a,b) result(y)
      real(dp), intent(in) :: a,b
      if(a<=0.0_dp .or. b<=0.0_dp)then
         y=nan_dp()
      else
         y=log_gamma(a)+log_gamma(b)-log_gamma(a+b)
      end if
   end function log_beta_fn

   elemental real(dp) function beta_fn(a,b) result(y)
      real(dp), intent(in) :: a,b
      y=exp(log_beta_fn(a,b))
   end function beta_fn

   real(dp) function gamma_quantile(p, shape, scale) result(x)
      real(dp), intent(in) :: p, shape, scale
      real(dp) :: lo,hi,mid
      integer :: iter
      if(shape<=0.0_dp .or. scale<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp)then
         x=nan_dp(); return
      else if(p==0.0_dp)then
         x=0.0_dp; return
      else if(p==1.0_dp)then
         x=inf_dp(); return
      end if
      lo=0.0_dp; hi=max(scale,shape*scale)
      do while(regularized_gamma_p(shape,hi/scale)<p)
         hi=2.0_dp*hi
         if(hi>huge(hi)/4.0_dp)exit
      end do
      do iter=1,100
         mid=0.5_dp*(lo+hi)
         if(regularized_gamma_p(shape,mid/scale)<p)then
            lo=mid
         else
            hi=mid
         end if
      end do
      x=0.5_dp*(lo+hi)
   end function gamma_quantile

   pure real(dp) function hurwitz_zeta(s,a) result(z)
      real(dp), intent(in) :: s,a
      integer, parameter :: n=20
      real(dp), parameter :: bern(8)=[1.0_dp/6.0_dp,-1.0_dp/30.0_dp,1.0_dp/42.0_dp, &
         -1.0_dp/30.0_dp,5.0_dp/66.0_dp,-691.0_dp/2730.0_dp,7.0_dp/6.0_dp,-3617.0_dp/510.0_dp]
      integer :: k,j
      real(dp) :: x,rfact,term,fact
      if(s<=1.0_dp .or. a<=0.0_dp)then
         z=nan_dp(); return
      end if
      z=0.0_dp
      do k=0,n-1
         z=z+(a+real(k,dp))**(-s)
      end do
      x=a+real(n,dp)
      z=z+x**(1.0_dp-s)/(s-1.0_dp)+0.5_dp*x**(-s)
      rfact=s; fact=2.0_dp
      do j=1,size(bern)
         if(j>1)then
            rfact=rfact*(s+real(2*j-3,dp))*(s+real(2*j-2,dp))
            fact=fact*real(2*j-1,dp)*real(2*j,dp)
         end if
         term=bern(j)/fact*rfact*x**(-s-real(2*j-1,dp))
         z=z+term
         if(abs(term)<epsilon(1.0_dp)*abs(z))exit
      end do
   end function hurwitz_zeta

   pure real(dp) function zeta_fn(s) result(z)
      real(dp), intent(in) :: s
      if(s>1.0_dp)then
         z=hurwitz_zeta(s,1.0_dp)
      else
         z=nan_dp()
      end if
   end function zeta_fn

   elemental real(dp) function lambert_w0(x) result(w)
      real(dp), intent(in) :: x
      real(dp) :: e,t,p,den
      integer :: k
      if(x < -1.0_dp/exp(1.0_dp))then
         w=nan_dp(); return
      else if(x==0.0_dp)then
         w=0.0_dp; return
      else if(abs(x)<1.0e-5_dp)then
         w=x-x*x+1.5_dp*x**3-(8.0_dp/3.0_dp)*x**4
      else if(x<1.0_dp)then
         p=sqrt(max(0.0_dp,2.0_dp*(exp(1.0_dp)*x+1.0_dp)))
         w=-1.0_dp+p-p*p/3.0_dp+11.0_dp*p**3/72.0_dp
      else
         w=log(x)-log(max(log(x),0.5_dp))
      end if
      do k=1,30
         e=exp(w); t=w*e-x
         den=e*(w+1.0_dp)-(w+2.0_dp)*t/(2.0_dp*w+2.0_dp)
         if(den==0.0_dp)exit
         p=t/den; w=w-p
         if(abs(p)<=8.0_dp*epsilon(w)*max(1.0_dp,abs(w)))exit
      end do
   end function lambert_w0

   pure real(dp) function expint_e1(x) result(v)
      real(dp), intent(in) :: x
      real(dp), parameter :: euler=0.5772156649015328606065120900824024_dp
      real(dp) :: term,sumv,del,a,b,c,d,h
      integer :: k
      if(x<0.0_dp)then
         v=nan_dp(); return
      else if(x==0.0_dp)then
         v=inf_dp(); return
      else if(x<=1.0_dp)then
         term=-x; sumv=0.0_dp
         do k=1,200
            if(k>1)term=term*(-x)*real(k-1,dp)/(real(k,dp)*real(k,dp))
            del=-term/real(k,dp)
            sumv=sumv+del
            if(abs(del)<epsilon(1.0_dp)*max(1.0_dp,abs(sumv)))exit
         end do
         v=-euler-log(x)+sumv
      else
         b=x+1.0_dp; c=1.0_dp/tiny(1.0_dp); d=1.0_dp/b; h=d
         do k=1,200
            a=-real(k*k,dp)
            b=b+2.0_dp
            d=1.0_dp/max(abs(a*d+b),tiny(1.0_dp))*sign(1.0_dp,a*d+b)
            c=b+a/c
            if(abs(c)<tiny(1.0_dp))c=tiny(1.0_dp)
            del=c*d; h=h*del
            if(abs(del-1.0_dp)<1.0e-14_dp)exit
         end do
         v=h*exp(-x)
      end if
   end function expint_e1

   pure real(dp) function expint_ei(x) result(v)
      real(dp), intent(in) :: x
      real(dp), parameter :: euler=0.5772156649015328606065120900824024_dp
      real(dp) :: term,sumv,add
      integer :: k
      if(x==0.0_dp)then
         v=-inf_dp(); return
      else if(x<0.0_dp)then
         v=-expint_e1(-x); return
      end if
      if(x<40.0_dp)then
         term=x; sumv=term
         do k=2,1000
            term=term*x/real(k,dp)
            add=term/real(k,dp)
            sumv=sumv+add
            if(abs(add)<epsilon(1.0_dp)*abs(sumv))exit
         end do
         v=euler+log(x)+sumv
      else
         term=1.0_dp; sumv=1.0_dp
         do k=1,100
            term=term*real(k,dp)/x
            if(abs(term)>abs(sumv) .and. k>x)exit
            sumv=sumv+term
            if(abs(term)<epsilon(1.0_dp)*abs(sumv))exit
         end do
         v=exp(x)*sumv/x
      end if
   end function expint_ei

   pure real(dp) function lerch_phi(z,s,a) result(v)
      real(dp), intent(in) :: z,s,a
      real(dp) :: term,sumv
      integer :: k
      if(abs(z)>=1.0_dp .or. a<=0.0_dp)then
         v=nan_dp(); return
      end if
      sumv=0.0_dp; term=1.0_dp
      do k=0,100000
         if(k>0)term=term*z
         sumv=sumv+term/(a+real(k,dp))**s
         if(abs(term/(a+real(k,dp))**s)<1.0e-14_dp*max(1.0_dp,abs(sumv)))exit
      end do
      v=sumv
   end function lerch_phi

   pure real(dp) function harmonic_generalized(n,s) result(h)
      integer, intent(in) :: n
      real(dp), intent(in) :: s
      integer :: k
      h=0.0_dp
      do k=1,max(0,n)
         h=h+real(k,dp)**(-s)
      end do
   end function harmonic_generalized

   pure real(dp) function kendall_tau(x,y) result(tau)
      real(dp), intent(in) :: x(:),y(:)
      integer :: i,j,n
      real(dp) :: c,d
      n=min(size(x),size(y)); c=0.0_dp; d=0.0_dp
      do i=1,n-1
         do j=i+1,n
            if((x(i)-x(j))*(y(i)-y(j))>0.0_dp)c=c+1.0_dp
            if((x(i)-x(j))*(y(i)-y(j))<0.0_dp)d=d+1.0_dp
         end do
      end do
      if(c+d>0.0_dp)then
         tau=(c-d)/(c+d)
      else
         tau=0.0_dp
      end if
   end function kendall_tau

end module gamlss_special
