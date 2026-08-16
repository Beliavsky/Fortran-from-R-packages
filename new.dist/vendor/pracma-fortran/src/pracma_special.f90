! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use pracma_kinds, only : dp, pi_dp, eps_dp, sqrt_two_pi_dp
   use pracma_status
   implicit none
   private

   public :: sinc, psinc, agmean, gammaz, gammainc, incgam, psi
   public :: erf_pracma, erfc_pracma, erfinv, erfcinv, erfcx, erfz, erfi
   public :: expint, expint_Ei, li, Si, Ci, fresnel
   public :: lambertWp, lambertWn, zeta, eta, polylog
   public :: ellipke, ellipj, nthroot, einsteinF, bernoulli
   public :: factorial2, golden_ratio, humps

contains

   elemental real(dp) function sinc(x)
      real(dp),intent(in)::x
      if(abs(x)<=sqrt(eps_dp))then
         sinc=1.0_dp-x*x/6.0_dp+x**4/120.0_dp
      else
         sinc=sin(x)/x
      end if
   end function sinc

   elemental real(dp) function psinc(x)
      real(dp),intent(in)::x
      psinc=sinc(pi_dp*x)
   end function psinc

   subroutine agmean(a,b,agm,iterations,tolerance,status)
      real(dp),intent(in)::a,b
      real(dp),intent(out)::agm
      integer,intent(out),optional::iterations,status
      real(dp),intent(in),optional::tolerance
      real(dp)::x,y,xn,yn,tol
      integer::iter,istat
      x=a; y=b; tol=100.0_dp*eps_dp; istat=pracma_ok
      if(present(tolerance))tol=tolerance
      if(a<0.0_dp .or. b<0.0_dp)then
         agm=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
         if(present(iterations))iterations=0
         if(present(status))status=istat
         return
      end if
      do iter=1,1000
         xn=0.5_dp*(x+y); yn=sqrt(x*y)
         if(abs(xn-yn)<=tol*max(1.0_dp,abs(xn)))exit
         x=xn; y=yn
      end do
      agm=0.5_dp*(xn+yn)
      if(iter>1000)istat=pracma_not_converged
      if(present(iterations))iterations=min(iter,1000)
      if(present(status))status=istat
   end subroutine agmean

   recursive function gammaz(z) result(g)
      complex(dp),intent(in)::z
      complex(dp)::g,t,s
      real(dp),parameter::coef(9)=[0.99999999999980993_dp,676.5203681218851_dp, &
         -1259.1392167224028_dp,771.32342877765313_dp,-176.61502916214059_dp, &
         12.507343278686905_dp,-0.13857109526572012_dp,9.9843695780195716e-6_dp, &
         1.5056327351493116e-7_dp]
      integer::i
      if(real(z,dp)<0.5_dp)then
         g=pi_dp/(sin(pi_dp*z)*gammaz(1.0_dp-z))
      else
         t=z-1.0_dp; s=cmplx(coef(1),0.0_dp,dp)
         do i=2,size(coef)
            s=s+coef(i)/(t+real(i-1,dp))
         end do
         g=sqrt_two_pi_dp*(t+7.5_dp)**(t+0.5_dp)*exp(-(t+7.5_dp))*s
      end if
   end function gammaz

   function gammainc(x,a,upper,regularized,status) result(v)
      real(dp),intent(in)::x,a
      logical,intent(in),optional::upper,regularized
      integer,intent(out),optional::status
      real(dp)::v,p,q
      logical::up,reg
      integer::istat
      up=.false.; reg=.true.
      if(present(upper))up=upper
      if(present(regularized))reg=regularized
      call regularized_gamma(a,x,p,q,istat)
      if(up)then; v=q; else; v=p; end if
      if(.not.reg .and. istat==pracma_ok)v=v*gamma(a)
      if(present(status))status=istat
   end function gammainc

   function incgam(a,x,status) result(v)
      real(dp),intent(in)::a,x
      integer,intent(out),optional::status
      real(dp)::v
      integer::istat
      v=gammainc(x,a,.false.,.false.,istat)
      if(present(status))status=istat
   end function incgam

   subroutine regularized_gamma(a,x,p,q,status)
      real(dp),intent(in)::a,x
      real(dp),intent(out)::p,q
      integer,intent(out)::status
      real(dp)::sumv,del,ap,b,c,d,h,an,logpref
      integer::n
      if(a<=0.0_dp .or. x<0.0_dp)then
         p=ieee_value(0.0_dp,ieee_quiet_nan); q=p; status=pracma_invalid_argument; return
      end if
      if(x<=tiny(1.0_dp))then
         p=0.0_dp; q=1.0_dp; status=pracma_ok; return
      end if
      logpref=a*log(x)-x-log_gamma(a)
      if(x<a+1.0_dp)then
         ap=a; sumv=1.0_dp/a; del=sumv
         do n=1,10000
            ap=ap+1.0_dp; del=del*x/ap; sumv=sumv+del
            if(abs(del)<=abs(sumv)*1.0e-15_dp)exit
         end do
         p=sumv*exp(logpref); q=max(0.0_dp,1.0_dp-p)
      else
         b=x+1.0_dp-a; c=1.0_dp/tiny(1.0_dp); d=1.0_dp/b; h=d
         do n=1,10000
            an=-real(n,dp)*(real(n,dp)-a)
            b=b+2.0_dp; d=an*d+b
            if(abs(d)<tiny(1.0_dp))d=tiny(1.0_dp)
            c=b+an/c
            if(abs(c)<tiny(1.0_dp))c=tiny(1.0_dp)
            d=1.0_dp/d; del=d*c; h=h*del
            if(abs(del-1.0_dp)<=1.0e-15_dp)exit
         end do
         q=exp(logpref)*h; p=max(0.0_dp,1.0_dp-q)
      end if
      status=merge(pracma_ok,pracma_not_converged,n<10000)
   end subroutine regularized_gamma

   recursive function psi(x) result(v)
      real(dp),intent(in)::x
      real(dp)::v,y,r,inv,inv2
      if(x<=0.0_dp .and. abs(x-real(nint(x),dp))<=eps_dp)then
         v=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      if(x<0.0_dp)then
         v=psi(1.0_dp-x)-pi_dp/tan(pi_dp*x)
         return
      end if
      y=x; r=0.0_dp
      do while(y<8.0_dp)
         r=r-1.0_dp/y; y=y+1.0_dp
      end do
      inv=1.0_dp/y; inv2=inv*inv
      v=r+log(y)-0.5_dp*inv-inv2*(1.0_dp/12.0_dp-inv2*(1.0_dp/120.0_dp- &
          inv2*(1.0_dp/252.0_dp-inv2/240.0_dp)))
   end function psi

   elemental real(dp) function erf_pracma(x)
      real(dp),intent(in)::x
      erf_pracma=erf(x)
   end function erf_pracma

   elemental real(dp) function erfc_pracma(x)
      real(dp),intent(in)::x
      erfc_pracma=erfc(x)
   end function erfc_pracma

   elemental real(dp) function erfinv(x)
      real(dp),intent(in)::x
      real(dp)::w,p,y,err
      integer::i
      if(abs(x)>1.0_dp)then
         erfinv=ieee_value(0.0_dp,ieee_quiet_nan); return
      else if(abs(x-1.0_dp)<=eps_dp)then
         erfinv=huge(1.0_dp); return
      else if(abs(x+1.0_dp)<=eps_dp)then
         erfinv=-huge(1.0_dp); return
      end if
      w=-log((1.0_dp-x)*(1.0_dp+x))
      if(w<5.0_dp)then
         w=w-2.5_dp
         p=2.81022636e-08_dp
         p=3.43273939e-07_dp+p*w; p=-3.5233877e-06_dp+p*w
         p=-4.39150654e-06_dp+p*w; p=0.00021858087_dp+p*w
         p=-0.00125372503_dp+p*w; p=-0.00417768164_dp+p*w
         p=0.246640727_dp+p*w; p=1.50140941_dp+p*w
      else
         w=sqrt(w)-3.0_dp
         p=-0.000200214257_dp
         p=0.000100950558_dp+p*w; p=0.00134934322_dp+p*w
         p=-0.00367342844_dp+p*w; p=0.00573950773_dp+p*w
         p=-0.0076224613_dp+p*w; p=0.00943887047_dp+p*w
         p=1.00167406_dp+p*w; p=2.83297682_dp+p*w
      end if
      y=p*x
      do i=1,3
         err=erf(y)-x
         y=y-err/(2.0_dp/sqrt(pi_dp)*exp(-y*y))
      end do
      erfinv=y
   end function erfinv

   elemental real(dp) function erfcinv(x)
      real(dp),intent(in)::x
      erfcinv=erfinv(1.0_dp-x)
   end function erfcinv

   elemental real(dp) function erfcx(x)
      real(dp),intent(in)::x
      if(x<26.0_dp)then
         erfcx=exp(x*x)*erfc(x)
      else
         erfcx=(1.0_dp/sqrt(pi_dp))/x*(1.0_dp+1.0_dp/(2.0_dp*x*x)+3.0_dp/(4.0_dp*x**4))
      end if
   end function erfcx

   function erfz(z) result(v)
      complex(dp),intent(in)::z
      complex(dp)::v,term,sumv
      integer::n
      term=z; sumv=z
      do n=1,200
         term=term*(-z*z)/real(n,dp)
         sumv=sumv+term/real(2*n+1,dp)
         if(abs(term/real(2*n+1,dp))<=eps_dp*max(1.0_dp,abs(sumv)))exit
      end do
      v=2.0_dp/sqrt(pi_dp)*sumv
   end function erfz

   function erfi(z) result(v)
      complex(dp),intent(in)::z
      complex(dp)::v
      v=cmplx(0.0_dp,-1.0_dp,dp)*erfz(cmplx(0.0_dp,1.0_dp,dp)*z)
   end function erfi

   function expint(x,status) result(v)
      real(dp),intent(in)::x
      integer,intent(out),optional::status
      real(dp)::v,term,sumv
      integer::k,istat
      if(x<0.0_dp)then
         v=-expint_Ei(-x); istat=pracma_ok
      else if(x<=tiny(1.0_dp))then
         v=huge(1.0_dp); istat=pracma_invalid_argument
      else if(x<=1.0_dp)then
         term=-x; sumv=-0.5772156649015328606_dp-log(x)-term
         do k=2,500
            term=term*(-x)*real(k-1,dp)/(real(k,dp)*real(k,dp))
            sumv=sumv-term
            if(abs(term)<=eps_dp*max(1.0_dp,abs(sumv)))exit
         end do
         v=sumv; istat=pracma_ok
      else
         term=1.0_dp; sumv=1.0_dp
         do k=1,200
            term=-term*real(k,dp)/x
            if(abs(term)>abs(sumv) .and. k>2)exit
            sumv=sumv+term
         end do
         v=exp(-x)/x*sumv; istat=pracma_ok
      end if
      if(present(status))status=istat
   end function expint

   function expint_Ei(x) result(v)
      real(dp),intent(in)::x
      real(dp)::v,term,sumv
      integer::k
      if(abs(x)<=tiny(1.0_dp))then
         v=-huge(1.0_dp); return
      end if
      if(x<0.0_dp)then
         v=-expint(-x); return
      end if
      if(x<=40.0_dp)then
         term=x; sumv=x
         do k=2,1000
            term=term*x/real(k,dp)
            sumv=sumv+term/real(k,dp)
            if(abs(term/real(k,dp))<=eps_dp*max(1.0_dp,abs(sumv)))exit
         end do
         v=0.5772156649015328606_dp+log(x)+sumv
      else
         term=1.0_dp; sumv=1.0_dp
         do k=1,200
            term=term*real(k,dp)/x
            if(abs(term)>abs(sumv) .and. k>2)exit
            sumv=sumv+term
         end do
         v=exp(x)/x*sumv
      end if
   end function expint_Ei

   function li(x) result(v)
      real(dp),intent(in)::x
      real(dp)::v
      if(x<=0.0_dp .or. abs(x-1.0_dp)<=eps_dp)then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         v=expint_Ei(log(x))
      end if
   end function li

   function Si(x) result(v)
      real(dp),intent(in)::x
      real(dp)::v,term,sumv
      integer::k
      term=x; sumv=x
      do k=1,500
         term=term*(-x*x)/real((2*k)*(2*k+1),dp)
         sumv=sumv+term/real(2*k+1,dp)
         if(abs(term/real(2*k+1,dp))<=eps_dp*max(1.0_dp,abs(sumv)))exit
      end do
      v=sumv
   end function Si

   function Ci(x) result(v)
      real(dp),intent(in)::x
      real(dp)::v,term,sumv
      integer::k
      if(x<=0.0_dp)then
         v=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      term=-x*x/2.0_dp; sumv=term/2.0_dp
      do k=2,500
         term=term*(-x*x)/real((2*k-1)*(2*k),dp)
         sumv=sumv+term/real(2*k,dp)
         if(abs(term/real(2*k,dp))<=eps_dp*max(1.0_dp,abs(sumv)))exit
      end do
      v=0.5772156649015328606_dp+log(x)+sumv
   end function Ci

   subroutine fresnel(x,c,s)
      real(dp),intent(in)::x
      real(dp),intent(out)::c,s
      real(dp)::termc,terms,sumc,sums,z
      integer::k
      z=0.5_dp*pi_dp*x*x
      termc=x; sumc=x; terms=pi_dp*x**3/6.0_dp; sums=terms
      do k=1,200
         termc=termc*(-z*z)*real(4*k-3,dp)/(real(2*k,dp)*real(2*k-1,dp)*real(4*k+1,dp))
         terms=terms*(-z*z)*real(4*k-1,dp)/(real(2*k+1,dp)*real(2*k,dp)*real(4*k+3,dp))
         sumc=sumc+termc; sums=sums+terms
         if(max(abs(termc),abs(terms))<=eps_dp*max(1.0_dp,max(abs(sumc),abs(sums))))exit
      end do
      c=sumc; s=sums
   end subroutine fresnel

   function lambertWp(x,status) result(w)
      real(dp),intent(in)::x
      integer,intent(out),optional::status
      real(dp)::w,ew,f,den,d
      integer::i,istat
      if(x<-1.0_dp/exp(1.0_dp))then
         w=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      if(x<1.0_dp)then
         w=x
         if(x<-0.2_dp)w=-1.0_dp+sqrt(2.0_dp*(1.0_dp+exp(1.0_dp)*x))
      else
         w=log(x)-log(log(x))
      end if
      istat=pracma_not_converged
      do i=1,100
         ew=exp(w); f=w*ew-x
         den=ew*(w+1.0_dp)-(w+2.0_dp)*f/(2.0_dp*w+2.0_dp)
         d=f/den; w=w-d
         if(abs(d)<=1.0e-14_dp*(1.0_dp+abs(w)))then
            istat=pracma_ok; exit
         end if
      end do
      if(present(status))status=istat
   end function lambertWp

   function lambertWn(x,status) result(w)
      real(dp),intent(in)::x
      integer,intent(out),optional::status
      real(dp)::w,ew,f,den,d
      integer::i,istat
      if(x<-1.0_dp/exp(1.0_dp) .or. x>=0.0_dp)then
         w=ieee_value(0.0_dp,ieee_quiet_nan); istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      w=log(-x)-log(-log(-x))
      if(x<-0.2_dp)w=-1.0_dp-sqrt(2.0_dp*(1.0_dp+exp(1.0_dp)*x))
      istat=pracma_not_converged
      do i=1,100
         ew=exp(w); f=w*ew-x
         den=ew*(w+1.0_dp)-(w+2.0_dp)*f/(2.0_dp*w+2.0_dp)
         d=f/den; w=w-d
         if(abs(d)<=1.0e-14_dp*(1.0_dp+abs(w)))then
            istat=pracma_ok; exit
         end if
      end do
      if(present(status))status=istat
   end function lambertWn

   function eta(s) result(v)
      real(dp),intent(in)::s
      real(dp)::v,term,sumv
      integer::k
      sumv=0.0_dp
      do k=1,1000000
         term=(-1.0_dp)**(k-1)/real(k,dp)**s
         sumv=sumv+term
         if(abs(term)<=1.0e-14_dp*max(1.0_dp,abs(sumv)))exit
      end do
      v=sumv
   end function eta

   recursive function zeta(s) result(v)
      real(dp),intent(in)::s
      real(dp)::v
      if(abs(s-1.0_dp)<=eps_dp)then
         v=huge(1.0_dp)
      else if(s>0.0_dp)then
         v=eta(s)/(1.0_dp-2.0_dp**(1.0_dp-s))
      else
         v=2.0_dp**s*pi_dp**(s-1.0_dp)*sin(0.5_dp*pi_dp*s)*gamma(1.0_dp-s)*zeta(1.0_dp-s)
      end if
   end function zeta

   function polylog(s,z) result(v)
      real(dp),intent(in)::s
      complex(dp),intent(in)::z
      complex(dp)::v,term,sumv
      integer::k
      if(abs(z)>=1.0_dp)then
         v=cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp); return
      end if
      term=z; sumv=z
      do k=2,1000000
         term=term*z
         sumv=sumv+term/real(k,dp)**s
         if(abs(term)/real(k,dp)**s<=1.0e-14_dp*max(1.0_dp,abs(sumv)))exit
      end do
      v=sumv
   end function polylog

   subroutine ellipke(m,k,e,status)
      real(dp),intent(in)::m
      real(dp),intent(out)::k,e
      integer,intent(out),optional::status
      real(dp)::a,b,c,sumc,pow2
      integer::iter,istat
      if(m<0.0_dp .or. m>1.0_dp)then
         k=ieee_value(0.0_dp,ieee_quiet_nan); e=k; istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      if(abs(m-1.0_dp)<=eps_dp)then
         k=huge(1.0_dp); e=1.0_dp; istat=pracma_ok
         if(present(status))status=istat
         return
      end if
      a=1.0_dp; b=sqrt(1.0_dp-m); sumc=0.0_dp; pow2=1.0_dp
      do iter=1,100
         c=0.5_dp*(a-b)
         if(abs(c)<=eps_dp*a)exit
         sumc=sumc+pow2*c*c; pow2=2.0_dp*pow2
         b=sqrt(a*b); a=a-c
      end do
      k=pi_dp/(2.0_dp*a)
      e=k*(1.0_dp-sumc)
      istat=pracma_ok
      if(present(status))status=istat
   end subroutine ellipke

   subroutine ellipj(u,m,sn,cn,dn,status)
      real(dp),intent(in)::u,m
      real(dp),intent(out)::sn,cn,dn
      integer,intent(out),optional::status
      real(dp)::a(40),c(40),b,phi,t
      integer::n,j,istat
      if(m<0.0_dp .or. m>1.0_dp)then
         sn=ieee_value(0.0_dp,ieee_quiet_nan); cn=sn; dn=sn
         istat=pracma_invalid_argument
         if(present(status))status=istat
         return
      end if
      a(1)=1.0_dp; b=sqrt(1.0_dp-m); c(1)=sqrt(m); n=1
      do while(abs(c(n))>eps_dp*a(n) .and. n<39)
         n=n+1; a(n)=0.5_dp*(a(n-1)+b); c(n)=0.5_dp*(a(n-1)-b); b=sqrt(a(n-1)*b)
      end do
      phi=2.0_dp**(n-1)*a(n)*u
      do j=n,2,-1
         t=c(j)*sin(phi)/a(j)
         phi=0.5_dp*(asin(max(-1.0_dp,min(1.0_dp,t)))+phi)
      end do
      sn=sin(phi); cn=cos(phi); dn=sqrt(max(0.0_dp,1.0_dp-m*sn*sn)); istat=pracma_ok
      if(present(status))status=istat
   end subroutine ellipj

   elemental real(dp) function nthroot(x,n)
      real(dp),intent(in)::x
      integer,intent(in)::n
      if(n==0 .or. (x<0.0_dp .and. modulo(n,2)==0))then
         nthroot=ieee_value(0.0_dp,ieee_quiet_nan)
      else
         nthroot=sign(abs(x)**(1.0_dp/real(n,dp)),x)
      end if
   end function nthroot

   elemental real(dp) function einsteinF(x)
      real(dp),intent(in)::x
      if(abs(x)<=sqrt(eps_dp))then
         einsteinF=1.0_dp+x/2.0_dp+x*x/12.0_dp-x**4/720.0_dp
      else
         einsteinF=x/(1.0_dp-exp(-x))
      end if
   end function einsteinF

   function bernoulli(n) result(b)
      integer,intent(in)::n
      real(dp),allocatable::b(:)
      real(dp),allocatable::a(:)
      integer::m,j
      allocate(b(n+1),a(n+1)); b=0.0_dp
      do m=0,n
         a(m+1)=1.0_dp/real(m+1,dp)
         do j=m,1,-1
            a(j)=real(j,dp)*(a(j)-a(j+1))
         end do
         b(m+1)=a(1)
      end do
      if(n>=1)b(2)=-0.5_dp
   end function bernoulli

   pure real(dp) function factorial2(n)
      integer,intent(in)::n
      integer::k
      if(n<-1)then
         factorial2=ieee_value(0.0_dp,ieee_quiet_nan); return
      end if
      factorial2=1.0_dp
      do k=max(1,n),2,-2
         factorial2=factorial2*real(k,dp)
      end do
   end function factorial2

   pure real(dp) function golden_ratio()
      golden_ratio=0.5_dp*(1.0_dp+sqrt(5.0_dp))
   end function golden_ratio

   elemental real(dp) function humps(x)
      real(dp),intent(in)::x
      humps=1.0_dp/((x-0.3_dp)**2+0.01_dp)+1.0_dp/((x-0.9_dp)**2+0.04_dp)-6.0_dp
   end function humps

end module pracma_special
