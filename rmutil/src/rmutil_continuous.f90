! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use rmutil_kinds, only : dp, pi
   use rmutil_special, only : normal_pdf, normal_cdf, normal_quantile, &
      regularized_gamma_p, gamma_quantile
   use rmutil_integrate, only : integrate_adaptive
   implicit none
   private
   public :: pinvgauss, dinvgauss, qinvgauss, rinvgauss
   public :: plaplace, dlaplace, qlaplace, rlaplace
   public :: plevy, dlevy, qlevy, rlevy
   public :: ppareto, dpareto, qpareto, rpareto
   public :: psimplex, dsimplex, qsimplex, rsimplex
   public :: ptwosidedpower, dtwosidedpower, qtwosidedpower, rtwosidedpower
   public :: pboxcox, dboxcox, qboxcox, rboxcox
   public :: pburr, dburr, qburr, rburr
   public :: pgextval, dgextval, qgextval, rgextval
   public :: pggamma, dggamma, qggamma, rggamma
   public :: pginvgauss, dginvgauss, qginvgauss, rginvgauss
   public :: pglogis, dglogis, qglogis, rglogis
   public :: pgweibull, dgweibull, qgweibull, rgweibull
   public :: phjorth, dhjorth, qhjorth, rhjorth
   public :: ppowexp, dpowexp, qpowexp, rpowexp
   public :: pskewlaplace, dskewlaplace, qskewlaplace, rskewlaplace

   abstract interface
      function local_cdf(x) result(p)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: p
      end function local_cdf
   end interface

contains

   elemental real(dp) function pinvgauss(q,m,s) result(p)
      real(dp), intent(in) :: q,m,s
      real(dp) :: t,v
      if (q <= 0.0_dp .or. m <= 0.0_dp .or. s <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      t = q/m
      v = sqrt(q*s)
      p = normal_cdf((t-1.0_dp)/v) + exp(2.0_dp/(m*s))*normal_cdf(-(t+1.0_dp)/v)
      p = min(1.0_dp,max(0.0_dp,p))
   end function pinvgauss

   elemental real(dp) function dinvgauss(y,m,s,log_p) result(v)
      real(dp), intent(in) :: y,m,s
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg = .false.; if (present(log_p)) lg = log_p
      if (y <= 0.0_dp .or. m <= 0.0_dp .or. s <= 0.0_dp) then
         v = merge(-huge(1.0_dp),0.0_dp,lg)
         return
      end if
      lv = -(y-m)**2/(2.0_dp*y*s*m*m) - &
         (log(2.0_dp*pi*s)+3.0_dp*log(y))/2.0_dp
      v = merge(lv,exp(lv),lg)
   end function dinvgauss

   real(dp) function qinvgauss(p,m,s) result(x)
      real(dp), intent(in) :: p,m,s
      real(dp) :: lo, hi
      if (p <= 0.0_dp) then
         x = 0.0_dp; return
      else if (p >= 1.0_dp) then
         x = ieee_value(0.0_dp,ieee_positive_inf); return
      end if
      lo = tiny(1.0_dp)
      hi = max(20.0_dp,m+10.0_dp*sqrt(max(s*m**3,tiny(1.0_dp))))
      do while (pinvgauss(hi,m,s) < p)
         hi = 2.0_dp*hi
      end do
      x = bisect(p,cdf,lo,hi)
   contains
      function cdf(z) result(v)
         real(dp), intent(in) :: z
         real(dp) :: v
         v = pinvgauss(z,m,s)
      end function cdf
   end function qinvgauss

   function rinvgauss(n,m,s) result(x)
      integer, intent(in) :: n
      real(dp), intent(in) :: m,s
      real(dp), allocatable :: x(:)
      real(dp) :: u
      integer :: i
      allocate(x(n))
      do i=1,n; call random_number(u); x(i)=qinvgauss(u,m,s); end do
   end function rinvgauss

   elemental real(dp) function plaplace(q,m,s) result(p)
      real(dp), intent(in) :: q,m,s
      real(dp) :: u,t
      if (s <= 0.0_dp) then; p=0.0_dp; return; end if
      u=(q-m)/s; t=exp(-abs(u))/2.0_dp
      p=merge(t,1.0_dp-t,u<0.0_dp)
   end function plaplace

   elemental real(dp) function dlaplace(y,m,s,log_p) result(v)
      real(dp), intent(in) :: y,m,s
      logical, intent(in), optional :: log_p
      real(dp) :: lv
      logical :: lg
      lg=.false.; if(present(log_p))lg=log_p
      lv=-abs(y-m)/s-log(2.0_dp*s)
      v=merge(lv,exp(lv),lg)
   end function dlaplace

   elemental real(dp) function qlaplace(p,m,s) result(x)
      real(dp), intent(in) :: p,m,s
      if(p<0.5_dp)then
         x=s*log(2.0_dp*p)+m
      else
         x=-s*log(2.0_dp*(1.0_dp-p))+m
      end if
   end function qlaplace

   function rlaplace(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qlaplace(u,m,s); end do
   end function rlaplace

   elemental real(dp) function plevy(q,m,s) result(p)
      real(dp),intent(in)::q,m,s
      if(q<=m .or. s<=0.0_dp)then; p=0.0_dp; else
         p=2.0_dp*(1.0_dp-normal_cdf(1.0_dp/sqrt((q-m)/s)))
      end if
   end function plevy

   elemental real(dp) function dlevy(y,m,s,log_p) result(v)
      real(dp),intent(in)::y,m,s; logical,intent(in),optional::log_p
      logical::lg; real(dp)::lv
      lg=.false.; if(present(log_p))lg=log_p
      if(y<=m .or. s<=0.0_dp)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      lv=0.5_dp*log(s/(2.0_dp*pi))-1.5_dp*log(y-m)-s/(2.0_dp*(y-m))
      v=merge(lv,exp(lv),lg)
   end function dlevy

   elemental real(dp) function qlevy(p,m,s) result(x)
      real(dp),intent(in)::p,m,s
      x=s/normal_quantile(1.0_dp-p/2.0_dp)**2+m
   end function qlevy

   function rlevy(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qlevy(u,m,s); end do
   end function rlevy

   elemental real(dp) function ppareto(q,m,s) result(p)
      real(dp),intent(in)::q,m,s
      if(q<0.0_dp)then; p=0.0_dp; else; p=1.0_dp-(1.0_dp+q/(m*(s-1.0_dp)))**(-s); end if
   end function ppareto

   elemental real(dp) function dpareto(y,m,s,log_p) result(v)
      real(dp),intent(in)::y,m,s; logical,intent(in),optional::log_p
      real(dp)::mm,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0.0_dp)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      mm=m*(s-1.0_dp); lv=log(s)-(s+1.0_dp)*log(1.0_dp+y/mm)-log(mm)
      v=merge(lv,exp(lv),lg)
   end function dpareto

   elemental real(dp) function qpareto(p,m,s) result(x)
      real(dp),intent(in)::p,m,s
      x=((1.0_dp-p)**(-1.0_dp/s)-1.0_dp)*m*(s-1.0_dp)
   end function qpareto

   function rpareto(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qpareto(u,m,s); end do
   end function rpareto

   elemental real(dp) function dsimplex(y,m,s,log_p) result(v)
      real(dp),intent(in)::y,m,s; logical,intent(in),optional::log_p
      real(dp)::lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<=0.0_dp .or. y>=1.0_dp)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      lv=-((y-m)/(m*(1.0_dp-m)))**2/(2.0_dp*y*(1.0_dp-y)*s) - &
         (log(2.0_dp*pi*s)+3.0_dp*(log(y)+log(1.0_dp-y)))/2.0_dp
      v=merge(lv,exp(lv),lg)
   end function dsimplex

   real(dp) function psimplex(q,m,s) result(p)
      real(dp),intent(in)::q,m,s
      if(q<=0.0_dp)then; p=0.0_dp; return; else if(q>=1.0_dp)then; p=1.0_dp; return; end if
      p=integrate_adaptive(integrand,0.0_dp,1.0_dp,1.0e-9_dp)
      p=min(1.0_dp,max(0.0_dp,p))
   contains
      function integrand(u) result(v)
         real(dp),intent(in)::u; real(dp)::v,y
         if(u<=0.0_dp)then; v=0.0_dp; return; end if
         y=q*u
         v=q*dsimplex(y,m,s)
      end function integrand
   end function psimplex

   real(dp) function qsimplex(p,m,s) result(x)
      real(dp),intent(in)::p,m,s
      x=bisect(p,cdf,1.0e-12_dp,1.0_dp-1.0e-12_dp)
   contains
      function cdf(z) result(v); real(dp),intent(in)::z; real(dp)::v; v=psimplex(z,m,s); end function cdf
   end function qsimplex

   function rsimplex(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qsimplex(u,m,s); end do
   end function rsimplex

   elemental real(dp) function ptwosidedpower(q,m,s) result(p)
      real(dp),intent(in)::q,m,s
      if(q<m)then; p=m*(q/m)**s; else; p=1.0_dp-(1.0_dp-m)*((1.0_dp-q)/(1.0_dp-m))**s; end if
   end function ptwosidedpower

   elemental real(dp) function dtwosidedpower(y,m,s,log_p) result(v)
      real(dp),intent(in)::y,m,s; logical,intent(in),optional::log_p
      real(dp)::lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      lv=log(s)+(s-1.0_dp)*merge(log(y)-log(m),log(1.0_dp-y)-log(1.0_dp-m),y<m)
      v=merge(lv,exp(lv),lg)
   end function dtwosidedpower

   elemental real(dp) function qtwosidedpower(p,m,s) result(x)
      real(dp),intent(in)::p,m,s
      if(p<m)then; x=m*(p/m)**(1.0_dp/s); else
         x=1.0_dp-(1.0_dp-m)*((1.0_dp-p)/(1.0_dp-m))**(1.0_dp/s)
      end if
   end function qtwosidedpower

   function rtwosidedpower(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qtwosidedpower(u,m,s); end do
   end function rtwosidedpower

   elemental real(dp) function pboxcox(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      real(dp)::normv,den,num
      normv=sign(1.0_dp,f)*normal_cdf((0.0_dp-m)/sqrt(s))
      num=normal_cdf((q**f/f-m)/sqrt(s))-merge(normv,0.0_dp,f>0.0_dp)
      den=1.0_dp-merge(1.0_dp,0.0_dp,f<0.0_dp)-normv
      p=num/den
   end function pboxcox

   elemental real(dp) function dboxcox(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::normv,den,z,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      normv=sign(1.0_dp,f)*normal_cdf(-m/sqrt(s))
      den=1.0_dp-merge(1.0_dp,0.0_dp,f<0.0_dp)-normv
      z=(y**f/f-m)/sqrt(s)
      lv=(f-1.0_dp)*log(y)+log(normal_pdf(z))-0.5_dp*log(s)-log(den)
      v=merge(lv,exp(lv),lg)
   end function dboxcox

   real(dp) function qboxcox(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      real(dp)::hi
      hi=max(20.0_dp,m+10.0_dp*sqrt(s)+1.0_dp)
      do while(pboxcox(hi,m,s,f)<p); hi=2.0_dp*hi; end do
      x=bisect(p,cdf,tiny(1.0_dp),hi)
   contains
      function cdf(z) result(v); real(dp),intent(in)::z; real(dp)::v; v=pboxcox(z,m,s,f); end function cdf
   end function qboxcox

   function rboxcox(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qboxcox(u,m,s,f); end do
   end function rboxcox

   elemental real(dp) function pburr(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      p=1.0_dp-(1.0_dp+(q/m)**s)**(-f)
   end function pburr

   elemental real(dp) function dburr(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::y1,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      y1=y/m; lv=log(f*s)+(s-1.0_dp)*log(y1)-log(m)-(f+1.0_dp)*log(1.0_dp+y1**s)
      v=merge(lv,exp(lv),lg)
   end function dburr

   elemental real(dp) function qburr(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      x=((1.0_dp-p)**(-1.0_dp/f)-1.0_dp)**(1.0_dp/s)*m
   end function qburr

   function rburr(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qburr(u,m,s,f); end do
   end function rburr

   elemental real(dp) function weibull_cdf(x,shape,scale) result(p)
      real(dp),intent(in)::x,shape,scale
      if(x<=0.0_dp)then; p=0.0_dp; else; p=1.0_dp-exp(-(x/scale)**shape); end if
   end function weibull_cdf

   elemental real(dp) function weibull_quantile(p,shape,scale) result(x)
      real(dp),intent(in)::p,shape,scale
      x=scale*(-log(1.0_dp-p))**(1.0_dp/shape)
   end function weibull_quantile

   elemental real(dp) function pgextval(q,s,m,f) result(p)
      real(dp),intent(in)::q,s,m,f
      real(dp)::normv; logical::ind
      normv=sign(1.0_dp,f)*exp(-m**(-s)); ind=f>0.0_dp
      p=(weibull_cdf(exp(q**f/f),s,m)-merge(1.0_dp,0.0_dp,ind)+ &
         merge(normv,0.0_dp,ind))/(1.0_dp-merge(1.0_dp,0.0_dp,ind)+normv)
   end function pgextval

   elemental real(dp) function dgextval(y,s,m,f,log_p) result(v)
      real(dp),intent(in)::y,s,m,f; logical,intent(in),optional::log_p
      real(dp)::normv,y1,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      normv=sign(1.0_dp,f)*exp(-m**(-s)); y1=exp(y**f/f)
      lv=(f-1.0_dp)*log(y)+log(y1)+log(s)+(s-1.0_dp)*log(y1)-s*log(m) - &
         (y1/m)**s-log(1.0_dp-merge(1.0_dp,0.0_dp,f>0.0_dp)+normv)
      v=merge(lv,exp(lv),lg)
   end function dgextval

   elemental real(dp) function qgextval(p,s,m,f) result(x)
      real(dp),intent(in)::p,s,m,f
      real(dp)::normv,target; logical::ind
      normv=sign(1.0_dp,f)*exp(-m**(-s)); ind=f>0.0_dp
      target=p*(1.0_dp-merge(1.0_dp,0.0_dp,ind)+normv)+ &
         merge(1.0_dp,0.0_dp,ind)-merge(normv,0.0_dp,ind)
      x=(f*log(weibull_quantile(target,s,m)))**(1.0_dp/f)
   end function qgextval

   function rgextval(n,s,m,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::s,m,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qgextval(u,s,m,f); end do
   end function rgextval

   elemental real(dp) function pggamma(q,s,m,f) result(p)
      real(dp),intent(in)::q,s,m,f
      p=regularized_gamma_p(s,q**f/(m/s)**f)
   end function pggamma

   elemental real(dp) function dggamma(y,s,m,f,log_p) result(v)
      real(dp),intent(in)::y,s,m,f; logical,intent(in),optional::log_p
      real(dp)::z,scale,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      z=y**f; scale=(m/s)**f
      lv=log(f)+(f-1.0_dp)*log(y)+(s-1.0_dp)*log(z)-z/scale- &
         log_gamma(s)-s*log(scale)
      v=merge(lv,exp(lv),lg)
   end function dggamma

   elemental real(dp) function qggamma(p,s,m,f) result(x)
      real(dp),intent(in)::p,s,m,f
      x=gamma_quantile(p,s,(m/s)**f)**(1.0_dp/f)
   end function qggamma

   function rggamma(n,s,m,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::s,m,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qggamma(u,s,m,f); end do
   end function rggamma

   real(dp) function bessel_k_nu(z,nu) result(kv)
      real(dp),intent(in)::z,nu
      real(dp)::tmax
      tmax=max(12.0_dp,log(2.0_dp/max(z,1.0e-12_dp))+8.0_dp)
      tmax=min(tmax,40.0_dp)
      kv=integrate_adaptive(integrand,0.0_dp,tmax,2.0e-10_dp,26)
   contains
      function integrand(t) result(v)
         real(dp),intent(in)::t; real(dp)::v,ct,at
         ct=cosh(t); at=-z*ct+log(cosh(nu*t))
         if(at < log(tiny(1.0_dp)))then; v=0.0_dp; else; v=exp(at); end if
      end function integrand
   end function bessel_k_nu

   real(dp) function dginvgauss(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::lv,kn; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<=0.0_dp)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      kn=bessel_k_nu(1.0_dp/(s*m),abs(f))
      lv=(f-1.0_dp)*log(y)-(1.0_dp/y+y/(m*m))/(2.0_dp*s)- &
         f*log(m)-log(2.0_dp*kn)
      v=merge(lv,exp(lv),lg)
   end function dginvgauss

   real(dp) function pginvgauss(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      if(q<=0.0_dp)then; p=0.0_dp; return; end if
      p=integrate_adaptive(integrand,0.0_dp,1.0_dp,2.0e-8_dp,25)
      p=min(1.0_dp,max(0.0_dp,p))
   contains
      function integrand(u) result(v)
         real(dp),intent(in)::u; real(dp)::v,y
         if(u<=0.0_dp)then; v=0.0_dp; return; end if
         y=q*u; v=q*dginvgauss(y,m,s,f)
      end function integrand
   end function pginvgauss

   real(dp) function qginvgauss(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      real(dp)::hi
      hi=max(20.0_dp,10.0_dp*m)
      do while(pginvgauss(hi,m,s,f)<p); hi=2.0_dp*hi; end do
      x=bisect(p,cdf,tiny(1.0_dp),hi)
   contains
      function cdf(z) result(v); real(dp),intent(in)::z; real(dp)::v; v=pginvgauss(z,m,s,f); end function cdf
   end function qginvgauss

   function rginvgauss(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qginvgauss(u,m,s,f); end do
   end function rginvgauss

   elemental real(dp) function pglogis(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      p=(1.0_dp+exp(-sqrt(3.0_dp)*(q-m)/(s*pi)))**(-f)
   end function pglogis

   elemental real(dp) function dglogis(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::y1,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      y1=exp(-sqrt(3.0_dp)*(y-m)/(s*pi))
      lv=0.5_dp*log(3.0_dp)+log(f*y1)-log(pi*s)-(f+1.0_dp)*log(1.0_dp+y1)
      v=merge(lv,exp(lv),lg)
   end function dglogis

   elemental real(dp) function qglogis(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      x=-log(p**(-1.0_dp/f)-1.0_dp)*s*pi/sqrt(3.0_dp)+m
   end function qglogis

   function rglogis(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qglogis(u,m,s,f); end do
   end function rglogis

   elemental real(dp) function pgweibull(q,s,m,f) result(p)
      real(dp),intent(in)::q,s,m,f
      p=(1.0_dp-exp(-(q/m)**s))**f
   end function pgweibull

   elemental real(dp) function dgweibull(y,s,m,f,log_p) result(v)
      real(dp),intent(in)::y,s,m,f; logical,intent(in),optional::log_p
      real(dp)::y1,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      y1=exp(-(y/m)**s)
      lv=log(s*f)+(s-1.0_dp)*log(y)+(f-1.0_dp)*log(1.0_dp-y1)+log(y1)-s*log(m)
      v=merge(lv,exp(lv),lg)
   end function dgweibull

   elemental real(dp) function qgweibull(p,s,m,f) result(x)
      real(dp),intent(in)::p,s,m,f
      x=m*(-log(1.0_dp-p**(1.0_dp/f)))**(1.0_dp/s)
   end function qgweibull

   function rgweibull(n,s,m,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::s,m,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qgweibull(u,s,m,f); end do
   end function rgweibull

   elemental real(dp) function phjorth(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      p=1.0_dp-(1.0_dp+s*q)**(-f/s)*exp(-0.5_dp*(q/m)**2)
   end function phjorth

   elemental real(dp) function dhjorth(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      lv=-(f/s)*log(1.0_dp+s*y)-0.5_dp*(y/m)**2+log(y/(m*m)+f/(1.0_dp+s*y))
      v=merge(lv,exp(lv),lg)
   end function dhjorth

   real(dp) function qhjorth(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      real(dp)::hi
      hi=max(20.0_dp,10.0_dp*m)
      do while(phjorth(hi,m,s,f)<p); hi=2.0_dp*hi; end do
      x=bisect(p,cdf,tiny(1.0_dp),hi)
   contains
      function cdf(z) result(v); real(dp),intent(in)::z; real(dp)::v; v=phjorth(z,m,s,f); end function cdf
   end function qhjorth

   function rhjorth(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qhjorth(u,m,s,f); end do
   end function rhjorth

   elemental real(dp) function dpowexp(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::ss,b,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      ss=sqrt(s); b=1.0_dp+1.0_dp/(2.0_dp*f)
      lv=-(abs(y-m)/ss)**(2.0_dp*f)/2.0_dp-log(ss)-log_gamma(b)-b*log(2.0_dp)
      v=merge(lv,exp(lv),lg)
   end function dpowexp

   elemental real(dp) function ppowexp(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      real(dp)::a,t,side
      if(q==m)then; p=0.5_dp; return; end if
      a=1.0_dp/(2.0_dp*f)
      t=(abs(q-m)/sqrt(s))**(2.0_dp*f)/2.0_dp
      side=0.5_dp*regularized_gamma_p(a,t)
      p=merge(0.5_dp-side,0.5_dp+side,q<m)
   end function ppowexp

   elemental real(dp) function qpowexp(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      real(dp)::a,pr,t,mag
      if(p==0.5_dp)then; x=m; return; end if
      a=1.0_dp/(2.0_dp*f)
      pr=abs(2.0_dp*p-1.0_dp)
      t=gamma_quantile(pr,a,1.0_dp)
      mag=sqrt(s)*(2.0_dp*t)**(1.0_dp/(2.0_dp*f))
      x=merge(m-mag,m+mag,p<0.5_dp)
   end function qpowexp

   function rpowexp(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qpowexp(u,m,s,f); end do
   end function rpowexp

   elemental real(dp) function pskewlaplace(q,m,s,f) result(p)
      real(dp),intent(in)::q,m,s,f
      real(dp)::u
      u=(q-m)/s
      if(u>0.0_dp)then; p=1.0_dp-exp(-f*abs(u))/(1.0_dp+f*f); else
         p=f*f*exp(-abs(u)/f)/(1.0_dp+f*f)
      end if
   end function pskewlaplace

   elemental real(dp) function dskewlaplace(y,m,s,f,log_p) result(v)
      real(dp),intent(in)::y,m,s,f; logical,intent(in),optional::log_p
      real(dp)::lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      lv=log(f)+merge(-f*(y-m),(y-m)/f,y>m)/s-log((1.0_dp+f*f)*s)
      v=merge(lv,exp(lv),lg)
   end function dskewlaplace

   elemental real(dp) function qskewlaplace(p,m,s,f) result(x)
      real(dp),intent(in)::p,m,s,f
      if(p<0.5_dp)then
         x=f*s*log((1.0_dp+f*f)*p/(f*f))+m
      else
         x=-s*log((1.0_dp+f*f)*(1.0_dp-p))/f+m
      end if
   end function qskewlaplace

   function rskewlaplace(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      real(dp),allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qskewlaplace(u,m,s,f); end do
   end function rskewlaplace

   real(dp) function bisect(p,f,lo0,hi0) result(x)
      real(dp),intent(in)::p,lo0,hi0
      procedure(local_cdf)::f
      real(dp)::lo,hi,mid
      integer::i
      lo=lo0; hi=hi0
      do i=1,160
         mid=0.5_dp*(lo+hi)
         if(f(mid)<p)then; lo=mid; else; hi=mid; end if
         if(hi-lo<=2.0e-12_dp*max(1.0_dp,abs(mid)))exit
      end do
      x=0.5_dp*(lo+hi)
   end function bisect

end module rmutil_continuous
