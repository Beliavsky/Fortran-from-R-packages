! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_base
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp, pi, sqrt2, log2pi
   use gamlss_special, only : normal_pdf, normal_cdf, normal_quantile, log1p_v, expm1_v, &
      regularized_gamma_p, regularized_beta, gamma_quantile, log_beta_fn, zeta_fn, &
      harmonic_generalized, lambert_w0
   use gamlss_random, only : random_normal, random_exponential, random_gamma, random_beta, &
      random_poisson, random_binomial
   implicit none
   private
   public :: dnorm_v, pnorm_v, qnorm_v, rnorm_v, dlogis_v, plogis_v, qlogis_v, rlogis_v
   public :: dexp_v, pexp_v, qexp_v, rexp_v, dgamma_v, pgamma_v, qgamma_v, rgamma_v
   public :: dbeta_v, pbeta_v, qbeta_v, rbeta_v, dpois_v, ppois_v, qpois_v, rpois_v
   public :: dbinom_v, pbinom_v, qbinom_v, rbinom_v, dnbinom_v, pnbinom_v, qnbinom_v, rnbinom_v
   public :: dgumbel, pgumbel, qgumbel, rgumbel
   public :: dfrechet, pfrechet, qfrechet, rfrechet
   public :: drayleigh, prayleigh, qrayleigh, rrayleigh
   public :: dpareto4, ppareto4, qpareto4, rpareto4
   public :: dpareto1, ppareto1, qpareto1, dpareto2, ppareto2, qpareto2
   public :: dgenray, pgenray, qgenray, rgenray
   public :: dexpgeom, pexpgeom, qexpgeom, rexpgeom
   public :: dexppois, pexppois, qexppois, rexppois
   public :: dexplog, pexplog, qexplog, rexplog
   public :: dbetabinom, dbetabinom_ab, pbetabinom, rbetabinom
   public :: ddirmultinomial, dzeta, pzeta, qzeta, rzeta, dzipf, pzipf, qzipf, rzipf
   public :: dtriangle, ptriangle, qtriangle, rtriangle
   public :: dinvgaussian, pinvgaussian, qinvgaussian, rinvgaussian

contains

   elemental real(dp) function qnan() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function
   elemental real(dp) function pinf() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function

   elemental real(dp) function dnorm_v(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(sigma<=0.0_dp)then
      v=qnan()
      return
      end if
      ld=-0.5_dp*((x-mu)/sigma)**2-log(sigma)-0.5_dp*log2pi
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pnorm_v(q,mu,sigma) result(v)
      real(dp),intent(in)::q,mu,sigma
      if(sigma<=0.0_dp)then
      v=qnan()
      else
      v=normal_cdf((q-mu)/sigma)
      end if
   end function
   elemental real(dp) function qnorm_v(p,mu,sigma) result(v)
      real(dp),intent(in)::p,mu,sigma
      if(sigma<=0.0_dp)then
      v=qnan()
      else
      v=mu+sigma*normal_quantile(p)
      end if
   end function
   real(dp) function rnorm_v(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma
      v=mu+sigma*random_normal()
   end function

   elemental real(dp) function dlogis_v(x,location,scale,log_density) result(v)
      real(dp),intent(in)::x,location,scale
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then
      v=qnan()
      return
      end if
      z=(x-location)/scale
      ld=-z-2.0_dp*log1p_v(exp(-z))-log(scale)
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function plogis_v(q,location,scale) result(v)
      real(dp),intent(in)::q,location,scale
      real(dp)::z
      if(scale<=0.0_dp)then
      v=qnan()
      return
      end if
      z=(q-location)/scale
      if(z>=0.0_dp)then
      v=1.0_dp/(1.0_dp+exp(-z))
      else
      v=exp(z)/(1.0_dp+exp(z))
      end if
   end function
   elemental real(dp) function qlogis_v(p,location,scale) result(v)
      real(dp),intent(in)::p,location,scale
      if(scale<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp)then
      v=qnan()
      else
      v=location+scale*(log(p)-log1p_v(-p))
      end if
   end function
   real(dp) function rlogis_v(location,scale) result(v)
      real(dp),intent(in)::location,scale
      real(dp)::u
      call random_number(u)
      v=qlogis_v(u,location,scale)
   end function

   elemental real(dp) function dexp_v(x,rate,log_density) result(v)
      real(dp),intent(in)::x,rate
      logical,intent(in),optional::log_density
      real(dp)::ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(rate<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<0.0_dp)then
      ld=-pinf()
      else
      ld=log(rate)-rate*x
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pexp_v(q,rate) result(v)
      real(dp),intent(in)::q,rate
      if(rate<=0.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=-expm1_v(-rate*q)
      end if
   end function
   elemental real(dp) function qexp_v(p,rate) result(v)
      real(dp),intent(in)::p,rate
      if(rate<=0.0_dp .or. p<0.0_dp .or. p>1.0_dp)then
      v=qnan()
      else
      v=-log1p_v(-p)/rate
      end if
   end function
   real(dp) function rexp_v(rate) result(v)
      real(dp),intent(in)::rate
      v=random_exponential(rate)
   end function

   elemental real(dp) function dgamma_v(x,shape,scale,log_density) result(v)
      real(dp),intent(in)::x,shape,scale
      logical,intent(in),optional::log_density
      real(dp)::ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(shape<=0.0_dp .or. scale<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<0.0_dp)then
      ld=-pinf()
      else if(x==0.0_dp)then
         if(shape<1.0_dp)then
      ld=pinf()
      else if(shape==1.0_dp)then
      ld=-log(scale)
      else
      ld=-pinf()
      end if
      else
      ld=(shape-1.0_dp)*log(x)-x/scale-log_gamma(shape)-shape*log(scale)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pgamma_v(q,shape,scale) result(v)
      real(dp),intent(in)::q,shape,scale
      if(shape<=0.0_dp .or. scale<=0.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=regularized_gamma_p(shape,q/scale)
      end if
   end function
   real(dp) function qgamma_v(p,shape,scale) result(v)
      real(dp),intent(in)::p,shape,scale
      v=gamma_quantile(p,shape,scale)
   end function
   real(dp) function rgamma_v(shape,scale) result(v)
      real(dp),intent(in)::shape,scale
      v=random_gamma(shape,scale)
   end function

   elemental real(dp) function dbeta_v(x,a,b,log_density) result(v)
      real(dp),intent(in)::x,a,b
      logical,intent(in),optional::log_density
      real(dp)::ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(a<=0.0_dp .or. b<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp .or. x>=1.0_dp)then
      ld=-pinf()
      else
      ld=(a-1.0_dp)*log(x)+(b-1.0_dp)*log1p_v(-x)-log_beta_fn(a,b)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pbeta_v(q,a,b) result(v)
      real(dp),intent(in)::q,a,b
      if(q<=0.0_dp)then
      v=0.0_dp
      else if(q>=1.0_dp)then
      v=1.0_dp
      else
      v=regularized_beta(q,a,b)
      end if
   end function
   real(dp) function qbeta_v(p,a,b) result(v)
      real(dp),intent(in)::p,a,b
      real(dp)::lo,hi,mid
      integer::i
      if(p<0.0_dp.or.p>1.0_dp.or.a<=0.0_dp.or.b<=0.0_dp)then
      v=qnan()
      return
      end if
      lo=0.0_dp
      hi=1.0_dp
      do i=1,100
      mid=(lo+hi)/2
      if(regularized_beta(mid,a,b)<p)then
      lo=mid
      else
      hi=mid
      end if
      end do
      v=(lo+hi)/2
   end function
   real(dp) function rbeta_v(a,b) result(v)
      real(dp),intent(in)::a,b
      v=random_beta(a,b)
   end function

   elemental real(dp) function dpois_v(x,lambda,log_density) result(v)
      integer,intent(in)::x
      real(dp),intent(in)::lambda
      logical,intent(in),optional::log_density
      real(dp)::ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(lambda<0.0_dp)then
      v=qnan()
      return
      end if
      if(x<0)then
      ld=-pinf()
      else if(lambda==0.0_dp)then
      ld=merge(0.0_dp,-pinf(),x==0)
      else
      ld=-lambda+real(x,dp)*log(lambda)-log_gamma(real(x+1,dp))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function ppois_v(q,lambda) result(v)
      integer,intent(in)::q
      real(dp),intent(in)::lambda
      if(lambda<0.0_dp)then
      v=qnan()
      else if(q<0)then
      v=0.0_dp
      else
      v=1.0_dp-regularized_gamma_p(real(q+1,dp),lambda)
      end if
   end function
   integer function qpois_v(p,lambda) result(q)
      real(dp),intent(in)::p,lambda
      integer::lo,hi,mid
      if(p<0.0_dp.or.p>1.0_dp.or.lambda<0.0_dp)then
      q=-1
      return
      end if
      if(p==1.0_dp)then
      q=huge(q)
      return
      end if
      lo=0
      hi=max(1,nint(lambda+10.0_dp*sqrt(max(lambda,1.0_dp))))
      do while(ppois_v(hi,lambda)<p)
      lo=hi+1
      hi=min(huge(hi)/2,hi*2+1)
      end do
      do while(lo<hi)
      mid=lo+(hi-lo)/2
      if(ppois_v(mid,lambda)>=p)then
      hi=mid
      else
      lo=mid+1
      end if
      end do
      q=lo
   end function
   integer function rpois_v(lambda) result(q)
      real(dp),intent(in)::lambda
      q=random_poisson(lambda)
   end function

   elemental real(dp) function dbinom_v(x,n,p,log_density) result(v)
      integer,intent(in)::x,n
      real(dp),intent(in)::p
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(n<0.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      return
      end if
      if(x<0.or.x>n)then
      ld=-pinf()
      else if(p==0.0_dp)then
      ld=merge(0.0_dp,-pinf(),x==0)
      else if(p==1.0_dp)then
      ld=merge(0.0_dp,-pinf(),x==n)
      else
      ld=log_gamma(real(n+1,dp))-log_gamma(real(x+1,dp))-log_gamma(real(n-x+1,dp))+real(x,dp)*log(p)+real(n-x,dp)*log1p_v(-p)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pbinom_v(q,n,p) result(v)
      integer,intent(in)::q,n
      real(dp),intent(in)::p
      if(n<0.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(q<0)then
      v=0.0_dp
      else if(q>=n)then
      v=1.0_dp
      else
      v=regularized_beta(1.0_dp-p,real(n-q,dp),real(q+1,dp))
      end if
   end function
   integer function qbinom_v(prob,n,p) result(q)
      real(dp),intent(in)::prob,p
      integer,intent(in)::n
      integer::lo,hi,mid
      if(prob<0.0_dp.or.prob>1.0_dp.or.n<0.or.p<0.0_dp.or.p>1.0_dp)then
      q=-1
      return
      end if
      lo=0
      hi=n
      do while(lo<hi)
      mid=(lo+hi)/2
      if(pbinom_v(mid,n,p)>=prob)then
      hi=mid
      else
      lo=mid+1
      end if
      end do
      q=lo
   end function
   integer function rbinom_v(n,p) result(q)
      integer,intent(in)::n
      real(dp),intent(in)::p
      q=random_binomial(n,p)
   end function

   elemental real(dp) function dnbinom_v(x,size,prob,log_density) result(v)
      integer,intent(in)::x
      real(dp),intent(in)::size,prob
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(size<=0.0_dp.or.prob<=0.0_dp.or.prob>1.0_dp)then
      v=qnan()
      return
      end if
      if(x<0)then
      ld=-pinf()
      else
      ld=log_gamma(real(x,dp)+size)-log_gamma(size)-log_gamma(real(x+1,dp))+size*log(prob)+real(x,dp)*log1p_v(-prob)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pnbinom_v(q,size,prob) result(v)
      integer,intent(in)::q
      real(dp),intent(in)::size,prob
      if(q<0)then
      v=0.0_dp
      else
      v=regularized_beta(prob,size,real(q+1,dp))
      end if
   end function
   integer function qnbinom_v(p,size,prob) result(q)
      real(dp),intent(in)::p,size,prob
      integer::lo,hi,mid
      if(p<0.0_dp.or.p>1.0_dp.or.size<=0.0_dp.or.prob<=0.0_dp.or.prob>1.0_dp)then
      q=-1
      return
      end if
      if(p==1.0_dp)then
      q=huge(q)
      return
      end if
      lo=0
      hi=max(1,nint(size*(1.0_dp-prob)/prob+10.0_dp*sqrt(size*(1.0_dp-prob))/prob))
      do while(pnbinom_v(hi,size,prob)<p)
      lo=hi+1
      hi=hi*2+1
      end do
      do while(lo<hi)
      mid=(lo+hi)/2
      if(pnbinom_v(mid,size,prob)>=p)then
      hi=mid
      else
      lo=mid+1
      end if
      end do
      q=lo
   end function
   integer function rnbinom_v(size,prob) result(q)
      real(dp),intent(in)::size,prob
      real(dp)::lambda
      lambda=random_gamma(size,(1.0_dp-prob)/prob)
      q=random_poisson(lambda)
   end function

   elemental real(dp) function dgumbel(x,location,scale,log_density) result(v)
      real(dp),intent(in)::x,location,scale
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::z,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then
      v=qnan()
      return
      end if
      z=(x-location)/scale
      ld=-z-exp(-z)-log(scale)
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pgumbel(q,location,scale) result(v)
      real(dp),intent(in)::q,location,scale
      if(scale<=0.0_dp)then
      v=qnan()
      else
      v=exp(-exp(-(q-location)/scale))
      end if
   end function
   elemental real(dp) function qgumbel(p,location,scale) result(v)
      real(dp),intent(in)::p,location,scale
      if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(p==0.0_dp)then
      v=-pinf()
      else if(p==1.0_dp)then
      v=pinf()
      else
      v=location-scale*log(-log(p))
      end if
   end function
   real(dp) function rgumbel(location,scale) result(v)
      real(dp),intent(in)::location,scale
      real(dp)::u
      call random_number(u)
      v=qgumbel(u,location,scale)
   end function

   elemental real(dp) function dfrechet(x,location,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,location,scale,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::r,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=location)then
      ld=-pinf()
      else
      r=scale/(x-location)
      ld=log(shape)-r**shape+(shape+1.0_dp)*log(r)-log(scale)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pfrechet(q,location,scale,shape) result(v)
      real(dp),intent(in)::q,location,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      else if(q<=location)then
      v=0.0_dp
      else
      v=exp(-(scale/(q-location))**shape)
      end if
   end function
   elemental real(dp) function qfrechet(p,location,scale,shape) result(v)
      real(dp),intent(in)::p,location,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(p==0.0_dp)then
      v=location
      else if(p==1.0_dp)then
      v=pinf()
      else
      v=location+scale*(-log(p))**(-1.0_dp/shape)
      end if
   end function
   real(dp) function rfrechet(location,scale,shape) result(v)
      real(dp),intent(in)::location,scale,shape
      real(dp)::u
      call random_number(u)
      v=qfrechet(u,location,scale,shape)
   end function

   elemental real(dp) function drayleigh(x,scale,log_density) result(v)
      real(dp),intent(in)::x,scale
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp)then
      ld=-pinf()
      else
      ld=log(x)-0.5_dp*(x/scale)**2-2.0_dp*log(scale)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function prayleigh(q,scale) result(v)
      real(dp),intent(in)::q,scale
      if(scale<=0.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=-expm1_v(-0.5_dp*(q/scale)**2)
      end if
   end function
   elemental real(dp) function qrayleigh(p,scale) result(v)
      real(dp),intent(in)::p,scale
      if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else
      v=scale*sqrt(-2.0_dp*log1p_v(-p))
      end if
   end function
   real(dp) function rrayleigh(scale) result(v)
      real(dp),intent(in)::scale
      real(dp)::u
      call random_number(u)
      v=qrayleigh(u,scale)
   end function

   elemental real(dp) function dpareto4(x,location,scale,inequality,shape,log_density) result(v)
      real(dp),intent(in)::x,location,scale,inequality,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::z,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.inequality<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=location)then
      ld=-pinf()
      else
      z=(x-location)/scale
      ld=log(shape)-log(scale)-log(inequality)+(1.0_dp/inequality-1.0_dp)*log(z)-(shape+1.0_dp)*log1p_v(z**(1.0_dp/inequality))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function ppareto4(q,location,scale,inequality,shape) result(v)
      real(dp),intent(in)::q,location,scale,inequality,shape
      real(dp)::z
      if(scale<=0.0_dp.or.inequality<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      else if(q<=location)then
      v=0.0_dp
      else
      z=(q-location)/scale
      v=-expm1_v(-shape*log1p_v(z**(1.0_dp/inequality)))
      end if
   end function
   elemental real(dp) function qpareto4(p,location,scale,inequality,shape) result(v)
      real(dp),intent(in)::p,location,scale,inequality,shape
      if(scale<=0.0_dp.or.inequality<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(p==0.0_dp)then
      v=location
      else if(p==1.0_dp)then
      v=pinf()
      else
      v=location+scale*expm1_v(-log1p_v(-p)/shape)**inequality
      end if
   end function
   real(dp) function rpareto4(location,scale,inequality,shape) result(v)
      real(dp),intent(in)::location,scale,inequality,shape
      real(dp)::u
      call random_number(u)
      v=qpareto4(u,location,scale,inequality,shape)
   end function
   elemental real(dp) function dpareto1(x,scale,shape) result(v)
      real(dp),intent(in)::x,scale,shape
      v=dpareto4(x,scale,scale,1.0_dp,shape)
   end function
   elemental real(dp) function ppareto1(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      v=ppareto4(q,scale,scale,1.0_dp,shape)
   end function
   elemental real(dp) function qpareto1(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      v=qpareto4(p,scale,scale,1.0_dp,shape)
   end function
   elemental real(dp) function dpareto2(x,location,scale,shape) result(v)
      real(dp),intent(in)::x,location,scale,shape
      v=dpareto4(x,location,scale,1.0_dp,shape)
   end function
   elemental real(dp) function ppareto2(q,location,scale,shape) result(v)
      real(dp),intent(in)::q,location,scale,shape
      v=ppareto4(q,location,scale,1.0_dp,shape)
   end function
   elemental real(dp) function qpareto2(p,location,scale,shape) result(v)
      real(dp),intent(in)::p,location,scale,shape
      v=qpareto4(p,location,scale,1.0_dp,shape)
   end function

   elemental real(dp) function dgenray(x,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::t,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp)then
      ld=-pinf()
      else
      t=x/scale
      ld=log(2.0_dp*shape*x)-2.0_dp*log(scale)-t*t+(shape-1.0_dp)*log1p_v(-exp(-t*t))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pgenray(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=(-expm1_v(-(q/scale)**2))**shape
      end if
   end function
   elemental real(dp) function qgenray(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else
      v=scale*sqrt(-log1p_v(-p**(1.0_dp/shape)))
      end if
   end function
   real(dp) function rgenray(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      real(dp)::u
      call random_number(u)
      v=qgenray(u,scale,shape)
   end function

   elemental real(dp) function dexpgeom(x,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::t,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.shape>=1.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp)then
      ld=-pinf()
      else
      t=-x/scale
      ld=-log(scale)+log1p_v(-shape)+t-2.0_dp*log1p_v(-shape*exp(t))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pexpgeom(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      real(dp)::t
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.shape>=1.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      t=-q/scale
      v=-expm1_v(t)/(1.0_dp-shape*exp(t))
      end if
   end function
   elemental real(dp) function qexpgeom(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.shape>=1.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(p==1.0_dp)then
      v=pinf()
      else
      v=-scale*log((p-1.0_dp)/(p*shape-1.0_dp))
      end if
   end function
   real(dp) function rexpgeom(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      real(dp)::u
      call random_number(u)
      v=qexpgeom(u,scale,shape)
   end function

   elemental real(dp) function dexppois(x,rate,shape,log_density) result(v)
      real(dp),intent(in)::x,rate,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(rate<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp)then
      ld=-pinf()
      else
      ld=log(shape*rate)-log1p_v(-exp(-shape))-shape-rate*x+shape*exp(-rate*x)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pexppois(q,rate,shape) result(v)
      real(dp),intent(in)::q,rate,shape
      if(rate<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=(exp(shape*exp(-rate*q))-exp(shape))/(-expm1_v(shape))
      end if
   end function
   elemental real(dp) function qexppois(p,rate,shape) result(v)
      real(dp),intent(in)::p,rate,shape
      if(rate<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(p==1.0_dp)then
      v=pinf()
      else
      v=-log(log(p*(-expm1_v(shape))+exp(shape))/shape)/rate
      end if
   end function
   real(dp) function rexppois(rate,shape) result(v)
      real(dp),intent(in)::rate,shape
      real(dp)::u
      call random_number(u)
      v=qexppois(u,rate,shape)
   end function

   elemental real(dp) function dexplog(x,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::t,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.shape>=1.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp)then
      ld=-pinf()
      else
      t=-x/scale
      ld=-log(-log(shape))-log(scale)+log1p_v(-shape)+t-log1p_v(-(1.0_dp-shape)*exp(t))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pexplog(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.shape>=1.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=1.0_dp-log1p_v(-(1.0_dp-shape)*exp(-q/scale))/log(shape)
      end if
   end function
   elemental real(dp) function qexplog(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.shape>=1.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=qnan()
      else if(p==1.0_dp)then
      v=pinf()
      else
      v=-scale*(log1p_v(-shape**(1.0_dp-p))-log1p_v(-shape))
      end if
   end function
   real(dp) function rexplog(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      real(dp)::u
      call random_number(u)
      v=qexplog(u,scale,shape)
   end function

   elemental real(dp) function dbetabinom_ab(x,size,a,b,log_density) result(v)
      integer,intent(in)::x,size
      real(dp),intent(in)::a,b
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(size<0.or.a<=0.0_dp.or.b<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<0.or.x>size)then
      ld=-pinf()
      else
      ld = log_gamma(real(size + 1, dp)) - log_gamma(real(x + 1, dp)) - &
           log_gamma(real(size - x + 1, dp)) + &
           log_beta_fn(real(x, dp) + a, real(size - x, dp) + b) - &
           log_beta_fn(a, b)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function dbetabinom(x,size,prob,rho,log_density) result(v)
      integer,intent(in)::x,size
      real(dp),intent(in)::prob,rho
      logical,intent(in),optional::log_density
      real(dp)::a,b
      if(prob<=0.0_dp.or.prob>=1.0_dp.or.rho<=0.0_dp.or.rho>=1.0_dp)then
      v=qnan()
      return
      end if
      a=prob*(1.0_dp-rho)/rho
      b=(1.0_dp-prob)*(1.0_dp-rho)/rho
      v=dbetabinom_ab(x,size,a,b,log_density)
   end function
   elemental real(dp) function pbetabinom(q,size,a,b) result(v)
      integer,intent(in)::q,size
      real(dp),intent(in)::a,b
      integer::k
      if(q<0)then
      v=0.0_dp
      return
      else if(q>=size)then
      v=1.0_dp
      return
      end if
      v=0.0_dp
      do k=0,q
      v=v+dbetabinom_ab(k,size,a,b)
      end do
      v=min(1.0_dp,v)
   end function
   integer function rbetabinom(size,a,b) result(x)
      integer,intent(in)::size
      real(dp),intent(in)::a,b
      real(dp)::p
      p=random_beta(a,b)
      x=random_binomial(size,p)
   end function

   pure real(dp) function ddirmultinomial(counts,alpha,log_density) result(v)
      integer,intent(in)::counts(:)
      real(dp),intent(in)::alpha(:)
      logical,intent(in),optional::log_density
      real(dp)::ld,asum
      integer::n,i
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(size(counts)/=size(alpha).or.any(counts<0).or.any(alpha<=0.0_dp))then
      v=qnan()
      return
      end if
      n=sum(counts)
      asum=sum(alpha)
      ld=log_gamma(real(n+1,dp))+log_gamma(asum)-log_gamma(real(n,dp)+asum)
      do i=1,size(alpha)
      ld=ld-log_gamma(real(counts(i)+1,dp))+log_gamma(real(counts(i),dp)+alpha(i))-log_gamma(alpha(i))
      end do
      v=merge(ld,exp(ld),lg)
   end function

   elemental real(dp) function dzeta(x,shape,log_density) result(v)
      integer,intent(in)::x
      real(dp),intent(in)::shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<1)then
      ld=-pinf()
      else
      ld=(-shape-1.0_dp)*log(real(x,dp))-log(zeta_fn(shape+1.0_dp))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   real(dp) function pzeta(q,shape) result(v)
      integer,intent(in)::q
      real(dp),intent(in)::shape
      integer::k
      if(shape<=0.0_dp)then
      v=qnan()
      return
      else if(q<1)then
      v=0.0_dp
      return
      end if
      v=0.0_dp
      do k=1,q
      v=v+real(k,dp)**(-shape-1.0_dp)
      end do
      v=min(1.0_dp,v/zeta_fn(shape+1.0_dp))
   end function
   integer function qzeta(p,shape) result(q)
      real(dp),intent(in)::p,shape
      real(dp)::cum,term,norm
      integer::k
      if(p<0.0_dp.or.p>1.0_dp.or.shape<=0.0_dp)then
      q=-1
      return
      end if
      if(p==1.0_dp)then
      q=huge(q)
      return
      end if
      norm=zeta_fn(shape+1.0_dp)
      cum=0.0_dp
      do k=1,10000000
      term=real(k,dp)**(-shape-1.0_dp)/norm
      cum=cum+term
      if(cum>=p)then
      q=k
      return
      end if
      end do
      q=k
   end function
   integer function rzeta(shape) result(q)
      real(dp),intent(in)::shape
      real(dp) :: u
      call random_number(u)
      q=qzeta(u,shape)
   end function

   elemental real(dp) function dzipf(x,n,shape,log_density) result(v)
      integer,intent(in)::x,n
      real(dp),intent(in)::shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(n<1.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<1.or.x>n)then
      ld=-pinf()
      else
      ld=-shape*log(real(x,dp))-log(harmonic_generalized(n,shape))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   real(dp) function pzipf(q,n,shape) result(v)
      integer,intent(in)::q,n
      real(dp),intent(in)::shape
      if(shape<=0.0_dp.or.n<1)then
      v=qnan()
      else if(q<1)then
      v=0.0_dp
      else if(q>=n)then
      v=1.0_dp
      else
      v=harmonic_generalized(q,shape)/harmonic_generalized(n,shape)
      end if
   end function
   integer function qzipf(p,n,shape) result(q)
      real(dp),intent(in)::p,shape
      integer,intent(in)::n
      integer::lo,hi,mid
      if(p<0.0_dp.or.p>1.0_dp.or.n<1.or.shape<=0.0_dp)then
      q=-1
      return
      end if
      lo=1
      hi=n
      do while(lo<hi)
      mid=(lo+hi)/2
      if(pzipf(mid,n,shape)>=p)then
      hi=mid
      else
      lo=mid+1
      end if
      end do
      q=lo
   end function
   integer function rzipf(n,shape) result(q)
      integer,intent(in)::n
      real(dp),intent(in)::shape
      real(dp) :: u
      call random_number(u)
      q=qzipf(u,n,shape)
   end function

   elemental real(dp) function dtriangle(x,minv,maxv,mode,log_density) result(v)
      real(dp),intent(in)::x,minv,maxv,mode
      logical,intent(in),optional::log_density
      real(dp)::d
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(maxv<=minv.or.mode<minv.or.mode>maxv)then
      v=qnan()
      return
      end if
      if(x<minv.or.x>maxv)then
      d=0.0_dp
      else if(x<=mode)then
      d=2.0_dp*(x-minv)/((maxv-minv)*(mode-minv))
      else
      d=2.0_dp*(maxv-x)/((maxv-minv)*(maxv-mode))
      end if
      if(lg)then
      if(d>0.0_dp)then
      v=log(d)
      else
      v=-pinf()
      end if
      else
      v=d
      end if
   end function
   elemental real(dp) function ptriangle(q,minv,maxv,mode) result(v)
      real(dp),intent(in)::q,minv,maxv,mode
      if(maxv<=minv.or.mode<minv.or.mode>maxv)then
      v=qnan()
      else if(q<=minv)then
      v=0.0_dp
      else if(q>=maxv)then
      v=1.0_dp
      else if(q<=mode)then
      v=(q-minv)**2/((maxv-minv)*(mode-minv))
      else
      v=1.0_dp-(maxv-q)**2/((maxv-minv)*(maxv-mode))
      end if
   end function
   elemental real(dp) function qtriangle(p,minv,maxv,mode) result(v)
      real(dp),intent(in)::p,minv,maxv,mode
      real(dp)::fc
      if(p<0.0_dp.or.p>1.0_dp.or.maxv<=minv.or.mode<minv.or.mode>maxv)then
      v=qnan()
      return
      end if
      fc=(mode-minv)/(maxv-minv)
      if(p<fc)then
      v=minv+sqrt(p*(maxv-minv)*(mode-minv))
      else
      v=maxv-sqrt((1.0_dp-p)*(maxv-minv)*(maxv-mode))
      end if
   end function
   real(dp) function rtriangle(minv,maxv,mode) result(v)
      real(dp),intent(in)::minv,maxv,mode
      real(dp) :: u
      call random_number(u)
      v=qtriangle(u,minv,maxv,mode)
   end function

   elemental real(dp) function dinvgaussian(x,mean,shape,log_density) result(v)
      real(dp),intent(in)::x,mean,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(mean<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(x<=0.0_dp)then
      ld=-pinf()
      else
      ld=0.5_dp*(log(shape)-log(2.0_dp*pi)-3.0_dp*log(x))-shape*(x-mean)**2/(2.0_dp*mean*mean*x)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pinvgaussian(q,mean,shape) result(v)
      real(dp),intent(in)::q,mean,shape
      real(dp)::z1,z2
      if(mean<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      z1=sqrt(shape/q)*(q/mean-1.0_dp)
      z2=-sqrt(shape/q)*(q/mean+1.0_dp)
      v=normal_cdf(z1)+exp(min(700.0_dp,2.0_dp*shape/mean))*normal_cdf(z2)
      v=min(1.0_dp,v)
      end if
   end function
   real(dp) function qinvgaussian(p,mean,shape) result(v)
      real(dp),intent(in)::p,mean,shape
      real(dp)::lo,hi,mid
      integer::i
      if(p<0.0_dp.or.p>1.0_dp.or.mean<=0.0_dp.or.shape<=0.0_dp)then
      v=qnan()
      return
      end if
      if(p==0.0_dp)then
      v=0.0_dp
      return
      else if(p==1.0_dp)then
      v=pinf()
      return
      end if
      lo=0.0_dp
      hi=max(mean,1.0_dp)
      do while(pinvgaussian(hi,mean,shape)<p)
      hi=2.0_dp*hi
      end do
      do i=1,100
      mid=(lo+hi)/2
      if(pinvgaussian(mid,mean,shape)<p)then
      lo=mid
      else
      hi=mid
      end if
      end do
      v=(lo+hi)/2
   end function
   real(dp) function rinvgaussian(mean,shape) result(v)
      real(dp),intent(in)::mean,shape
      real(dp)::y,x,u
      y=random_normal()**2
      x=mean+mean*mean*y/(2.0_dp*shape)-mean/(2.0_dp*shape)*sqrt(4.0_dp*mean*shape*y+mean*mean*y*y)
      call random_number(u)
      if(u<=mean/(mean+x))then
      v=x
      else
      v=mean*mean/x
      end if
   end function

end module gamlss_base
