! Computational translation of continuous families from gamlss.dist 6.1-1.
! Original package: Copyright (C) gamlss.dist authors, GPL-2 | GPL-3.
! Translation: 2026. SPDX-License-Identifier: GPL-3.0-only
module gamlss_continuous
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp, pi
   use gamlss_special, only : log1p_v, expm1_v
   use gamlss_base, only : dnorm_v, pnorm_v, qnorm_v, rnorm_v, dlogis_v, plogis_v, qlogis_v, rlogis_v, &
      dexp_v, pexp_v, qexp_v, rexp_v, dgamma_v, pgamma_v, qgamma_v, rgamma_v, &
      dbeta_v, pbeta_v, qbeta_v, rbeta_v, dinvgaussian, pinvgaussian, qinvgaussian, rinvgaussian
   use gamlss_student_t, only : student_t_pdf, student_t_cdf, student_t_quantile, random_student_t
   implicit none
   private
   public :: dNO,pNO,qNO,rNO,dNO2,pNO2,qNO2,rNO2
   public :: dEXP,pEXP,qEXP,rEXP,dGA,pGA,qGA,rGA
   public :: dWEI,pWEI,qWEI,rWEI,dWEI2,pWEI2,qWEI2,rWEI2,dWEI3,pWEI3,qWEI3,rWEI3
   public :: dBE,pBE,qBE,rBE,dBEo,pBEo,qBEo,rBEo
   public :: dBEINF,pBEINF,qBEINF,rBEINF,dBEINF0,pBEINF0,qBEINF0,rBEINF0
   public :: dBEINF1,pBEINF1,qBEINF1,rBEINF1
   public :: dLO,pLO,qLO,rLO,dGU,pGU,qGU,rGU
   public :: dLOGNO,pLOGNO,qLOGNO,rLOGNO,dLOGNO2,pLOGNO2,qLOGNO2,rLOGNO2
   public :: dLNO,pLNO,qLNO,rLNO
   public :: dTF,pTF,qTF,rTF,dPE,pPE,qPE,rPE,dPE2,pPE2,qPE2,rPE2
   public :: dEGB2,pEGB2,qEGB2,rEGB2,dGB2,pGB2,qGB2,rGB2,dGB1,pGB1,qGB1,rGB1
   public :: dGG,pGG,qGG,rGG,dIG,pIG,qIG,rIG,dIGAMMA,pIGAMMA,qIGAMMA,rIGAMMA
   public :: dGP,pGP,qGP,rGP,dJSUo,pJSUo,qJSUo,rJSUo,dJSU,pJSU,qJSU,rJSU

contains

   elemental real(dp) function nanv() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function
   elemental real(dp) function infv() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function
   elemental logical function want_log(flag) result(v)
      logical,intent(in),optional::flag
      v=.false.
      if(present(flag))v=flag
   end function
   elemental logical function lower(flag) result(v)
      logical,intent(in),optional::flag
      v=.true.
      if(present(flag))v=flag
   end function
   elemental real(dp) function finish_prob(p,lower_tail,log_p) result(v)
      real(dp),intent(in)::p
      logical,intent(in),optional::lower_tail,log_p
      v=p
      if(.not.lower(lower_tail))v=1.0_dp-v
      v=max(0.0_dp,min(1.0_dp,v))
      if(present(log_p))then
         if(log_p)then
            if(v==0.0_dp)then
            v=-infv()
            else
            v=log(v)
            end if
         end if
      end if
   end function
   elemental real(dp) function input_prob(p,lower_tail,log_p) result(v)
      real(dp),intent(in)::p
      logical,intent(in),optional::lower_tail,log_p
      v=p
      if(present(log_p))then
         if(log_p)v=exp(v)
      end if
      if(.not.lower(lower_tail))v=1.0_dp-v
   end function

! Normal, NO2 (variance parameter), exponential, gamma.
   elemental real(dp) function dNO(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      v=dnorm_v(x,mu,sigma,want_log(log_density))
   end function
   elemental real(dp) function pNO(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      v=finish_prob(pnorm_v(q,mu,sigma),lower_tail,log_p)
   end function
   elemental real(dp) function qNO(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      v=qnorm_v(input_prob(p,lower_tail,log_p),mu,sigma)
   end function
   real(dp) function rNO(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma
      v=rnorm_v(mu,sigma)
   end function
   elemental real(dp) function dNO2(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      if(sigma<=0)then
      v=nanv()
      else
      v=dnorm_v(x,mu,sqrt(sigma),want_log(log_density))
      end if
   end function
   elemental real(dp) function pNO2(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(sigma<=0)then
      v=nanv()
      else
      v=finish_prob(pnorm_v(q,mu,sqrt(sigma)),lower_tail,log_p)
      end if
   end function
   elemental real(dp) function qNO2(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(sigma<=0)then
      v=nanv()
      else
      v=qnorm_v(input_prob(p,lower_tail,log_p),mu,sqrt(sigma))
      end if
   end function
   real(dp) function rNO2(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma
      if(sigma<=0)then
      v=nanv()
      else
      v=rnorm_v(mu,sqrt(sigma))
      end if
   end function
   elemental real(dp) function dEXP(x,mu,log_density) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_density
      if(mu<=0)then
      v=nanv()
      else if(x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else
      v=dexp_v(x,1/mu,want_log(log_density))
      end if
   end function
   elemental real(dp) function pEXP(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0)then
      v=nanv()
      else
      v=finish_prob(pexp_v(q,1/mu),lower_tail,log_p)
      end if
   end function
   elemental real(dp) function qEXP(p,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0)then
      v=nanv()
      else
      v=qexp_v(input_prob(p,lower_tail,log_p),1/mu)
      end if
   end function
   real(dp) function rEXP(mu) result(v)
      real(dp),intent(in)::mu
      if(mu<=0)then
      v=nanv()
      else
      v=rexp_v(1/mu)
      end if
   end function
   elemental real(dp) function dGA(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else if(x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         v=dgamma_v(x,1/(sigma*sigma),mu*sigma*sigma,want_log(log_density))
         end if
   end function
   elemental real(dp) function pGA(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=finish_prob(pgamma_v(q,1/(sigma*sigma),mu*sigma*sigma),lower_tail,log_p)
      end if
   end function
   real(dp) function qGA(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=qgamma_v(input_prob(p,lower_tail,log_p),1/(sigma*sigma),mu*sigma*sigma)
      end if
   end function
   real(dp) function rGA(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=rgamma_v(1/(sigma*sigma),mu*sigma*sigma)
      end if
   end function

! Weibull parameterizations.
   elemental real(dp) function dweib(x,scale,shape,lg) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in)::lg
      real(dp)::ld
      if(scale<=0.or.shape<=0)then
      v=nanv()
      return
      end if
      if(x<=0)then
      v=merge(-infv(),0.0_dp,lg)
      return
      end if
      ld=log(shape)-log(scale)+(shape-1)*log(x/scale)-(x/scale)**shape
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pweib(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      if(scale<=0.or.shape<=0)then
      v=nanv()
      else if(q<=0)then
      v=0
      else
      v=-expm1_v(-(q/scale)**shape)
      end if
   end function
   elemental real(dp) function qweib(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      if(scale<=0.or.shape<=0.or.p<0.or.p>1)then
      v=nanv()
      else if(p==1)then
      v=infv()
      else
      v=scale*(-log1p_v(-p))**(1/shape)
      end if
   end function
   real(dp) function rweib(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      real(dp)::u
      call random_number(u)
      v=qweib(u,scale,shape)
   end function
   elemental real(dp) function dWEI(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      v=dweib(x,mu,sigma,want_log(log_density))
   end function
   elemental real(dp) function pWEI(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      v=finish_prob(pweib(q,mu,sigma),lower_tail,log_p)
   end function
   elemental real(dp) function qWEI(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      v=qweib(input_prob(p,lower_tail,log_p),mu,sigma)
   end function
   real(dp) function rWEI(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   v=rweib(mu,sigma)
   end function
   elemental real(dp) function dWEI2(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=dweib(x,mu**(-1/sigma),sigma,want_log(log_density))
      end if
   end function
   elemental real(dp) function pWEI2(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=finish_prob(pweib(q,mu**(-1/sigma),sigma),lower_tail,log_p)
      end if
   end function
   elemental real(dp) function qWEI2(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=qweib(input_prob(p,lower_tail,log_p),mu**(-1/sigma),sigma)
      end if
   end function
   real(dp) function rWEI2(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=rweib(mu**(-1/sigma),sigma)
      end if
   end function
   elemental real(dp) function dWEI3(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::sc
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      sc=mu/gamma(1+1/sigma)
      v=dweib(x,sc,sigma,want_log(log_density))
      end if
   end function
   elemental real(dp) function pWEI3(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::sc
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      sc=mu/gamma(1+1/sigma)
      v=finish_prob(pweib(q,sc,sigma),lower_tail,log_p)
      end if
   end function
   elemental real(dp) function qWEI3(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::sc
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      sc=mu/gamma(1+1/sigma)
      v=qweib(input_prob(p,lower_tail,log_p),sc,sigma)
      end if
   end function
   real(dp) function rWEI3(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=rweib(mu/gamma(1+1/sigma),sigma)
      end if
   end function

! Beta and beta-inflated families.
   elemental subroutine be_shapes(mu,sigma,a,b,ok)
      real(dp),intent(in)::mu,sigma
      real(dp),intent(out)::a,b
      logical,intent(out)::ok
      ok=mu>0.and.mu<1.and.sigma>0.and.sigma<1
      if(ok)then
      a=mu*(1-sigma*sigma)/(sigma*sigma)
      b=a*(1-mu)/mu
      else
      a=1
      b=1
      end if
   end subroutine
   elemental real(dp) function dBE(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::a,b
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok)then
      v=nanv()
      else
      v=dbeta_v(x,a,b,want_log(log_density))
      end if
   end function
   elemental real(dp) function pBE(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok)then
      v=nanv()
      else
      v=finish_prob(pbeta_v(q,a,b),lower_tail,log_p)
      end if
   end function
   real(dp) function qBE(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok)then
      v=nanv()
      else
      v=qbeta_v(input_prob(p,lower_tail,log_p),a,b)
      end if
   end function
   real(dp) function rBE(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::a,b
   logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok)then
      v=nanv()
      else
      v=rbeta_v(a,b)
      end if
   end function
   elemental real(dp) function dBEo(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=dbeta_v(x,mu,sigma,want_log(log_density))
      end if
   end function
   elemental real(dp) function pBEo(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=finish_prob(pbeta_v(q,mu,sigma),lower_tail,log_p)
      end if
   end function
   real(dp) function qBEo(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      v=qbeta_v(input_prob(p,lower_tail,log_p),mu,sigma)
      end if
   end function
   real(dp) function rBEo(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   if(mu<=0.or.sigma<=0)then
   v=nanv()
   else
   v=rbeta_v(mu,sigma)
   end if
   end function

   elemental real(dp) function dBEINF(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::a,b,ld
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok.or.nu<=0.or.tau<=0)then
      v=nanv()
      return
      end if
      if(x<0.or.x>1)then
      ld=-infv()
      else if(x==0)then
      ld=log(nu)-log(1+nu+tau)
      else if(x==1)then
      ld=log(tau)-log(1+nu+tau)
      else
         ld=dbeta_v(x,a,b,.true.)-log(1+nu+tau)
         end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pBEINF(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,p
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok.or.nu<=0.or.tau<=0)then
      v=nanv()
      return
      end if
      if(q<0)then
      p=0
      else if(q<1)then
      p=(nu+pbeta_v(q,a,b))/(1+nu+tau)
      else
      p=1
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qBEINF(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,u,p0,p1
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      u=input_prob(p,lower_tail,log_p)
      if(.not.ok.or.nu<=0.or.tau<=0.or.u<0.or.u>1)then
      v=nanv()
      return
      end if
      p0=nu/(1+nu+tau)
      p1=(1+nu)/(1+nu+tau)
      if(u<=p0)then
      v=0
      else if(u>=p1)then
      v=1
      else
      v=qbeta_v(u*(1+nu+tau)-nu,a,b)
      end if
   end function
   real(dp) function rBEINF(mu,sigma,nu,tau) result(v)
   real(dp),intent(in)::mu,sigma,nu,tau
   real(dp)::u
   call random_number(u)
   v=qBEINF(u,mu,sigma,nu,tau)
   end function
   elemental real(dp) function dBEINF0(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      v=dBEINF(x,mu,sigma,nu,tiny(1.0_dp),log_density)
   end function
   elemental real(dp) function pBEINF0(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,p
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok.or.nu<=0)then
      v=nanv()
      return
      end if
      if(q<0)then
      p=0
      else if(q<1)then
      p=(nu+pbeta_v(q,a,b))/(1+nu)
      else
      p=1
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qBEINF0(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,u,p0
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      u=input_prob(p,lower_tail,log_p)
      if(.not.ok.or.nu<=0)then
      v=nanv()
      return
      end if
      p0=nu/(1+nu)
      if(u<=p0)then
      v=0
      else
      v=qbeta_v(u*(1+nu)-nu,a,b)
      end if
   end function
   real(dp) function rBEINF0(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   real(dp)::u
   call random_number(u)
   v=qBEINF0(u,mu,sigma,nu)
   end function
   elemental real(dp) function dBEINF1(x,mu,sigma,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,tau
      logical,intent(in),optional::log_density
      real(dp)::a,b,ld
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok.or.tau<=0)then
      v=nanv()
      return
      end if
      if(x<0.or.x>1)then
      ld=-infv()
      else if(x==1)then
      ld=log(tau)-log(1+tau)
      else if(x==0)then
      ld=-infv()
      else
      ld=dbeta_v(x,a,b,.true.)-log(1+tau)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pBEINF1(q,mu,sigma,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,p
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      if(.not.ok.or.tau<=0)then
      v=nanv()
      return
      end if
      if(q<=0)then
      p=0
      else if(q<1)then
      p=pbeta_v(q,a,b)/(1+tau)
      else
      p=1
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qBEINF1(p,mu,sigma,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,u,pc
      logical::ok
      call be_shapes(mu,sigma,a,b,ok)
      u=input_prob(p,lower_tail,log_p)
      if(.not.ok.or.tau<=0)then
      v=nanv()
      return
      end if
      pc=1/(1+tau)
      if(u>=pc)then
      v=1
      else
      v=qbeta_v(u*(1+tau),a,b)
      end if
   end function
   real(dp) function rBEINF1(mu,sigma,tau) result(v)
   real(dp),intent(in)::mu,sigma,tau
   real(dp)::u
   call random_number(u)
   v=qBEINF1(u,mu,sigma,tau)
   end function

! Logistic, Gumbel minimum, lognormal and Box-Cox normal.
   elemental real(dp) function dLO(x,mu,sigma,log_density) result(v)
   real(dp),intent(in)::x,mu,sigma
   logical,intent(in),optional::log_density
   v=dlogis_v(x,mu,sigma,want_log(log_density))
   end function
   elemental real(dp) function pLO(q,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::q,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
   v=finish_prob(plogis_v(q,mu,sigma),lower_tail,log_p)
   end function
   elemental real(dp) function qLO(p,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
   v=qlogis_v(input_prob(p,lower_tail,log_p),mu,sigma)
   end function
   real(dp) function rLO(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   v=rlogis_v(mu,sigma)
   end function
   elemental real(dp) function dGU(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(sigma<=0)then
      v=nanv()
      return
      end if
      z=(x-mu)/sigma
      ld=-log(sigma)+z-exp(z)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pGU(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(sigma<=0)then
      v=nanv()
      else
      p=-expm1_v(-exp((q-mu)/sigma))
      v=finish_prob(p,lower_tail,log_p)
      end if
   end function
   elemental real(dp) function qGU(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u
      u=input_prob(p,lower_tail,log_p)
      if(sigma<=0.or.u<0.or.u>1)then
      v=nanv()
      else
      v=mu+sigma*log(-log1p_v(-u))
      end if
   end function
   real(dp) function rGU(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   v=qGU(u,mu,sigma)
   end function
   elemental real(dp) function dLOGNO(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(sigma<=0)then
      v=nanv()
      else if(x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else
      ld=dnorm_v(log(x),mu,sigma,.true.)-log(x)
      v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function
   elemental real(dp) function pLOGNO(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(sigma<=0)then
      v=nanv()
      else
      p=merge(pnorm_v(log(q),mu,sigma),0.0_dp,q>0)
      v=finish_prob(p,lower_tail,log_p)
      end if
   end function
   elemental real(dp) function qLOGNO(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      v=exp(qnorm_v(input_prob(p,lower_tail,log_p),mu,sigma))
   end function
   real(dp) function rLOGNO(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   v=exp(rnorm_v(mu,sigma))
   end function
   elemental real(dp) function dLOGNO2(x,mu,sigma,log_density) result(v)
   real(dp),intent(in)::x,mu,sigma
   logical,intent(in),optional::log_density
      if(mu<=0)then
      v=nanv()
      else
      v=dLOGNO(x,log(mu),sigma,log_density)
      end if
      end function
   elemental real(dp) function pLOGNO2(q,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::q,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
      if(mu<=0)then
      v=nanv()
      else
      v=pLOGNO(q,log(mu),sigma,lower_tail,log_p)
      end if
      end function
   elemental real(dp) function qLOGNO2(p,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
      if(mu<=0)then
      v=nanv()
      else
      v=qLOGNO(p,log(mu),sigma,lower_tail,log_p)
      end if
      end function
   real(dp) function rLOGNO2(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   if(mu<=0)then
   v=nanv()
   else
   v=rLOGNO(log(mu),sigma)
   end if
   end function
   elemental real(dp) function boxcox(x,nu) result(z)
   real(dp),intent(in)::x,nu
   if(abs(nu)>1e-12_dp)then
   z=(x**nu-1)/nu
   else
   z=log(x)
   end if
   end function
   elemental real(dp) function invboxcox(z,nu) result(x)
   real(dp),intent(in)::z,nu
   if(abs(nu)>1e-12_dp)then
   if(1+nu*z>0)then
   x=(1+nu*z)**(1/nu)
   else
   x=nanv()
   end if
   else
   x=exp(z)
   end if
   end function
   elemental real(dp) function dLNO(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(sigma<=0.or.x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      z=boxcox(x,nu)
      ld=dnorm_v(z,mu,sigma,.true.)+(nu-1)*log(x)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pLNO(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(q<=0)then
      p=0
      else
      p=pnorm_v(boxcox(q,nu),mu,sigma)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   elemental real(dp) function qLNO(p,mu,sigma,nu,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu,sigma,nu
   logical,intent(in),optional::lower_tail,log_p
   v=invboxcox(qnorm_v(input_prob(p,lower_tail,log_p),mu,sigma),nu)
   end function
   real(dp) function rLNO(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   v=invboxcox(rnorm_v(mu,sigma),nu)
   end function

! Student t and exponential-power parameterizations.
   elemental real(dp) function dTF(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::d
      if(sigma<=0.or.nu<=0)then
      v=nanv()
      return
      end if
      if(nu>1e6_dp)then
      v=dnorm_v(x,mu,sigma,want_log(log_density))
      else
      d=student_t_pdf((x-mu)/sigma,nu)/sigma
      v=merge(log(d),d,want_log(log_density))
      end if
   end function
   elemental real(dp) function pTF(q,mu,sigma,nu,lower_tail,log_p) result(v)
   real(dp),intent(in)::q,mu,sigma,nu
   logical,intent(in),optional::lower_tail,log_p
   real(dp)::p
      if(nu>1e6_dp)then
      p=pnorm_v(q,mu,sigma)
      else
      p=student_t_cdf((q-mu)/sigma,nu)
      end if
      v=finish_prob(p,lower_tail,log_p)
      end function
   real(dp) function qTF(p,mu,sigma,nu,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu,sigma,nu
   logical,intent(in),optional::lower_tail,log_p
   real(dp)::u
      u=input_prob(p,lower_tail,log_p)
      if(nu>1e6_dp)then
      v=qnorm_v(u,mu,sigma)
      else
      v=mu+sigma*student_t_quantile(u,nu)
      end if
      end function
   real(dp) function rTF(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   if(nu>1e6_dp)then
   v=rnorm_v(mu,sigma)
   else
   v=mu+sigma*random_student_t(nu)
   end if
   end function
   elemental real(dp) function dPE2(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(sigma<=0.or.nu<=0)then
      v=nanv()
      return
      end if
      z=(x-mu)/sigma
      ld=-log(sigma)+log(nu)-log(2.0_dp)-abs(z)**nu-log_gamma(1/nu)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pPE2(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,p
      z=(q-mu)/sigma
      if(abs(z)<tiny(1.0_dp))then
      p=.5_dp
      else
      p=.5_dp*(1+sign(1.0_dp,z)*pgamma_v(abs(z)**nu,1/nu,1.0_dp))
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qPE2(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,s
      u=input_prob(p,lower_tail,log_p)
      if(u==.5_dp)then
      v=mu
      else
      s=qgamma_v(abs(2*u-1),1/nu,1.0_dp)
      v=mu+sigma*sign(s**(1/nu),u-.5_dp)
      end if
   end function
   real(dp) function rPE2(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   real(dp)::u
   call random_number(u)
   v=qPE2(u,mu,sigma,nu)
   end function
   elemental real(dp) function dPE(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::c,z,ld
      if(sigma<=0.or.nu<=0)then
      v=nanv()
      return
      end if
      c=exp(.5_dp*(-(2/nu)*log(2.0_dp)+log_gamma(1/nu)-log_gamma(3/nu)))
      z=(x-mu)/sigma
      ld=-log(sigma)+log(nu)-log(c)-.5_dp*(abs(z/c)**nu)-(.5_dp+1/nu)*log(2.0_dp)-log_gamma(1/nu)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pPE(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::c,z,p
      c=exp(.5_dp*(-(2/nu)*log(2.0_dp)+log_gamma(1/nu)-log_gamma(3/nu)))
      z=(q-mu)/sigma
      if(abs(z)<tiny(1.0_dp))then
      p=.5_dp
      else
      p=.5_dp*(1+sign(1.0_dp,z)*pgamma_v(.5_dp*abs(z/c)**nu,1/nu,1.0_dp))
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qPE(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,s,c
      u=input_prob(p,lower_tail,log_p)
      c=exp(.5_dp*(-(2/nu)*log(2.0_dp)+log_gamma(1/nu)-log_gamma(3/nu)))
      if(u==.5_dp)then
      v=mu
      else
      s=qgamma_v(abs(2*u-1),1/nu,1.0_dp)
      v=mu+sigma*sign((2*s)**(1/nu)*c,u-.5_dp)
      end if
   end function
   real(dp) function rPE(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   real(dp)::u
   call random_number(u)
   v=qPE(u,mu,sigma,nu)
   end function

! EGB2, GB1, GB2 and generalized gamma.
   elemental real(dp) function dEGB2(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(sigma==0.or.nu<=0.or.tau<=0)then
      v=nanv()
      return
      end if
      z=(x-mu)/sigma
      ld=nu*z-log(abs(sigma))-log_gamma(nu)-log_gamma(tau)+log_gamma(nu+tau)-(nu+tau)*log(1+exp(z))
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pEGB2(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,t,p
      z=(q-mu)/sigma
      t=1/(1+exp(-z))
      p=pbeta_v(t,nu,tau)
      if(sigma<0)p=1-p
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qEGB2(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,t
      u=input_prob(p,lower_tail,log_p)
      if(sigma<0)u=1-u
      t=qbeta_v(u,nu,tau)
      v=mu+sigma*log(t/(1-t))
   end function
   real(dp) function rEGB2(mu,sigma,nu,tau) result(v)
   real(dp),intent(in)::mu,sigma,nu,tau
   real(dp)::u
   call random_number(u)
   v=qEGB2(u,mu,sigma,nu,tau)
   end function
   elemental real(dp) function dGB2(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(mu<=0.or.sigma==0.or.nu<=0.or.tau<=0)then
      v=nanv()
      return
      end if
      if(x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      z=(x/mu)**sigma
      ld=nu*log(z)+log(abs(sigma))-log(x)-log_gamma(nu)-log_gamma(tau)+log_gamma(nu+tau)-(nu+tau)*log(1+z)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pGB2(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,p
      if(q<=0)then
      p=0
      else
      z=(q/mu)**sigma
      p=pbeta_v(z/(1+z),nu,tau)
      if(sigma<0)p=1-p
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qGB2(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,t
      u=input_prob(p,lower_tail,log_p)
      if(sigma<0)u=1-u
      t=qbeta_v(u,nu,tau)
      v=mu*(t/(1-t))**(1/sigma)
   end function
   real(dp) function rGB2(mu,sigma,nu,tau) result(v)
   real(dp),intent(in)::mu,sigma,nu,tau
   real(dp)::u
   call random_number(u)
   v=qGB2(u,mu,sigma,nu,tau)
   end function
   elemental real(dp) function dGB1(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::a,b,ld
      if(mu<=0.or.mu>=1.or.sigma<=0.or.sigma>=1.or.nu<=0.or.tau<=0)then
      v=nanv()
      return
      end if
      if(x<=0.or.x>=1)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      a=mu*(1/(sigma*sigma)-1)
      b=a*(1-mu)/mu
      ld=log(tau)+b*log(nu)+(tau*a-1)*log(x)+(b-1)*log(1-x**tau) &
         -log_gamma(a)-log_gamma(b)+log_gamma(a+b) &
         -(a+b)*log(nu+(1-nu)*x**tau)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pGB1(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,z,p
      if(q<=0)then
      p=0
      else if(q>=1)then
      p=1
      else
      a=mu*(1/(sigma*sigma)-1)
      b=a*(1-mu)/mu
      z=q**tau/(nu+(1-nu)*q**tau)
      p=pbeta_v(z,a,b)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qGB1(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,b,u,z
      a=mu*(1/(sigma*sigma)-1)
      b=a*(1-mu)/mu
      u=input_prob(p,lower_tail,log_p)
      z=qbeta_v(u,a,b)
      v=(nu/(1/z-(1-nu)))**(1/tau)
   end function
   real(dp) function rGB1(mu,sigma,nu,tau) result(v)
   real(dp),intent(in)::mu,sigma,nu,tau
   real(dp)::u
   call random_number(u)
   v=qGB1(u,mu,sigma,nu,tau)
   end function
   elemental real(dp) function dGG(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(mu<=0.or.sigma<=0.or.x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      if(abs(nu)<=1e-6_dp)then
      ld=-log(x)-.5_dp*log(2*pi)-log(sigma)-.5_dp*((log(x)-log(mu))/sigma)**2
      else
      z=(x/mu)**nu
      ld=dGA(z,1.0_dp,sigma*abs(nu),.true.)+log(abs(nu)*z/x)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pGG(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,p
      if(q<=0)then
      p=0
      else if(abs(nu)<=1e-6_dp)then
      p=pnorm_v(log(q),log(mu),sigma)
      else
      z=(q/mu)**nu
      p=pGA(z,1.0_dp,sigma*abs(nu))
      if(nu<0)p=1-p
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qGG(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,z
      u=input_prob(p,lower_tail,log_p)
      if(abs(nu)<=1e-6_dp)then
      v=exp(qnorm_v(u,log(mu),sigma))
      else
      if(nu<0)u=1-u
      z=qGA(u,1.0_dp,sigma*abs(nu))
      v=mu*z**(1/nu)
      end if
   end function
   real(dp) function rGG(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   real(dp)::u
   call random_number(u)
   v=qGG(u,mu,sigma,nu)
   end function

! Inverse Gaussian (GAMLSS sigma parameter) and inverse gamma.
   elemental real(dp) function dIG(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::shape
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      shape=1/(sigma*sigma)
      v=dinvgaussian(x,mu,shape,want_log(log_density))
      end if
   end function
   elemental real(dp) function pIG(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      else
      p=pinvgaussian(q,mu,1/(sigma*sigma))
      v=finish_prob(p,lower_tail,log_p)
      end if
   end function
   real(dp) function qIG(p,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
      v=qinvgaussian(input_prob(p,lower_tail,log_p),mu,1/(sigma*sigma))
      end function
   real(dp) function rIG(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   v=rinvgaussian(mu,1/(sigma*sigma))
   end function
   elemental real(dp) function dIGAMMA(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::a,ld
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      return
      end if
      if(x<=0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      a=1/(sigma*sigma)
      ld=a*log(mu)+a*log(a+1)-log_gamma(a)-(a+1)*log(x)-mu*(a+1)/x
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pIGAMMA(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::a,p
      if(q<=0)then
      p=0
      else
      a=1/(sigma*sigma)
      p=1-pgamma_v(mu*(a+1)/q,a,1.0_dp)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   real(dp) function qIGAMMA(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,a,g
      u=input_prob(p,lower_tail,log_p)
      a=1/(sigma*sigma)
      g=qgamma_v(1-u,a,1.0_dp)
      v=mu*(a+1)/g
   end function
   real(dp) function rIGAMMA(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   v=qIGAMMA(u,mu,sigma)
   end function

! GAMLSS Generalised Pareto (Lomax/Beta-prime form).
   elemental real(dp) function dGP(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      return
      end if
      if(x<0)then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      ld=log(sigma)-log(mu)-(sigma+1)*log1p_v(x/mu)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pGP(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(q<0)then
      p=0
      else
      p=1-(mu/(mu+q))**sigma
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   elemental real(dp) function qGP(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u
      u=input_prob(p,lower_tail,log_p)
      if(u==1)then
      v=infv()
      else
      v=mu*((1-u)**(-1/sigma)-1)
      end if
   end function
   real(dp) function rGP(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   v=qGP(u,mu,sigma)
   end function

! Johnson SU original and moment-standardized GAMLSS forms.
   elemental real(dp) function dJSUo(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::z,r,ld
      if(sigma<=0.or.tau<=0)then
      v=nanv()
      return
      end if
      z=(x-mu)/sigma
      r=nu+tau*asinh(z)
      ld=-log(sigma)+log(tau)-.5_dp*log(z*z+1)-.5_dp*log(2*pi)-.5_dp*r*r
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pJSUo(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::r
      r=nu+tau*asinh((q-mu)/sigma)
      v=finish_prob(pnorm_v(r,0.0_dp,1.0_dp),lower_tail,log_p)
   end function
   elemental real(dp) function qJSUo(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::r
      r=qnorm_v(input_prob(p,lower_tail,log_p),0.0_dp,1.0_dp)
      v=mu+sigma*sinh((r-nu)/tau)
   end function
   real(dp) function rJSUo(mu,sigma,nu,tau) result(v)
   real(dp),intent(in)::mu,sigma,nu,tau
   real(dp)::u
   call random_number(u)
   v=qJSUo(u,mu,sigma,nu,tau)
   end function
   elemental subroutine jsu_consts(nu,tau,rtau,w,omega,c,shift,ok)
      real(dp),intent(in)::nu,tau
      real(dp),intent(out)::rtau,w,omega,c,shift
      logical,intent(out)::ok
      ok=tau>0
      if(.not.ok)then
      rtau=0
      w=0
      omega=0
      c=0
      shift=0
      return
      end if
      rtau=1/tau
      if(rtau<1e-7_dp)then
      w=1
      else
      w=exp(rtau*rtau)
      end if
      omega=-nu*rtau
      c=(.5_dp*(w-1)*(w*cosh(2*omega)+1))**(-.5_dp)
      shift=c*sqrt(w)*sinh(omega)
   end subroutine
   elemental real(dp) function dJSU(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::rt,w,o,c,s,z,r,ld
      logical::ok
      call jsu_consts(nu,tau,rt,w,o,c,s,ok)
      if(.not.ok.or.sigma<=0)then
      v=nanv()
      return
      end if
      z=(x-(mu+sigma*s))/(c*sigma)
      r=-nu+asinh(z)/rt
      ld=-log(sigma)-log(c)-log(rt)-.5_dp*log(z*z+1)-.5_dp*log(2*pi)-.5_dp*r*r
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pJSU(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::rt,w,o,c,s,z,r
      logical::ok
      call jsu_consts(nu,tau,rt,w,o,c,s,ok)
      z=(q-(mu+sigma*s))/(c*sigma)
      r=-nu+asinh(z)/rt
      v=finish_prob(pnorm_v(r,0.0_dp,1.0_dp),lower_tail,log_p)
   end function
   elemental real(dp) function qJSU(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::rt,w,o,c,s,r,z
      logical::ok
      call jsu_consts(nu,tau,rt,w,o,c,s,ok)
      r=qnorm_v(input_prob(p,lower_tail,log_p),0.0_dp,1.0_dp)
      z=sinh(rt*(r+nu))
      v=mu+sigma*s+c*sigma*z
   end function
   real(dp) function rJSU(mu,sigma,nu,tau) result(v)
   real(dp),intent(in)::mu,sigma,nu,tau
   real(dp)::u
   call random_number(u)
   v=qJSU(u,mu,sigma,nu,tau)
   end function

end module gamlss_continuous
