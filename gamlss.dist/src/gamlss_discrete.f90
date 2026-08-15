! Computational translation of discrete families from gamlss.dist 6.1-1.
! Original package GPL-2 | GPL-3. Translation SPDX-License-Identifier: GPL-3.0-only
module gamlss_discrete
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp
   use gamlss_special, only : log_beta_fn, log1p_v
   use gamlss_base, only : dpois_v, ppois_v, qpois_v, rpois_v, dbinom_v, pbinom_v, qbinom_v, rbinom_v, &
      dnbinom_v, pnbinom_v, qnbinom_v, rnbinom_v
   implicit none
   private
   public :: dPO,pPO,qPO,rPO,dBI,pBI,qBI,rBI,dGEOM,pGEOM,qGEOM,rGEOM,dGEOMo,pGEOMo,qGEOMo,rGEOMo
   public :: dNBI,pNBI,qNBI,rNBI,dNBII,pNBII,qNBII,rNBII
   public :: dZIP,pZIP,qZIP,rZIP,dZIP2,pZIP2,qZIP2,rZIP2,dZAP,pZAP,qZAP,rZAP
   public :: dZINBI,pZINBI,qZINBI,rZINBI,dZANBI,pZANBI,qZANBI,rZANBI
   public :: dZIBI,pZIBI,qZIBI,rZIBI,dZABI,pZABI,qZABI,rZABI
   public :: dBB,pBB,qBB,rBB,dBNB,pBNB,qBNB,rBNB
   public :: dPIG,pPIG,qPIG,rPIG
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
      if(.not.lower(lower_tail))v=1-v
      v=max(0.0_dp,min(1.0_dp,v))
      if(present(log_p))then
      if(log_p)then
      if(v==0)then
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
      if(.not.lower(lower_tail))v=1-v
   end function
   elemental logical function is_int(x) result(ok)
   real(dp),intent(in)::x
   ok=abs(x-nint(x))<1e-10_dp
   end function

   elemental real(dp) function dPO(x,mu,log_density) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_density
      if(mu<=0)then
      v=nanv()
      else if(x<0.0_dp .or. .not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else
      v=dpois_v(nint(x),mu,want_log(log_density))
      end if
   end function
   elemental real(dp) function pPO(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0)then
      v=nanv()
      else
      v=finish_prob(ppois_v(int(floor(q)),mu),lower_tail,log_p)
      end if
   end function
   integer function qPO(p,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      if(mu<=0)then
      v=-huge(v)
      else
      v=qpois_v(input_prob(p,lower_tail,log_p),mu)
      end if
   end function
   integer function rPO(mu) result(v)
   real(dp),intent(in)::mu
   v=rpois_v(mu)
   end function

   elemental real(dp) function dBI(x,bd,mu,log_density) result(v)
      real(dp),intent(in)::x,mu
      integer,intent(in)::bd
      logical,intent(in),optional::log_density
      if(mu<0.or.mu>1.or.bd<0)then
      v=nanv()
      else if(x<0.0_dp .or. x>real(bd,dp) .or. .not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else
      v=dbinom_v(nint(x),bd,mu,want_log(log_density))
      end if
   end function
   elemental real(dp) function pBI(q,bd,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      if(mu<0.or.mu>1.or.bd<0)then
      v=nanv()
      else
      v=finish_prob(pbinom_v(int(floor(q)),bd,mu),lower_tail,log_p)
      end if
   end function
   integer function qBI(p,bd,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      v=qbinom_v(input_prob(p,lower_tail,log_p),bd,mu)
   end function
   integer function rBI(bd,mu) result(v)
   integer,intent(in)::bd
   real(dp),intent(in)::mu
   v=rbinom_v(bd,mu)
   end function

   elemental real(dp) function dGEOMo(x,mu,log_density) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(mu<=0.or.mu>1)then
      v=nanv()
      else if(x<0.or..not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         ld=log(mu)+x*log1p_v(-mu)
         v=merge(ld,exp(ld),want_log(log_density))
         end if
   end function
   elemental real(dp) function pGEOMo(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(q<0)then
      p=0
      else
      p=1-(1-mu)**(floor(q)+1)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qGEOMo(p,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u
      u=input_prob(p,lower_tail,log_p)
      if(u<=0)then
      v=0
      else if(u>=1)then
      v=huge(v)
      else
      v=max(0,ceiling(log1p_v(-u)/log1p_v(-mu)-1.0_dp))
      end if
   end function
   integer function rGEOMo(mu) result(v)
   real(dp),intent(in)::mu
   real(dp)::u
   call random_number(u)
   v=qGEOMo(u,mu)
   end function
   elemental real(dp) function dGEOM(x,mu,log_density) result(v)
   real(dp),intent(in)::x,mu
   logical,intent(in),optional::log_density
      if(mu<0)then
      v=nanv()
      else
      v=dGEOMo(x,1/(mu+1),log_density)
      end if
      end function
   elemental real(dp) function pGEOM(q,mu,lower_tail,log_p) result(v)
   real(dp),intent(in)::q,mu
   logical,intent(in),optional::lower_tail,log_p
      if(mu<0)then
      v=nanv()
      else
      v=pGEOMo(q,1/(mu+1),lower_tail,log_p)
      end if
      end function
   integer function qGEOM(p,mu,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu
   logical,intent(in),optional::lower_tail,log_p
   v=qGEOMo(p,1/(mu+1),lower_tail,log_p)
   end function
   integer function rGEOM(mu) result(v)
   real(dp),intent(in)::mu
   v=rGEOMo(1/(mu+1))
   end function

   elemental subroutine nbpars(mu,sigma,type2,size,prob,ok)
      real(dp),intent(in)::mu,sigma
      logical,intent(in)::type2
      real(dp),intent(out)::size,prob
      logical,intent(out)::ok
      ok=mu>0.and.sigma>0
      if(.not.ok)then
      size=1
      prob=.5
      return
      end if
      if(type2)then
      size=mu/sigma
      else
      size=1/sigma
      end if
      prob=size/(size+mu)
   end subroutine
   elemental real(dp) function dNBI(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::size,prob
      logical::ok
      call nbpars(mu,sigma,.false.,size,prob,ok)
      if(.not.ok)then
      v=nanv()
      else if(x<0.0_dp .or. .not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else if(sigma<1e-4_dp)then
      v=dPO(x,mu,log_density)
      else
      v=dnbinom_v(nint(x),size,prob,want_log(log_density))
      end if
   end function
   elemental real(dp) function pNBI(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::size,prob,p
      logical::ok
      call nbpars(mu,sigma,.false.,size,prob,ok)
      if(.not.ok)then
      v=nanv()
      else
      if(sigma<1e-4_dp)then
      p=ppois_v(int(floor(q)),mu)
      else
      p=pnbinom_v(int(floor(q)),size,prob)
      end if
      v=finish_prob(p,lower_tail,log_p)
      end if
   end function
   integer function qNBI(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::size,prob,u
      logical::ok
      call nbpars(mu,sigma,.false.,size,prob,ok)
      u=input_prob(p,lower_tail,log_p)
      if(sigma<1e-4_dp)then
      v=qpois_v(u,mu)
      else
      v=qnbinom_v(u,size,prob)
      end if
   end function
   integer function rNBI(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::size,prob
   logical::ok
      call nbpars(mu,sigma,.false.,size,prob,ok)
      if(sigma<1e-4_dp)then
      v=rpois_v(mu)
      else
      v=rnbinom_v(size,prob)
      end if
      end function
   elemental real(dp) function dNBII(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::size,prob
      logical::ok
      call nbpars(mu,sigma,.true.,size,prob,ok)
      if(.not.ok)then
      v=nanv()
      else if(x<0.0_dp .or. .not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      else if(sigma<1e-4_dp)then
      v=dPO(x,mu,log_density)
      else
      v=dnbinom_v(nint(x),size,prob,want_log(log_density))
      end if
   end function
   elemental real(dp) function pNBII(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::size,prob,p
      logical::ok
      call nbpars(mu,sigma,.true.,size,prob,ok)
      if(.not.ok)then
      v=nanv()
      else
      if(sigma<1e-4_dp)then
      p=ppois_v(int(floor(q)),mu)
      else
      p=pnbinom_v(int(floor(q)),size,prob)
      end if
      v=finish_prob(p,lower_tail,log_p)
      end if
   end function
   integer function qNBII(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::size,prob,u
      logical::ok
      call nbpars(mu,sigma,.true.,size,prob,ok)
      u=input_prob(p,lower_tail,log_p)
      if(sigma<1e-4_dp)then
      v=qpois_v(u,mu)
      else
      v=qnbinom_v(u,size,prob)
      end if
   end function
   integer function rNBII(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::size,prob
   logical::ok
      call nbpars(mu,sigma,.true.,size,prob,ok)
      if(sigma<1e-4_dp)then
      v=rpois_v(mu)
      else
      v=rnbinom_v(size,prob)
      end if
      end function

! Zero-inflated Poisson and mean-parameterized ZIP2.
   elemental real(dp) function dZIP(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::p,ld
      if(mu<=0.or.sigma<=0.or.sigma>=1)then
      v=nanv()
      return
      end if
      if(x<0.or..not.is_int(x))then
      ld=-infv()
      else if(x==0)then
      p=sigma+(1-sigma)*exp(-mu)
      ld=log(p)
      else
      ld=log(1-sigma)+dPO(x,mu,.true.)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pZIP(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(q<0)then
      p=0
      else
      p=sigma+(1-sigma)*ppois_v(int(floor(q)),mu)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qZIP(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,pnew
      u=input_prob(p,lower_tail,log_p)
      pnew=(u-sigma)/(1-sigma)
      if(pnew<=0)then
      v=0
      else
      v=qpois_v(pnew,mu)
      end if
   end function
   integer function rZIP(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   if(u<sigma)then
   v=0
   else
   v=rpois_v(mu)
   end if
   end function
   elemental real(dp) function dZIP2(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::mus
      if(sigma<=0.or.sigma>=1)then
      v=nanv()
      else
      mus=mu/(1-sigma)
      v=dZIP(x,mus,sigma,log_density)
      end if
   end function
   elemental real(dp) function pZIP2(q,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::q,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
   v=pZIP(q,mu/(1-sigma),sigma,lower_tail,log_p)
   end function
   integer function qZIP2(p,mu,sigma,lower_tail,log_p) result(v)
   real(dp),intent(in)::p,mu,sigma
   logical,intent(in),optional::lower_tail,log_p
   v=qZIP(p,mu/(1-sigma),sigma,lower_tail,log_p)
   end function
   integer function rZIP2(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   v=rZIP(mu/(1-sigma),sigma)
   end function

! Generic zero-altered base helpers, specialized below.
   elemental real(dp) function dZAP(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::p0,ld
      if(mu<=0.or.sigma<=0.or.sigma>=1)then
      v=nanv()
      return
      end if
      if(x<0.or..not.is_int(x))then
      ld=-infv()
      else if(x==0)then
      ld=log(sigma)
      else
      p0=exp(-mu)
      ld=log(1-sigma)+dPO(x,mu,.true.)-log1p_v(-p0)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pZAP(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p0,p
      if(q<0)then
      p=0
      else if(q<1)then
      p=sigma
      else
      p0=exp(-mu)
      p=sigma+(1-sigma)*(ppois_v(int(floor(q)),mu)-p0)/(1-p0)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qZAP(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,pnew,p0,target
      u=input_prob(p,lower_tail,log_p)
      if(u<=sigma)then
      v=0
      else
      pnew=(u-sigma)/(1-sigma)
      p0=exp(-mu)
      target=p0+(1-p0)*pnew
      v=max(1,qpois_v(target,mu))
      end if
   end function
   integer function rZAP(mu,sigma) result(v)
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   v=qZAP(u,mu,sigma)
   end function

   elemental real(dp) function dZINBI(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::base,ld
      if(nu<=0.or.nu>=1)then
      v=nanv()
      return
      end if
      base=dNBI(x,mu,sigma,.false.)
      if(x<0)then
      ld=-infv()
      else if(x==0)then
      ld=log(nu+(1-nu)*base)
      else
      ld=log(1-nu)+log(base)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pZINBI(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(q<0)then
      p=0
      else
      p=nu+(1-nu)*pNBI(q,mu,sigma)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qZINBI(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,pnew
      u=input_prob(p,lower_tail,log_p)
      pnew=(u-nu)/(1-nu)
      if(pnew<=0)then
      v=0
      else
      v=qNBI(pnew,mu,sigma)
      end if
   end function
   integer function rZINBI(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   real(dp)::u
   call random_number(u)
   v=qZINBI(u,mu,sigma,nu)
   end function
   elemental real(dp) function dZANBI(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::p0,ld
      if(nu<=0.or.nu>=1)then
      v=nanv()
      return
      end if
      if(x<0)then
      ld=-infv()
      else if(x==0)then
      ld=log(nu)
      else
      p0=dNBI(0.0_dp,mu,sigma)
      ld=log(1-nu)+dNBI(x,mu,sigma,.true.)-log1p_v(-p0)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pZANBI(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p0,p
      if(q<0)then
      p=0
      else if(q<1)then
      p=nu
      else
      p0=pNBI(0.0_dp,mu,sigma)
      p=nu+(1-nu)*(pNBI(q,mu,sigma)-p0)/(1-p0)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qZANBI(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,p0,t
      u=input_prob(p,lower_tail,log_p)
      if(u<=nu)then
      v=0
      else
      p0=pNBI(0.0_dp,mu,sigma)
      t=p0+(1-p0)*(u-nu)/(1-nu)
      v=max(1,qNBI(t,mu,sigma))
      end if
   end function
   integer function rZANBI(mu,sigma,nu) result(v)
   real(dp),intent(in)::mu,sigma,nu
   real(dp)::u
   call random_number(u)
   v=qZANBI(u,mu,sigma,nu)
   end function

   elemental real(dp) function dZIBI(x,bd,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::log_density
      real(dp)::base,ld
      base=dBI(x,bd,mu)
      if(x<0)then
      ld=-infv()
      else if(x==0)then
      ld=log(sigma+(1-sigma)*base)
      else
      ld=log(1-sigma)+log(base)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pZIBI(q,bd,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(q<0)then
      p=0
      else
      p=sigma+(1-sigma)*pBI(q,bd,mu)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qZIBI(p,bd,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,t
      u=input_prob(p,lower_tail,log_p)
      t=(u-sigma)/(1-sigma)
      if(t<=0)then
      v=0
      else
      v=qBI(t,bd,mu)
      end if
   end function
   integer function rZIBI(bd,mu,sigma) result(v)
   integer,intent(in)::bd
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   v=qZIBI(u,bd,mu,sigma)
   end function
   elemental real(dp) function dZABI(x,bd,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::log_density
      real(dp)::p0,ld
      if(x<0)then
      ld=-infv()
      else if(x==0)then
      ld=log(sigma)
      else
      p0=dBI(0.0_dp,bd,mu)
      ld=log(1-sigma)+dBI(x,bd,mu,.true.)-log1p_v(-p0)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   elemental real(dp) function pZABI(q,bd,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p0,p
      if(q<0)then
      p=0
      else if(q<1)then
      p=sigma
      else
      p0=pBI(0.0_dp,bd,mu)
      p=sigma+(1-sigma)*(pBI(q,bd,mu)-p0)/(1-p0)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qZABI(p,bd,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,p0,t
      u=input_prob(p,lower_tail,log_p)
      if(u<=sigma)then
      v=0
      else
      p0=pBI(0.0_dp,bd,mu)
      t=p0+(1-p0)*(u-sigma)/(1-sigma)
      v=max(1,qBI(t,bd,mu))
      end if
   end function
   integer function rZABI(bd,mu,sigma) result(v)
   integer,intent(in)::bd
   real(dp),intent(in)::mu,sigma
   real(dp)::u
   call random_number(u)
   v=qZABI(u,bd,mu,sigma)
   end function

! Beta-binomial BB.
   elemental real(dp) function dBB(x,mu,sigma,bd,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::log_density
      real(dp)::a,b,ld,s
      if(mu<0.or.mu>1.or.sigma<=0.or.bd<0.or.x<0.or.x>bd.or..not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      s=max(sigma,1e-10_dp)
      if(s<1e-4_dp)then
      v=dBI(x,bd,mu,log_density)
      return
      end if
      a=mu/s
      b=(1-mu)/s
      ld=log_gamma(real(bd+1,dp))-log_gamma(x+1)-log_gamma(real(bd,dp)-x+1)+log_beta_fn(x+a,real(bd,dp)-x+b)-log_beta_fn(a,b)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   real(dp) function pBB(q,mu,sigma,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      integer::k,qq
      real(dp)::p
      if(q<0)then
      p=0
      else
      qq=min(bd,int(floor(q)))
      p=0
      do k=0,qq
      p=p+dBB(real(k,dp),mu,sigma,bd)
      end do
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function
   integer function qBB(p,mu,sigma,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::u,c
      integer::k
      u=input_prob(p,lower_tail,log_p)
      c=0
      v=bd
      do k=0,bd
      c=c+dBB(real(k,dp),mu,sigma,bd)
      if(u<=c)then
      v=k
      return
      end if
      end do
   end function
   integer function rBB(mu,sigma,bd) result(v)
   real(dp),intent(in)::mu,sigma
   integer,intent(in)::bd
   real(dp)::u
   call random_number(u)
   v=qBB(u,mu,sigma,bd)
   end function

! Beta negative binomial BNB.
   elemental real(dp) function dBNB(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::m,n,k,ld
      if(mu<=0.or.sigma<=0.or.nu<=0.or.x<0.or..not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      m=1/sigma+1
      n=mu*nu/sigma
      k=1/nu
      ld=log_beta_fn(x+n,m+k)-log_beta_fn(n,m)-log_gamma(x+1)-log_gamma(k)+log_gamma(x+k)
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   real(dp) function pBNB(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      integer::j
      real(dp)::p
      if(q<0)then
      p=0
      else
      p=0
      do j=0,int(floor(q))
      p=p+dBNB(real(j,dp),mu,sigma,nu)
      end do
      end if
      v=finish_prob(min(1.0_dp,p),lower_tail,log_p)
   end function
   integer function qBNB(p,mu,sigma,nu,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::u,c
      integer::j,mx
      mx=10000
      if(present(max_value))mx=max_value
      u=input_prob(p,lower_tail,log_p)
      if(u>=1)then
      v=huge(v)
      return
      end if
      c=0
      v=mx
      do j=0,mx
      c=c+dBNB(real(j,dp),mu,sigma,nu)
      if(u<=c)then
      v=j
      return
      end if
      end do
   end function
   integer function rBNB(mu,sigma,nu,max_value) result(v)
   real(dp),intent(in)::mu,sigma,nu
   integer,intent(in),optional::max_value
   real(dp)::u
   call random_number(u)
   v=qBNB(u,mu,sigma,nu,max_value=max_value)
   end function

! Poisson-inverse-Gaussian recurrence ported from tofyPIG2.c/tocdf.c.
   elemental real(dp) function dPIG(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::t0,tprev,tcur,sumlog,ld
      integer::j,y
      if(mu<=0.or.sigma<=0)then
      v=nanv()
      return
      end if
      if(x<0.or..not.is_int(x))then
      v=merge(-infv(),0.0_dp,want_log(log_density))
      return
      end if
      y=nint(x)
      t0=mu/sqrt(1+2*sigma*mu)
      tprev=t0
      sumlog=0
      do j=1,y
         tcur=(sigma*(2*j-1)/mu+1/tprev)*t0*t0
         sumlog=sumlog+log(tprev)
         tprev=tcur
      end do
      ld=-log_gamma(x+1)+(1-sqrt(1+2*sigma*mu))/sigma+sumlog
      v=merge(ld,exp(ld),want_log(log_density))
   end function
   real(dp) function pPIG(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      integer::j
      if(q<0)then
      p=0
      else
      p=0
      do j=0,int(floor(q))
      p=p+dPIG(real(j,dp),mu,sigma)
      end do
      end if
      v=finish_prob(min(1.0_dp,p),lower_tail,log_p)
   end function
   integer function qPIG(p,mu,sigma,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::u,c
      integer::j,mx
      mx=10000
      if(present(max_value))mx=max_value
      u=input_prob(p,lower_tail,log_p)
      c=0
      v=mx
      do j=0,mx
      c=c+dPIG(real(j,dp),mu,sigma)
      if(u<=c)then
      v=j
      return
      end if
      end do
   end function
   integer function rPIG(mu,sigma,max_value) result(v)
   real(dp),intent(in)::mu,sigma
   integer,intent(in),optional::max_value
   real(dp)::u
   call random_number(u)
   v=qPIG(u,mu,sigma,max_value=max_value)
   end function

end module gamlss_discrete
