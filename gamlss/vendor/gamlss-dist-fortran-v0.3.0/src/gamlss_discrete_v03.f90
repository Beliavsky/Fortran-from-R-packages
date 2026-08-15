! Remaining discrete and zero-modified gamlss.dist families for v0.3.0.
! Upstream gamlss.dist GPL-2 | GPL-3. Translation SPDX-License-Identifier: GPL-3.0-only
module gamlss_discrete_v03
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp
   use gamlss_base, only : dbinom_v
   use gamlss_discrete, only : dPIG,pPIG,qPIG,rPIG,dBB,pBB,qBB,rBB,dBNB,pBNB,qBNB,rBNB
   use gamlss_discrete_v02, only : dSICHEL,pSICHEL,qSICHEL,rSICHEL,dZIPF,pZIPF,qZIPF,rZIPF
   implicit none
   private
   public :: dDBI,pDBI,qDBI,rDBI,dPIG2,pPIG2,qPIG2,rPIG2
   public :: dZIPIG,pZIPIG,qZIPIG,rZIPIG,dZAPIG,pZAPIG,qZAPIG,rZAPIG
   public :: dZISICHEL,pZISICHEL,qZISICHEL,rZISICHEL,dZASICHEL,pZASICHEL,qZASICHEL,rZASICHEL
   public :: dZIBB,pZIBB,qZIBB,rZIBB,dZABB,pZABB,qZABB,rZABB
   public :: dZIBNB,pZIBNB,qZIBNB,rZIBNB,dZABNB,pZABNB,qZABNB,rZABNB
   public :: dZAZIPF,pZAZIPF,qZAZIPF,rZAZIPF
contains
   elemental real(dp) function nanv() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nanv
   elemental real(dp) function infv() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function infv
   elemental logical function want_log(flag) result(v)
      logical,intent(in),optional::flag
      v=.false.; if(present(flag))v=flag
   end function want_log
   elemental logical function lower(flag) result(v)
      logical,intent(in),optional::flag
      v=.true.; if(present(flag))v=flag
   end function lower
   elemental logical function is_int(x) result(v)
      real(dp),intent(in)::x
      v=abs(x-real(nint(x),dp))<=16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x))
   end function is_int
   elemental real(dp) function input_prob(p,lower_tail,log_p) result(v)
      real(dp),intent(in)::p
      logical,intent(in),optional::lower_tail,log_p
      v=p
      if(present(log_p))then; if(log_p)v=exp(v); end if
      if(.not.lower(lower_tail))v=1.0_dp-v
   end function input_prob
   elemental real(dp) function finish_prob(p,lower_tail,log_p) result(v)
      real(dp),intent(in)::p
      logical,intent(in),optional::lower_tail,log_p
      v=max(0.0_dp,min(1.0_dp,p))
      if(.not.lower(lower_tail))v=1.0_dp-v
      if(present(log_p))then
         if(log_p)then
            if(v<=0.0_dp)then; v=-infv(); else; v=log(v); end if
         end if
      end if
   end function finish_prob

   real(dp) function dDBI(x,mu,sigma,bd,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,bd
      logical,intent(in),optional::log_density
      integer::k,n,j
      real(dp)::ld,lz,mx,t,s
      if(mu<0.0_dp.or.mu>1.0_dp.or.sigma<=0.0_dp.or.bd<0.0_dp.or..not.is_int(bd))then
         v=nanv(); return
      end if
      n=nint(bd)
      if(x<0.0_dp.or.x>bd.or..not.is_int(x))then
         v=merge(-infv(),0.0_dp,want_log(log_density)); return
      end if
      k=nint(x)
      if(abs(sigma-1.0_dp)<1.0e-3_dp)then
         v=dbinom_v(k,n,mu,log_density); return
      end if
      mx=-huge(1.0_dp)
      do j=0,n
         t=dbi_logweight(j,n,mu,sigma)
         if(t>mx)mx=t
      end do
      s=0.0_dp
      do j=0,n
         s=s+exp(dbi_logweight(j,n,mu,sigma)-mx)
      end do
      lz=mx+log(s); ld=dbi_logweight(k,n,mu,sigma)-lz
      v=merge(ld,exp(ld),want_log(log_density))
   end function dDBI
   pure real(dp) function dbi_logweight(k,n,mu,sigma) result(v)
      integer,intent(in)::k,n
      real(dp),intent(in)::mu,sigma
      real(dp)::rk,rn,rnk,a,b
      rk=real(k,dp); rn=real(n,dp); rnk=real(n-k,dp)
      if((mu<=0.0_dp.and.k>0).or.(mu>=1.0_dp.and.k<n))then; v=-huge(1.0_dp); return; end if
      a=0.0_dp; b=0.0_dp
      if(k>0)a=rk*log(rk)
      if(n-k>0)b=rnk*log(rnk)
      v=log_gamma(rn+1.0_dp)-log_gamma(rk+1.0_dp)-log_gamma(rnk+1.0_dp)
      v=v+(1.0_dp-1.0_dp/sigma)*(a+b-rn*log(max(rn,1.0_dp)))
      v=v+(rn/sigma)*log(max(rn,1.0_dp))
      if(k>0)v=v+(rk/sigma)*log(mu)
      if(n-k>0)v=v+(rnk/sigma)*log(1.0_dp-mu)
      ! Remove duplicated rn/sigma term introduced by the compact expression above.
      v=v-(rn/sigma)*log(max(rn,1.0_dp))
   end function dbi_logweight
   real(dp) function pDBI(q,mu,sigma,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,bd
      logical,intent(in),optional::lower_tail,log_p
      integer::j,k,n
      real(dp)::s
      if(bd<0.0_dp.or..not.is_int(bd))then; v=nanv(); return; end if
      n=nint(bd); k=min(n,int(floor(q))); s=0.0_dp
      if(k>=0)then; do j=0,k; s=s+dDBI(real(j,dp),mu,sigma,bd); end do; end if
      v=finish_prob(s,lower_tail,log_p)
   end function pDBI
   integer function qDBI(p,mu,sigma,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,c
      integer::j,n
      pp=input_prob(p,lower_tail,log_p); n=nint(bd); c=0.0_dp; v=n
      do j=0,n
         c=c+dDBI(real(j,dp),mu,sigma,bd)
         if(c>=pp)then; v=j; return; end if
      end do
   end function qDBI
   integer function rDBI(mu,sigma,bd) result(v)
      real(dp),intent(in)::mu,sigma,bd
      real(dp)::u
      call random_number(u); v=qDBI(u,mu,sigma,bd)
   end function rDBI

   real(dp) function dPIG2(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::alpha
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      alpha=1.0_dp/(sqrt(mu*mu+sigma*sigma)-mu)
      v=dPIG(x,mu,alpha,log_density)
   end function dPIG2
   real(dp) function pPIG2(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::alpha
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      alpha=1.0_dp/(sqrt(mu*mu+sigma*sigma)-mu)
      v=pPIG(q,mu,alpha,lower_tail,log_p)
   end function pPIG2
   integer function qPIG2(p,mu,sigma,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::alpha,pp
      alpha=1.0_dp/(sqrt(mu*mu+sigma*sigma)-mu); pp=input_prob(p,lower_tail,log_p)
      v=qPIG(pp,mu,alpha,max_value=max_value)
   end function qPIG2
   integer function rPIG2(mu,sigma,max_value) result(v)
      real(dp),intent(in)::mu,sigma
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qPIG2(u,mu,sigma,max_value=max_value)
   end function rPIG2

   real(dp) function dZIPIG(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      v=zi_pmf(x,nu,dPIG(x,mu,sigma),dPIG(0.0_dp,mu,sigma),log_density)
   end function dZIPIG
   real(dp) function pZIPIG(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(q<0.0_dp)then; raw=0.0_dp; else; raw=nu+(1.0_dp-nu)*pPIG(q,mu,sigma); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZIPIG
   integer function qZIPIG(p,mu,sigma,nu,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,pb
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=nu)then; v=0; else; pb=(pp-nu)/(1.0_dp-nu); v=qPIG(pb,mu,sigma,max_value=max_value); end if
   end function qZIPIG
   integer function rZIPIG(mu,sigma,nu,max_value) result(v)
      real(dp),intent(in)::mu,sigma,nu
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZIPIG(u,mu,sigma,nu,max_value=max_value)
   end function rZIPIG

   real(dp) function dZAPIG(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      v=za_pmf(x,nu,dPIG(x,mu,sigma),dPIG(0.0_dp,mu,sigma),log_density)
   end function dZAPIG
   real(dp) function pZAPIG(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::f0,raw
      f0=dPIG(0.0_dp,mu,sigma)
      if(q<0.0_dp)then; raw=0.0_dp; else if(q<1.0_dp)then; raw=nu; &
      else; raw=nu+(1.0_dp-nu)*(pPIG(q,mu,sigma)-f0)/(1.0_dp-f0); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZAPIG
   integer function qZAPIG(p,mu,sigma,nu,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,f0,pb
      pp=input_prob(p,lower_tail,log_p); f0=dPIG(0.0_dp,mu,sigma)
      if(pp<=nu)then; v=0; else; pb=f0+(1.0_dp-f0)*(pp-nu)/(1.0_dp-nu); &
      v=qPIG(pb,mu,sigma,max_value=max_value); end if
   end function qZAPIG
   integer function rZAPIG(mu,sigma,nu,max_value) result(v)
      real(dp),intent(in)::mu,sigma,nu
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZAPIG(u,mu,sigma,nu,max_value=max_value)
   end function rZAPIG

   real(dp) function dZISICHEL(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      v=zi_pmf(x,tau,dSICHEL(x,mu,sigma,nu),dSICHEL(0.0_dp,mu,sigma,nu),log_density)
   end function dZISICHEL
   real(dp) function pZISICHEL(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(q<0.0_dp)then; raw=0.0_dp; else; raw=tau+(1.0_dp-tau)*pSICHEL(q,mu,sigma,nu); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZISICHEL
   integer function qZISICHEL(p,mu,sigma,nu,tau,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,pb
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=tau)then; v=0; else; pb=(pp-tau)/(1.0_dp-tau); &
      v=qSICHEL(pb,mu,sigma,nu,max_value=max_value); end if
   end function qZISICHEL
   integer function rZISICHEL(mu,sigma,nu,tau,max_value) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZISICHEL(u,mu,sigma,nu,tau,max_value=max_value)
   end function rZISICHEL

   real(dp) function dZASICHEL(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      v=za_pmf(x,tau,dSICHEL(x,mu,sigma,nu),dSICHEL(0.0_dp,mu,sigma,nu),log_density)
   end function dZASICHEL
   real(dp) function pZASICHEL(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::f0,raw
      f0=dSICHEL(0.0_dp,mu,sigma,nu)
      if(q<0.0_dp)then; raw=0.0_dp; else if(q<1.0_dp)then; raw=tau; &
      else; raw=tau+(1.0_dp-tau)*(pSICHEL(q,mu,sigma,nu)-f0)/(1.0_dp-f0); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZASICHEL
   integer function qZASICHEL(p,mu,sigma,nu,tau,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,f0,pb
      pp=input_prob(p,lower_tail,log_p); f0=dSICHEL(0.0_dp,mu,sigma,nu)
      if(pp<=tau)then; v=0; else; pb=f0+(1.0_dp-f0)*(pp-tau)/(1.0_dp-tau); &
      v=qSICHEL(pb,mu,sigma,nu,max_value=max_value); end if
   end function qZASICHEL
   integer function rZASICHEL(mu,sigma,nu,tau,max_value) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZASICHEL(u,mu,sigma,nu,tau,max_value=max_value)
   end function rZASICHEL

   real(dp) function dZIBB(x,mu,sigma,nu,bd,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      integer,intent(in)::bd
      logical,intent(in),optional::log_density
      v=zi_pmf(x,nu,dBB(x,mu,sigma,bd),dBB(0.0_dp,mu,sigma,bd),log_density)
   end function dZIBB
   real(dp) function pZIBB(q,mu,sigma,nu,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(q<0.0_dp)then; raw=0.0_dp; else; raw=nu+(1.0_dp-nu)*pBB(q,mu,sigma,bd); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZIBB
   integer function qZIBB(p,mu,sigma,nu,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,pb
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=nu)then; v=0; else; pb=(pp-nu)/(1.0_dp-nu); v=qBB(pb,mu,sigma,bd); end if
   end function qZIBB
   integer function rZIBB(mu,sigma,nu,bd) result(v)
      real(dp),intent(in)::mu,sigma,nu
      integer,intent(in)::bd
      real(dp)::u
      call random_number(u); v=qZIBB(u,mu,sigma,nu,bd)
   end function rZIBB

   real(dp) function dZABB(x,mu,sigma,nu,bd,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      integer,intent(in)::bd
      logical,intent(in),optional::log_density
      v=za_pmf(x,nu,dBB(x,mu,sigma,bd),dBB(0.0_dp,mu,sigma,bd),log_density)
   end function dZABB
   real(dp) function pZABB(q,mu,sigma,nu,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::f0,raw
      f0=dBB(0.0_dp,mu,sigma,bd)
      if(q<0.0_dp)then; raw=0.0_dp; else if(q<1.0_dp)then; raw=nu; &
      else; raw=nu+(1.0_dp-nu)*(pBB(q,mu,sigma,bd)-f0)/(1.0_dp-f0); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZABB
   integer function qZABB(p,mu,sigma,nu,bd,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      integer,intent(in)::bd
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,f0,pb
      pp=input_prob(p,lower_tail,log_p); f0=dBB(0.0_dp,mu,sigma,bd)
      if(pp<=nu)then; v=0; else; pb=f0+(1.0_dp-f0)*(pp-nu)/(1.0_dp-nu); v=qBB(pb,mu,sigma,bd); end if
   end function qZABB
   integer function rZABB(mu,sigma,nu,bd) result(v)
      real(dp),intent(in)::mu,sigma,nu
      integer,intent(in)::bd
      real(dp)::u
      call random_number(u); v=qZABB(u,mu,sigma,nu,bd)
   end function rZABB

   real(dp) function dZIBNB(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      v=zi_pmf(x,tau,dBNB(x,mu,sigma,nu),dBNB(0.0_dp,mu,sigma,nu),log_density)
   end function dZIBNB
   real(dp) function pZIBNB(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(q<0.0_dp)then; raw=0.0_dp; else; raw=tau+(1.0_dp-tau)*pBNB(q,mu,sigma,nu); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZIBNB
   integer function qZIBNB(p,mu,sigma,nu,tau,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,pb
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=tau)then; v=0; else; pb=(pp-tau)/(1.0_dp-tau); v=qBNB(pb,mu,sigma,nu,max_value=max_value); end if
   end function qZIBNB
   integer function rZIBNB(mu,sigma,nu,tau,max_value) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZIBNB(u,mu,sigma,nu,tau,max_value=max_value)
   end function rZIBNB

   real(dp) function dZABNB(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      v=za_pmf(x,tau,dBNB(x,mu,sigma,nu),dBNB(0.0_dp,mu,sigma,nu),log_density)
   end function dZABNB
   real(dp) function pZABNB(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::f0,raw
      f0=dBNB(0.0_dp,mu,sigma,nu)
      if(q<0.0_dp)then; raw=0.0_dp; else if(q<1.0_dp)then; raw=tau; &
      else; raw=tau+(1.0_dp-tau)*(pBNB(q,mu,sigma,nu)-f0)/(1.0_dp-f0); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZABNB
   integer function qZABNB(p,mu,sigma,nu,tau,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,f0,pb
      pp=input_prob(p,lower_tail,log_p); f0=dBNB(0.0_dp,mu,sigma,nu)
      if(pp<=tau)then; v=0; else; pb=f0+(1.0_dp-f0)*(pp-tau)/(1.0_dp-tau); &
      v=qBNB(pb,mu,sigma,nu,max_value=max_value); end if
   end function qZABNB
   integer function rZABNB(mu,sigma,nu,tau,max_value) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZABNB(u,mu,sigma,nu,tau,max_value=max_value)
   end function rZABNB

   real(dp) function dZAZIPF(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::p,ld
      if(sigma<=0.0_dp.or.sigma>=1.0_dp)then; v=nanv(); return; end if
      if(x<0.0_dp.or.(x>0.0_dp.and..not.is_int(x)))then
         v=merge(-infv(),0.0_dp,want_log(log_density)); return
      end if
      if(x==0.0_dp)then; p=sigma; else; p=(1.0_dp-sigma)*dZIPF(x,mu); end if
      ld=merge(log(p),-infv(),p>0.0_dp); v=merge(ld,p,want_log(log_density))
   end function dZAZIPF
   real(dp) function pZAZIPF(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(q<0.0_dp)then; raw=0.0_dp; else if(q<1.0_dp)then; raw=sigma; &
      else; raw=sigma+(1.0_dp-sigma)*pZIPF(q,mu); end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pZAZIPF
   integer function qZAZIPF(p,mu,sigma,lower_tail,log_p,max_value) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::max_value
      real(dp)::pp,pb
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=sigma)then; v=0; else; pb=(pp-sigma)/(1.0_dp-sigma); v=qZIPF(pb,mu,max_value=max_value); end if
   end function qZAZIPF
   integer function rZAZIPF(mu,sigma,max_value) result(v)
      real(dp),intent(in)::mu,sigma
      integer,intent(in),optional::max_value
      real(dp)::u
      call random_number(u); v=qZAZIPF(u,mu,sigma,max_value=max_value)
   end function rZAZIPF

   real(dp) function zi_pmf(x,pzero,pbase,p0,log_density) result(v)
      real(dp),intent(in)::x,pzero,pbase,p0
      logical,intent(in),optional::log_density
      real(dp)::p
      if(pzero<=0.0_dp.or.pzero>=1.0_dp)then; v=nanv(); return; end if
      if(x<0.0_dp.or..not.is_int(x))then; p=0.0_dp; &
      else if(nint(x)==0)then; p=pzero+(1.0_dp-pzero)*p0; else; p=(1.0_dp-pzero)*pbase; end if
      if(want_log(log_density))then; if(p>0.0_dp)then; v=log(p); else; v=-infv(); end if; else; v=p; end if
   end function zi_pmf
   real(dp) function za_pmf(x,pzero,pbase,p0,log_density) result(v)
      real(dp),intent(in)::x,pzero,pbase,p0
      logical,intent(in),optional::log_density
      real(dp)::p
      if(pzero<=0.0_dp.or.pzero>=1.0_dp)then; v=nanv(); return; end if
      if(x<0.0_dp.or..not.is_int(x))then; p=0.0_dp; &
      else if(nint(x)==0)then; p=pzero; else; p=(1.0_dp-pzero)*pbase/(1.0_dp-p0); end if
      if(want_log(log_density))then; if(p>0.0_dp)then; v=log(p); else; v=-infv(); end if; else; v=p; end if
   end function za_pmf
end module gamlss_discrete_v03
