! Remaining continuous gamlss.dist families for v0.3.0.
! Upstream gamlss.dist GPL-2 | GPL-3. Translation SPDX-License-Identifier: GPL-3.0-only
module gamlss_continuous_v03
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp, pi
   use gamlss_base, only : dnorm_v, pnorm_v, qnorm_v, qbeta_v
   use gamlss_special, only : normal_cdf, normal_quantile, regularized_beta
   use gamlss_student_t, only : student_t_pdf, student_t_cdf, student_t_quantile
   use gamlss_v02_numerics, only : adaptive_integral, bisection_root
   use gamlss_continuous_v02, only : dST3,pST3,qST3,rST3
   implicit none
   private
   public :: dST3C,pST3C,qST3C,rST3C,dSN1,pSN1,qSN1,rSN1,dSN2,pSN2,qSN2,rSN2
   public :: dSST,pSST,qSST,rSST,dGT,pGT,qGT,rGT,dexGAUS,pexGAUS,qexGAUS,rexGAUS
   public :: dPARETO,pPARETO,qPARETO,rPARETO,dPARETO1,pPARETO1,qPARETO1,rPARETO1
   public :: dPARETO2,pPARETO2,qPARETO2,rPARETO2,dPARETO2o,pPARETO2o,qPARETO2o,rPARETO2o
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

   elemental real(dp) function dST3C(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      v=dST3(x,mu,sigma,nu,tau,log_density)
   end function dST3C
   elemental real(dp) function pST3C(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      v=pST3(q,mu,sigma,nu,tau,lower_tail,log_p)
   end function pST3C
   real(dp) function qST3C(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      v=qST3(p,mu,sigma,nu,tau,lower_tail,log_p)
   end function qST3C
   real(dp) function rST3C(mu,sigma,nu,tau) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      v=rST3(mu,sigma,nu,tau)
   end function rST3C

   elemental real(dp) function dSN1(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::z,ld,cdf
      if(sigma<=0.0_dp)then; v=nanv(); return; end if
      z=(x-mu)/sigma
      cdf=normal_cdf(nu*z)
      if(cdf<=0.0_dp)then
         ld=-infv()
      else
         ld=log(2.0_dp)-log(sigma)-0.5_dp*log(2.0_dp*pi)-0.5_dp*z*z+log(cdf)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSN1
   real(dp) function pSN1(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,raw,a
      if(sigma<=0.0_dp)then; v=nanv(); return; end if
      z=(q-mu)/sigma
      if(z<=-12.0_dp)then
         raw=0.0_dp
      else if(z>=12.0_dp)then
         raw=1.0_dp
      else
         a=-12.0_dp
         raw=adaptive_integral(std_pdf,a,z,5.0e-10_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   contains
      function std_pdf(t) result(y)
         real(dp),intent(in)::t
         real(dp)::y
         y=dSN1(t,0.0_dp,1.0_dp,nu)
      end function std_pdf
   end function pSN1
   real(dp) function qSN1(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,lo,hi
      pp=input_prob(p,lower_tail,log_p)
      if(pp<0.0_dp.or.pp>1.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      if(pp==0.0_dp)then; v=-infv(); return; end if
      if(pp==1.0_dp)then; v=infv(); return; end if
      lo=mu-8.0_dp*sigma; hi=mu+8.0_dp*sigma
      do while(pSN1(lo,mu,sigma,nu)>pp); lo=mu-2.0_dp*(mu-lo); end do
      do while(pSN1(hi,mu,sigma,nu)<pp); hi=mu+2.0_dp*(hi-mu); end do
      v=bisection_root(cdf,pp,lo,hi,2.0e-9_dp,120)
   contains
      function cdf(x) result(y)
         real(dp),intent(in)::x
         real(dp)::y
         y=pSN1(x,mu,sigma,nu)
      end function cdf
   end function qSN1
   real(dp) function rSN1(mu,sigma,nu) result(v)
      real(dp),intent(in)::mu,sigma,nu
      real(dp)::u
      call random_number(u); v=qSN1(u,mu,sigma,nu)
   end function rSN1

   elemental real(dp) function dSN2(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::z,ld
      if(sigma<=0.0_dp.or.nu<=0.0_dp)then; v=nanv(); return; end if
      z=(x-mu)/sigma
      if(x<mu)then; ld=-0.5_dp*(nu*z)**2; else; ld=-0.5_dp*(z/nu)**2; end if
      ld=ld-log(sigma)+log(nu)-log(1.0_dp+nu*nu)+0.5_dp*log(2.0_dp/pi)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSN2
   elemental real(dp) function pSN2(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,k,raw
      if(sigma<=0.0_dp.or.nu<=0.0_dp)then; v=nanv(); return; end if
      z=(q-mu)/sigma; k=nu*nu
      if(q<mu)then
         raw=2.0_dp*normal_cdf(nu*z)/(1.0_dp+k)
      else
         raw=(1.0_dp+k*(2.0_dp*normal_cdf(z/nu)-1.0_dp))/(1.0_dp+k)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pSN2
   elemental real(dp) function qSN2(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,k
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=0.0_dp)then; v=-infv(); return; end if
      if(pp>=1.0_dp)then; v=infv(); return; end if
      if(sigma<=0.0_dp.or.nu<=0.0_dp)then; v=nanv(); return; end if
      k=nu*nu
      if(pp<1.0_dp/(1.0_dp+k))then
         v=mu+(sigma/nu)*normal_quantile(0.5_dp*pp*(1.0_dp+k))
      else
         v=mu+sigma*nu*normal_quantile(0.5_dp*(1.0_dp+(pp*(1.0_dp+k)-1.0_dp)/k))
      end if
   end function qSN2
   real(dp) function rSN2(mu,sigma,nu) result(v)
      real(dp),intent(in)::mu,sigma,nu
      real(dp)::u
      call random_number(u); v=qSN2(u,mu,sigma,nu)
   end function rSN2

   elemental real(dp) function dSST(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::m1,m2,s1,mu1,sigma1
      if(sigma<=0.0_dp.or.nu<=0.0_dp.or.tau<=2.0_dp)then; v=nanv(); return; end if
      call sst_transform(mu,sigma,nu,tau,m1,m2,s1,mu1,sigma1)
      v=dST3(x,mu1,sigma1,nu,tau,log_density)
   end function dSST
   elemental real(dp) function pSST(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::m1,m2,s1,mu1,sigma1
      if(sigma<=0.0_dp.or.nu<=0.0_dp.or.tau<=2.0_dp)then; v=nanv(); return; end if
      call sst_transform(mu,sigma,nu,tau,m1,m2,s1,mu1,sigma1)
      v=pST3(q,mu1,sigma1,nu,tau,lower_tail,log_p)
   end function pSST
   real(dp) function qSST(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::m1,m2,s1,mu1,sigma1
      if(sigma<=0.0_dp.or.nu<=0.0_dp.or.tau<=2.0_dp)then; v=nanv(); return; end if
      call sst_transform(mu,sigma,nu,tau,m1,m2,s1,mu1,sigma1)
      v=qST3(p,mu1,sigma1,nu,tau,lower_tail,log_p)
   end function qSST
   real(dp) function rSST(mu,sigma,nu,tau) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      real(dp)::u
      call random_number(u); v=qSST(u,mu,sigma,nu,tau)
   end function rSST
   elemental subroutine sst_transform(mu,sigma,nu,tau,m1,m2,s1,mu1,sigma1)
      real(dp),intent(in)::mu,sigma,nu,tau
      real(dp),intent(out)::m1,m2,s1,mu1,sigma1
      m1=2.0_dp*sqrt(tau)*(nu*nu-1.0_dp)/(tau-1.0_dp)
      m1=m1/(exp(log_gamma(0.5_dp)+log_gamma(tau/2.0_dp)-log_gamma((tau+1.0_dp)/2.0_dp))*nu)
      m2=tau*(nu**3+nu**(-3))/((tau-2.0_dp)*(nu+1.0_dp/nu))
      s1=sqrt(m2-m1*m1); mu1=mu-sigma*m1/s1; sigma1=sigma/s1
   end subroutine sst_transform

   elemental real(dp) function dGT(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::z,zt,ld
      if(sigma<=0.0_dp.or.nu<=0.0_dp.or.tau<=0.0_dp)then; v=nanv(); return; end if
      z=(x-mu)/sigma; zt=abs(z)**tau
      if(nu<1.0e6_dp)then
         ld=log(tau)-log(2.0_dp*sigma)-log(nu)/tau-log_gamma(1.0_dp/tau)-log_gamma(nu)
         ld=ld+log_gamma(nu+1.0_dp/tau)-(nu+1.0_dp/tau)*log(1.0_dp+zt/nu)
      else
         ld=log(tau)-log(2.0_dp*sigma)-log_gamma(1.0_dp/tau)-zt
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function dGT
   real(dp) function pGT(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,w,b,raw
      if(sigma<=0.0_dp.or.nu<=0.0_dp.or.tau<=0.0_dp)then; v=nanv(); return; end if
      z=(q-mu)/sigma
      if(abs(z)==0.0_dp)then
         raw=0.5_dp
      else if(nu<1.0e6_dp)then
         w=abs(z)**tau; b=regularized_beta(w/(nu+w),1.0_dp/tau,nu)
         raw=0.5_dp+0.5_dp*sign(b,z)
      else
         ! generalized-error limit; integrate only for this rare branch
         if(z<0.0_dp)then
            raw=0.5_dp-adaptive_integral(lim_pdf,z,0.0_dp,2.0e-10_dp,20)
         else
            raw=0.5_dp+adaptive_integral(lim_pdf,0.0_dp,z,2.0e-10_dp,20)
         end if
      end if
      v=finish_prob(raw,lower_tail,log_p)
   contains
      function lim_pdf(t) result(y)
         real(dp),intent(in)::t
         real(dp)::y
         y=dGT(mu+sigma*t,mu,sigma,nu,tau)*sigma
      end function lim_pdf
   end function pGT
   real(dp) function qGT(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,s,b,u,w
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=0.0_dp)then; v=-infv(); return; end if
      if(pp>=1.0_dp)then; v=infv(); return; end if
      if(sigma<=0.0_dp.or.nu<=0.0_dp.or.tau<=0.0_dp)then; v=nanv(); return; end if
      if(nu>=1.0e6_dp)then
         v=qgt_numeric(pp,mu,sigma,nu,tau); return
      end if
      s=merge(-1.0_dp,1.0_dp,pp<0.5_dp)
      b=2.0_dp*abs(pp-0.5_dp)
      u=qbeta_v(b,1.0_dp/tau,nu)
      w=nu*u/max(1.0e-300_dp,1.0_dp-u)
      v=mu+s*sigma*w**(1.0_dp/tau)
   end function qGT
   real(dp) function qgt_numeric(pp,mu,sigma,nu,tau) result(v)
      real(dp),intent(in)::pp,mu,sigma,nu,tau
      real(dp)::lo,hi
      lo=mu-8.0_dp*sigma; hi=mu+8.0_dp*sigma
      do while(pGT(lo,mu,sigma,nu,tau)>pp); lo=mu-2.0_dp*(mu-lo); end do
      do while(pGT(hi,mu,sigma,nu,tau)<pp); hi=mu+2.0_dp*(hi-mu); end do
      v=bisection_root(cdf,pp,lo,hi,2.0e-9_dp,100)
   contains
      function cdf(x) result(y)
         real(dp),intent(in)::x
         real(dp)::y
         y=pGT(x,mu,sigma,nu,tau)
      end function cdf
   end function qgt_numeric
   real(dp) function rGT(mu,sigma,nu,tau) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      real(dp)::u
      call random_number(u); v=qGT(u,mu,sigma,nu,tau)
   end function rGT

   elemental real(dp) function dexGAUS(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::z,pc,ld
      if(sigma<=0.0_dp.or.nu<=0.0_dp)then; v=nanv(); return; end if
      if(nu<=0.05_dp*sigma)then
         v=dnorm_v(x,mu,sigma,log_density); return
      end if
      z=x-mu-sigma*sigma/nu; pc=normal_cdf(z/sigma)
      if(pc<=0.0_dp)then; ld=-infv(); else; ld=-log(nu)-(z+sigma*sigma/(2.0_dp*nu))/nu+log(pc); end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function dexGAUS
   elemental real(dp) function pexGAUS(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::z,raw,e
      if(sigma<=0.0_dp.or.nu<=0.0_dp)then; v=nanv(); return; end if
      if(nu<=0.05_dp*sigma)then
         raw=normal_cdf((q-mu)/sigma)
      else
         z=q-mu-sigma*sigma/nu
         e=((mu+sigma*sigma/nu)**2-mu*mu-2.0_dp*q*sigma*sigma/nu)/(2.0_dp*sigma*sigma)
         raw=normal_cdf((q-mu)/sigma)-normal_cdf(z/sigma)*exp(e)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pexGAUS
   real(dp) function qexGAUS(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,lo,hi
      pp=input_prob(p,lower_tail,log_p)
      if(pp<=0.0_dp)then; v=-infv(); return; end if
      if(pp>=1.0_dp)then; v=infv(); return; end if
      lo=mu-8.0_dp*sigma; hi=mu+8.0_dp*sigma+12.0_dp*nu
      do while(pexGAUS(lo,mu,sigma,nu)>pp); lo=lo-4.0_dp*sigma; end do
      do while(pexGAUS(hi,mu,sigma,nu)<pp); hi=hi+max(4.0_dp*sigma,4.0_dp*nu); end do
      v=bisection_root(cdf,pp,lo,hi,2.0e-9_dp,120)
   contains
      function cdf(x) result(y)
         real(dp),intent(in)::x
         real(dp)::y
         y=pexGAUS(x,mu,sigma,nu)
      end function cdf
   end function qexGAUS
   real(dp) function rexGAUS(mu,sigma,nu) result(v)
      real(dp),intent(in)::mu,sigma,nu
      real(dp)::u
      call random_number(u); v=qexGAUS(u,mu,sigma,nu)
   end function rexGAUS

   elemental real(dp) function dPARETO(x,mu,log_density) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(mu<=0.0_dp)then; v=nanv(); return; end if
      if(x<=1.0_dp)then; v=merge(-infv(),0.0_dp,want_log(log_density)); return; end if
      ld=log(mu)-(mu+1.0_dp)*log(x); v=merge(ld,exp(ld),want_log(log_density))
   end function dPARETO
   elemental real(dp) function pPARETO(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(mu<=0.0_dp)then; v=nanv(); return; end if
      raw=merge(1.0_dp-q**(-mu),0.0_dp,q>1.0_dp); v=finish_prob(raw,lower_tail,log_p)
   end function pPARETO
   elemental real(dp) function qPARETO(p,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp
      pp=input_prob(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then; v=nanv(); else if(pp>=1.0_dp)then; v=infv(); &
      else; v=(1.0_dp-pp)**(-1.0_dp/mu); end if
   end function qPARETO
   real(dp) function rPARETO(mu) result(v)
      real(dp),intent(in)::mu; real(dp)::u
      call random_number(u); v=qPARETO(u,mu)
   end function rPARETO

   elemental real(dp) function dPARETO1(x,mu,log_density) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(mu<=0.0_dp)then; v=nanv(); return; end if
      if(x<=0.0_dp)then; v=merge(-infv(),0.0_dp,want_log(log_density)); return; end if
      ld=log(mu)-(mu+1.0_dp)*log(1.0_dp+x); v=merge(ld,exp(ld),want_log(log_density))
   end function dPARETO1
   elemental real(dp) function pPARETO1(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(mu<=0.0_dp)then; v=nanv(); return; end if
      raw=merge(1.0_dp-(1.0_dp+q)**(-mu),0.0_dp,q>0.0_dp); v=finish_prob(raw,lower_tail,log_p)
   end function pPARETO1
   elemental real(dp) function qPARETO1(p,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp
      pp=input_prob(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then; v=nanv(); else if(pp>=1.0_dp)then; v=infv(); &
      else; v=(1.0_dp-pp)**(-1.0_dp/mu)-1.0_dp; end if
   end function qPARETO1
   real(dp) function rPARETO1(mu) result(v)
      real(dp),intent(in)::mu; real(dp)::u
      call random_number(u); v=qPARETO1(u,mu)
   end function rPARETO1

   elemental real(dp) function dPARETO2(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      if(x<=0.0_dp)then; v=merge(-infv(),0.0_dp,want_log(log_density)); return; end if
      ld=-log(sigma)+log(mu)/sigma-(1.0_dp/sigma+1.0_dp)*log(x+mu)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dPARETO2
   elemental real(dp) function pPARETO2(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      raw=merge(1.0_dp-(mu/(mu+q))**(1.0_dp/sigma),0.0_dp,q>0.0_dp)
      v=finish_prob(raw,lower_tail,log_p)
   end function pPARETO2
   elemental real(dp) function qPARETO2(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp
      pp=input_prob(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then; v=nanv(); &
      else if(pp>=1.0_dp)then; v=infv(); else; v=mu*((1.0_dp-pp)**(-sigma)-1.0_dp); end if
   end function qPARETO2
   real(dp) function rPARETO2(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma; real(dp)::u
      call random_number(u); v=qPARETO2(u,mu,sigma)
   end function rPARETO2

   elemental real(dp) function dPARETO2o(x,mu,sigma,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_density
      real(dp)::ld
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      if(x<=0.0_dp)then; v=merge(-infv(),0.0_dp,want_log(log_density)); return; end if
      ld=log(sigma)+sigma*log(mu)-(sigma+1.0_dp)*log(x+mu)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dPARETO2o
   elemental real(dp) function pPARETO2o(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::raw
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      raw=merge(1.0_dp-(mu/(mu+q))**sigma,0.0_dp,q>0.0_dp); v=finish_prob(raw,lower_tail,log_p)
   end function pPARETO2o
   elemental real(dp) function qPARETO2o(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp
      pp=input_prob(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.pp<0.0_dp.or.pp>1.0_dp)then; v=nanv(); &
      else if(pp>=1.0_dp)then; v=infv(); else; v=mu*((1.0_dp-pp)**(-1.0_dp/sigma)-1.0_dp); end if
   end function qPARETO2o
   real(dp) function rPARETO2o(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma; real(dp)::u
      call random_number(u); v=qPARETO2o(u,mu,sigma)
   end function rPARETO2o
end module gamlss_continuous_v03
