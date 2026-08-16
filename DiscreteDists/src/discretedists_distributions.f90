module discretedists_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use discretedists_kinds, only : dp
   use discretedists_numerics, only : lambert_wm1
   use compoissonreg_distributions, only : pcmp, qcmp, rcmp, ecmp, vcmp
   use compoissonreg_normalizer, only : cmp_logz_hybrid
   implicit none
   private

   integer, parameter, public :: DD_QINF = huge(0)

   public :: dberg,pberg,qberg,rberg
   public :: dcompo,pcompo,qcompo,rcompo
   public :: dcompo2,pcompo2,qcompo2,rcompo2,mu_phi_2_lambda_nu_compo2
   public :: ddbh,pdbh,qdbh,rdbh
   public :: ddgeii,pdgeii,qdgeii,rdgeii
   public :: ddikum,pdikum,qdikum,rdikum
   public :: ddld,pdld,qdld,rdld
   public :: ddmolbe,pdmolbe,qdmolbe,rdmolbe
   public :: ddperks,pdperks,qdperks,rdperks
   public :: ddspa,pdspa,qdspa,rdspa
   public :: dggeo,pggeo,qggeo,rggeo
   public :: dhyperpo,phyperpo,qhyperpo,rhyperpo
   public :: dhyperpo2,phyperpo2,qhyperpo2,rhyperpo2
   public :: dpoisxl,ppoisxl,qpoisxl,rpoisxl
   public :: f11, ar, obtaining_lambda, mean_var_hp, mean_var_hp2, simulate_hp
   public :: compo_mean_exact, compo_variance_exact

contains

   pure real(dp) function nan_dp() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nan_dp

   pure real(dp) function pinf_dp() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function pinf_dp

   pure logical function want(flag,default) result(v)
      logical,intent(in),optional::flag
      logical,intent(in)::default
      v=default;if(present(flag))v=flag
   end function want

   pure real(dp) function prob_input(p,lower_tail,log_p) result(t)
      real(dp),intent(in)::p
      logical,intent(in),optional::lower_tail,log_p
      logical::lt,lp
      lt=want(lower_tail,.true.);lp=want(log_p,.false.)
      if(lp)then;t=exp(p);else;t=p;end if
      if(.not.lt)t=1.0_dp-t
   end function prob_input

   pure real(dp) function prob_output(p,lower_tail,log_p) result(t)
      real(dp),intent(in)::p
      logical,intent(in),optional::lower_tail,log_p
      logical::lt,lp
      lt=want(lower_tail,.true.);lp=want(log_p,.false.)
      t=p;if(.not.lt)t=1.0_dp-t
      t=min(1.0_dp,max(0.0_dp,t))
      if(lp)then
         if(t==0.0_dp)then;t=-pinf_dp();else;t=log(t);end if
      end if
   end function prob_output

   impure elemental real(dp) function dberg(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::p,r
      integer::k
      if(mu<0.0_dp.or.sigma<0.0_dp.or.sigma<=abs(mu-1.0_dp).or.x<0.0_dp)then
         p=0.0_dp
      else
         k=nint(x)
         if(k==0)then
            p=(1.0_dp-mu+sigma)/(1.0_dp+mu+sigma)
         else
            r=mu+sigma-1.0_dp
            p=4.0_dp*mu*r**real(k-1,dp)/(mu+sigma+1.0_dp)**real(k+1,dp)
         end if
      end if
      if(want(log_p,.false.))then
         v=log(max(p,tiny(1.0_dp)))
      else
         v=p
      end if
   end function dberg

   impure elemental real(dp) function pberg(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p,r
      integer::k
      if(mu<0.0_dp.or.sigma<0.0_dp.or.sigma<=abs(mu-1.0_dp))then
         p=0.0_dp
      else if(q<0.0_dp)then
         p=0.0_dp
      else
         k=int(floor(q));r=(mu+sigma-1.0_dp)/(mu+sigma+1.0_dp)
         p=(1.0_dp-mu+sigma)/(1.0_dp+mu+sigma)+ &
           2.0_dp*mu/(1.0_dp+mu+sigma)*(1.0_dp-r**real(k,dp))
      end if
      v=prob_output(p,lower_tail,log_p)
   end function pberg

   impure elemental integer function qberg(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,p0,r,z
      if(mu<0.0_dp.or.sigma<0.0_dp.or.sigma<=abs(mu-1.0_dp))then;q=-1;return;end if
      t=prob_input(p,lower_tail,log_p)
      if(t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      p0=(1.0_dp-mu+sigma)/(1.0_dp+mu+sigma)
      if(t<=p0)then;q=0;return;end if
      r=(mu+sigma-1.0_dp)/(mu+sigma+1.0_dp)
      z=(log1p_safe(-t)+log(mu+sigma+1.0_dp)-log(2.0_dp*mu))/log(r)
      q=max(1,ceiling(z))
   end function qberg

   subroutine rberg(n,mu,sigma,x)
      integer,intent(in)::n
      real(dp),intent(in)::mu,sigma
      integer,intent(out)::x(n)
      real(dp)::u
      integer::i
      do i=1,n;call random_number(u);x(i)=qberg(u,mu,sigma);end do
   end subroutine rberg

   pure real(dp) function log1p_safe(x) result(v)
      real(dp),intent(in)::x
      if(abs(x)<1.0e-8_dp)then
         v=x-x*x/2.0_dp+x**3/3.0_dp-x**4/4.0_dp
      else
         v=log(1.0_dp+x)
      end if
   end function log1p_safe

   impure elemental real(dp) function dcompo(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::lp
      if(mu<=0.0_dp.or.sigma<0.0_dp.or.x<0.0_dp)then
         lp=-pinf_dp()
      else
         lp=x*log(mu)-sigma*log_gamma(x+1.0_dp)-cmp_logz_hybrid(mu,sigma)
      end if
      if(want(log_p,.false.))then;v=lp;else;v=exp(lp);end if
   end function dcompo

   impure elemental real(dp) function pcompo(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(mu<=0.0_dp.or.sigma<0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp;else;p=pcmp(int(floor(q)),mu,sigma);end if
      v=prob_output(p,lower_tail,log_p)
   end function pcompo

   impure elemental integer function qcompo(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      q=qcmp(t,mu,sigma)
   end function qcompo

   subroutine rcompo(n,mu,sigma,x)
      integer,intent(in)::n
      real(dp),intent(in)::mu,sigma
      integer,intent(out)::x(n)
      call rcmp(n,mu,sigma,x)
   end subroutine rcompo

   pure subroutine mu_phi_2_lambda_nu_compo2(mu,phi,lambda,nu)
      real(dp),intent(in)::mu,phi
      real(dp),intent(out)::lambda,nu
      real(dp)::base
      nu=exp(phi);base=mu+(nu-1.0_dp)/(2.0_dp*nu)
      if(base<=0.0_dp)then;lambda=nan_dp();else;lambda=base**nu;end if
   end subroutine mu_phi_2_lambda_nu_compo2

   impure elemental real(dp) function dcompo2(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::lambda,nu
      call mu_phi_2_lambda_nu_compo2(mu,sigma,lambda,nu)
      v=dcompo(x,lambda,nu,log_p)
   end function dcompo2

   impure elemental real(dp) function pcompo2(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::lambda,nu
      call mu_phi_2_lambda_nu_compo2(mu,sigma,lambda,nu)
      v=pcompo(q,lambda,nu,lower_tail,log_p)
   end function pcompo2

   impure elemental integer function qcompo2(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::lambda,nu
      call mu_phi_2_lambda_nu_compo2(mu,sigma,lambda,nu)
      q=qcompo(p,lambda,nu,lower_tail,log_p)
   end function qcompo2

   subroutine rcompo2(n,mu,sigma,x)
      integer,intent(in)::n
      real(dp),intent(in)::mu,sigma
      integer,intent(out)::x(n)
      real(dp)::lambda,nu
      call mu_phi_2_lambda_nu_compo2(mu,sigma,lambda,nu)
      call rcmp(n,lambda,nu,x)
   end subroutine rcompo2

   impure elemental real(dp) function ddbh(x,mu,log_p) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_p
      real(dp)::lp,a
      if(mu<=0.0_dp.or.mu>1.0_dp.or.x<0.0_dp)then;lp=-pinf_dp()
      else
         a=1.0_dp/(x+1.0_dp)-mu/(x+2.0_dp)
         if(a<=0.0_dp)then;lp=-pinf_dp();else;lp=log(a)+x*log(mu);end if
      end if
      if(want(log_p,.false.))then;v=lp;else;v=exp(lp);end if
   end function ddbh

   impure elemental real(dp) function pdbh(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      integer::k
      if(mu<=0.0_dp.or.mu>1.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp;else
         k=int(floor(q));p=1.0_dp-mu**real(k+1,dp)/real(k+2,dp)
      end if
      v=prob_output(p,lower_tail,log_p)
   end function pdbh

   impure elemental integer function qdbh(p,mu,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t
      integer::lo,hi,mid
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.mu>1.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      lo=0;hi=1
      do while(pdbh(real(hi,dp),mu)<t.and.hi<huge(1)/2);hi=2*hi+1;end do
      do while(lo<hi)
         mid=lo+(hi-lo)/2
         if(pdbh(real(mid,dp),mu)>=t)then;hi=mid;else;lo=mid+1;end if
      end do
      q=lo
   end function qdbh

   subroutine rdbh(n,mu,x)
      integer,intent(in)::n;real(dp),intent(in)::mu;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdbh(u,mu);end do
   end subroutine rdbh

   impure elemental real(dp) function ddgeii(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::p
      integer::k
      if(mu<=0.0_dp.or.mu>=1.0_dp.or.sigma<=0.0_dp.or.x<0.0_dp.or.abs(x-nint(x))>sqrt(epsilon(1.0_dp)))then
         p=0.0_dp
      else
         k=nint(x);p=(1.0_dp-mu**real(k+1,dp))**sigma-(1.0_dp-mu**real(k,dp))**sigma
      end if
      if(want(log_p,.false.))then
         if(p<=0.0_dp)then;v=-pinf_dp();else;v=log(p);end if
      else;v=p;end if
   end function ddgeii

   impure elemental real(dp) function pdgeii(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      integer::k
      if(mu<=0.0_dp.or.mu>=1.0_dp.or.sigma<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp;else;k=int(floor(q));p=(1.0_dp-mu**real(k+1,dp))**sigma;end if
      v=prob_output(p,lower_tail,log_p)
   end function pdgeii

   impure elemental integer function qdgeii(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,z
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.mu>=1.0_dp.or.sigma<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      z=log(max(tiny(1.0_dp),1.0_dp-t**(1.0_dp/sigma)))/log(mu)-1.0_dp
      q=max(0,ceiling(z))
   end function qdgeii

   subroutine rdgeii(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdgeii(u,mu,sigma);end do
   end subroutine rdgeii

   impure elemental real(dp) function ddikum(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::p
      if(x<0.0_dp.or.mu<=0.0_dp.or.sigma<=0.0_dp)then;p=0.0_dp
      else;p=(1.0_dp-(x+2.0_dp)**(-mu))**sigma-(1.0_dp-(x+1.0_dp)**(-mu))**sigma;end if
      if(want(log_p,.false.))then
         if(p<=0.0_dp)then;v=-pinf_dp();else;v=log(p);end if
      else;v=p;end if
   end function ddikum

   impure elemental real(dp) function pdikum(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp;else;k=int(floor(q));p=(1.0_dp-(real(k,dp)+2.0_dp)**(-mu))**sigma;end if
      v=prob_output(p,lower_tail,log_p)
   end function pdikum

   impure elemental integer function qdikum(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,z
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      z=(1.0_dp-t**(1.0_dp/sigma))**(-1.0_dp/mu)-2.0_dp
      q=max(0,ceiling(z))
   end function qdikum

   subroutine rdikum(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdikum(u,mu,sigma);end do
   end subroutine rdikum

   impure elemental real(dp) function ddld(x,mu,log_p) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_p
      real(dp)::lp,a,e
      if(mu<=0.0_dp.or.x<0.0_dp)then;lp=-pinf_dp()
      else
         e=exp(-mu);a=mu*(1.0_dp-2.0_dp*e)+(1.0_dp-e)*(1.0_dp+mu*x)
         if(a<=0.0_dp)then;lp=-pinf_dp();else;lp=-mu*x-log(1.0_dp+mu)+log(a);end if
      end if
      if(want(log_p,.false.))then;v=lp;else;v=exp(lp);end if
   end function ddld

   impure elemental real(dp) function pdld(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p,t
      integer::k
      if(mu<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else
         ! Correct discrete CDF: upstream pDLD is shifted by one count.
         k=int(floor(q));t=real(k+1,dp)
         p=1.0_dp-((1.0_dp+mu+mu*t)/(1.0_dp+mu))*exp(-mu*t)
      end if
      v=prob_output(p,lower_tail,log_p)
   end function pdld

   impure elemental integer function qdld(p,mu,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,z,w,val
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      z=-(1.0_dp+mu)*exp(-1.0_dp-mu)*(1.0_dp-t);w=lambert_wm1(z)
      val=-1.0_dp-1.0_dp/mu-w/mu
      q=max(0,int(floor(val+16.0_dp*epsilon(1.0_dp))))
      do while(q>0.and.pdld(real(q-1,dp),mu)>=t);q=q-1;end do
      do while(pdld(real(q,dp),mu)<t);q=q+1;end do
   end function qdld

   subroutine rdld(n,mu,x)
      integer,intent(in)::n;real(dp),intent(in)::mu;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdld(u,mu);end do
   end subroutine rdld

   pure real(dp) function dmolbe_surv_edge(k,mu,sigma) result(s)
      integer,intent(in)::k
      real(dp),intent(in)::mu,sigma
      real(dp)::t
      t=(1.0_dp+real(k,dp)/mu)*exp(-real(k,dp)/mu)
      s=sigma*t/(1.0_dp-(1.0_dp-sigma)*t)
   end function dmolbe_surv_edge

   impure elemental real(dp) function ddmolbe(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::p,t1,t2,b1,b2
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.x<0.0_dp.or.abs(x-nint(x))>sqrt(epsilon(1.0_dp)))then;p=0.0_dp
      else
         k=nint(x);t1=(1.0_dp+real(k,dp)/mu)*exp(-real(k,dp)/mu)
         t2=(1.0_dp+real(k+1,dp)/mu)*exp(-real(k+1,dp)/mu)
         b1=1.0_dp-(1.0_dp-sigma)*t1;b2=1.0_dp-(1.0_dp-sigma)*t2
         p=sigma*(t1-t2)/(b1*b2)
      end if
      if(want(log_p,.false.))then
         if(p<=0.0_dp)then;v=-pinf_dp();else;v=log(p);end if
      else;v=p;end if
   end function ddmolbe

   impure elemental real(dp) function pdmolbe(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p,t
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else
         k=int(floor(q));t=(1.0_dp+real(k+1,dp)/mu)*exp(-real(k+1,dp)/mu)
         p=(1.0_dp-t)/(1.0_dp-(1.0_dp-sigma)*t)
      end if
      v=prob_output(p,lower_tail,log_p)
   end function pdmolbe

   impure elemental integer function qdmolbe(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t
      integer::lo,hi,mid
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      lo=0;hi=1
      do while(pdmolbe(real(hi,dp),mu,sigma)<t.and.hi<huge(1)/2);hi=2*hi+1;end do
      do while(lo<hi);mid=lo+(hi-lo)/2;if(pdmolbe(real(mid,dp),mu,sigma)>=t)then;hi=mid;else;lo=mid+1;end if;end do
      q=lo
   end function qdmolbe

   subroutine rdmolbe(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdmolbe(u,mu,sigma);end do
   end subroutine rdmolbe

   impure elemental real(dp) function ddperks(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::lp
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.x<0.0_dp.or.abs(x-nint(x))>sqrt(epsilon(1.0_dp)))then;lp=-pinf_dp()
      else
         lp=log(mu)+log(1.0_dp+mu)+log(exp(sigma)-1.0_dp)+sigma*x- &
            log(1.0_dp+mu*exp(sigma*x))-log(1.0_dp+mu*exp(sigma*(x+1.0_dp)))
      end if
      if(want(log_p,.false.))then;v=lp;else;v=exp(lp);end if
   end function ddperks

   impure elemental real(dp) function pdperks(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p,z
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else;k=int(floor(q));z=exp(sigma*real(k+1,dp));p=mu*(z-1.0_dp)/(1.0_dp+mu*z);end if
      v=prob_output(p,lower_tail,log_p)
   end function pdperks

   impure elemental integer function qdperks(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,z
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      z=log((t+mu)/(mu*(1.0_dp-t)))/sigma-1.0_dp;q=max(0,ceiling(z))
   end function qdperks

   subroutine rdperks(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdperks(u,mu,sigma);end do
   end subroutine rdperks

   impure elemental real(dp) function ddspa(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::p,l,t0,t1
      if(x<0.0_dp.or.mu<=0.0_dp.or.sigma<=0.0_dp.or.sigma>=1.0_dp)then;p=0.0_dp
      else
         l=log(sigma);t0=x**mu;t1=(x+1.0_dp)**mu
         p=exp(l*t0)*(1.0_dp-t0*l)-exp(l*t1)*(1.0_dp-t1*l)
      end if
      if(want(log_p,.false.))then
         if(p<=0.0_dp)then;v=-pinf_dp();else;v=log(p);end if
      else;v=p;end if
   end function ddspa

   impure elemental real(dp) function pdspa(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p,t,l
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.sigma>=1.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else;k=int(floor(q));t=real(k+1,dp)**mu;l=log(sigma);p=1.0_dp-exp(l*t)*(1.0_dp-t*l);end if
      v=prob_output(p,lower_tail,log_p)
   end function pdspa

   impure elemental integer function qdspa(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t
      integer::lo,hi,mid
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.sigma>=1.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      lo=0;hi=1
      do while(pdspa(real(hi,dp),mu,sigma)<t.and.hi<huge(1)/2);hi=2*hi+1;end do
      do while(lo<hi);mid=lo+(hi-lo)/2;if(pdspa(real(mid,dp),mu,sigma)>=t)then;hi=mid;else;lo=mid+1;end if;end do
      q=lo
   end function qdspa

   subroutine rdspa(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qdspa(u,mu,sigma);end do
   end subroutine rdspa

   impure elemental real(dp) function dggeo(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::lp
      if(mu<=0.0_dp.or.mu>=1.0_dp.or.sigma<=0.0_dp.or.x<0.0_dp)then;lp=-pinf_dp()
      else
         lp=log(sigma)+x*log(mu)+log(1.0_dp-mu)-log(1.0_dp-(1.0_dp-sigma)*mu**(x+1.0_dp))- &
            log(1.0_dp-(1.0_dp-sigma)*mu**x)
      end if
      if(want(log_p,.false.))then;v=lp;else;v=exp(lp);end if
   end function dggeo

   impure elemental real(dp) function pggeo(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p,t
      integer::k
      if(mu<=0.0_dp.or.mu>=1.0_dp.or.sigma<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else;k=int(floor(q));t=mu**real(k+1,dp);p=(1.0_dp-t)/(1.0_dp-(1.0_dp-sigma)*t);end if
      v=prob_output(p,lower_tail,log_p)
   end function pggeo

   impure elemental integer function qggeo(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,z,r
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.mu>=1.0_dp.or.sigma<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      r=(1.0_dp-t)/(1.0_dp-t*(1.0_dp-sigma));z=log(r)/log(mu)-1.0_dp;q=max(0,ceiling(z))
   end function qggeo

   subroutine rggeo(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qggeo(u,mu,sigma);end do
   end subroutine rggeo

   pure real(dp) function f11(z,c,maxiter_series,tol) result(value)
      real(dp),intent(in)::z,c
      integer,intent(in),optional::maxiter_series
      real(dp),intent(in),optional::tol
      integer::n,m
      real(dp)::fac,temp,series,l,eps
      m=10000;if(present(maxiter_series))m=maxiter_series
      eps=1.0e-10_dp;if(present(tol))eps=tol
      fac=1.0_dp;temp=1.0_dp;series=1.0_dp;l=c
      if(c<=0.0_dp)then;value=nan_dp();return;end if
      do n=1,m
         fac=fac*z/l;series=temp+fac
         if(abs(series-temp)<=eps*max(1.0_dp,abs(series)))exit
         temp=series;l=l+1.0_dp
      end do
      value=series
   end function f11

   pure real(dp) function ar(a,r) result(v)
      real(dp),intent(in)::a,r
      if(a<=0.0_dp.or.a+r<=0.0_dp)then;v=nan_dp();else;v=exp(log_gamma(a+r)-log_gamma(a));end if
   end function ar

   pure real(dp) function dhyperpo_core(x,lambda,sigma) result(p)
      real(dp),intent(in)::x,lambda,sigma
      if(x<0.0_dp.or.lambda<=0.0_dp.or.sigma<=0.0_dp)then;p=0.0_dp;return;end if
      p=exp(x*log(lambda)-log_gamma(sigma+x)+log_gamma(sigma)-log(f11(lambda,sigma)))
   end function dhyperpo_core

   impure elemental real(dp) function dhyperpo(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::p
      p=dhyperpo_core(x,mu,sigma)
      if(want(log_p,.false.))then;if(p<=0.0_dp)then;v=-pinf_dp();else;v=log(p);end if;else;v=p;end if
   end function dhyperpo

   impure elemental real(dp) function phyperpo(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      integer::k,j
      real(dp)::p,term,z
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else
         k=int(floor(q));z=f11(mu,sigma);term=1.0_dp/z;p=term
         do j=1,k;term=term*mu/(sigma+real(j-1,dp));p=p+term;end do
         p=min(1.0_dp,p)
      end if
      v=prob_output(p,lower_tail,log_p)
   end function phyperpo

   impure elemental integer function qhyperpo(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t,term,cum,z
      integer::j
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.sigma<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      z=f11(mu,sigma);term=1.0_dp/z;cum=term
      if(t<=cum)then;q=0;return;end if
      do j=1,1000000
         term=term*mu/(sigma+real(j-1,dp));cum=cum+term
         if(t<=cum.or.term<=tiny(1.0_dp))then;q=j;return;end if
      end do
      q=DD_QINF
   end function qhyperpo

   subroutine rhyperpo(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qhyperpo(u,mu,sigma);end do
   end subroutine rhyperpo

   impure real(dp) function obtaining_lambda(media,gamma,tol,max_iter) result(lambda)
      real(dp),intent(in)::media,gamma
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::max_iter
      real(dp)::lo,hi,fl,fm,mid,eps
      integer::it,m
      if(media<=0.0_dp.or.gamma<=0.0_dp)then;lambda=nan_dp();return;end if
      if(abs(gamma-1.0_dp)<=epsilon(1.0_dp))then;lambda=media;return;end if
      eps=1.0e-10_dp;if(present(tol))eps=tol;m=1000;if(present(max_iter))m=max_iter
      lo=min(media,max(media+gamma-1.0_dp,gamma*media));hi=max(media,min(media+gamma-1.0_dp,gamma*media))
      fl=hp_mean_equation(lo,media,gamma)
      if(fl*hp_mean_equation(hi,media,gamma)>0.0_dp)then;lambda=media;return;end if
      mid=0.5_dp*(lo+hi)
      do it=1,m
         mid=0.5_dp*(lo+hi);fm=hp_mean_equation(mid,media,gamma)
         if(abs(fm)<eps.or.0.5_dp*(hi-lo)<eps)exit
         if(fl*fm<0.0_dp)then;hi=mid;else;lo=mid;fl=fm;end if
      end do
      lambda=mid
   end function obtaining_lambda

   pure real(dp) function hp_mean_equation(x,media,gamma) result(v)
      real(dp),intent(in)::x,media,gamma
      v=x-(gamma-1.0_dp)*(1.0_dp-1.0_dp/f11(x,gamma))-media
   end function hp_mean_equation

   impure elemental real(dp) function dhyperpo2(x,mu,sigma,log_p) result(v)
      real(dp),intent(in)::x,mu,sigma
      logical,intent(in),optional::log_p
      real(dp)::lambda
      lambda=obtaining_lambda(mu,sigma);v=dhyperpo(x,lambda,sigma,log_p)
   end function dhyperpo2

   impure elemental real(dp) function phyperpo2(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::lambda
      lambda=obtaining_lambda(mu,sigma);v=phyperpo(q,lambda,sigma,lower_tail,log_p)
   end function phyperpo2

   impure elemental integer function qhyperpo2(p,mu,sigma,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu,sigma
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::lambda
      lambda=obtaining_lambda(mu,sigma);q=qhyperpo(p,lambda,sigma,lower_tail,log_p)
   end function qhyperpo2

   subroutine rhyperpo2(n,mu,sigma,x)
      integer,intent(in)::n;real(dp),intent(in)::mu,sigma;integer,intent(out)::x(n)
      real(dp)::lambda
      lambda=obtaining_lambda(mu,sigma);call rhyperpo(n,lambda,sigma,x)
   end subroutine rhyperpo2

   subroutine mean_var_hp(mu,sigma,mean,variance)
      real(dp),intent(in)::mu,sigma
      real(dp),intent(out)::mean,variance
      real(dp)::z
      z=f11(mu,sigma);mean=mu-(sigma-1.0_dp)*(z-1.0_dp)/z
      variance=mu+(mu-(sigma-1.0_dp))*mean-mean*mean
   end subroutine mean_var_hp

   subroutine mean_var_hp2(mu,sigma,mean,variance)
      real(dp),intent(in)::mu,sigma
      real(dp),intent(out)::mean,variance
      real(dp)::lambda
      lambda=obtaining_lambda(mu,sigma);mean=mu
      variance=lambda+(lambda-(sigma-1.0_dp))*mean-mean*mean
   end subroutine mean_var_hp2

   integer function simulate_hp(sigma,mu) result(y)
      real(dp),intent(in)::sigma,mu
      real(dp)::u
      call random_number(u);y=qhyperpo(u,mu,sigma)
   end function simulate_hp

   impure elemental real(dp) function dpoisxl(x,mu,log_p) result(v)
      real(dp),intent(in)::x,mu
      logical,intent(in),optional::log_p
      real(dp)::lp
      if(x<0.0_dp.or.mu<=0.0_dp)then;lp=-pinf_dp()
      else;lp=2.0_dp*log(mu)+log(x+mu*mu+3.0_dp*(1.0_dp+mu))-(4.0_dp+x)*log(1.0_dp+mu);end if
      if(want(log_p,.false.))then;v=lp;else;v=exp(lp);end if
   end function dpoisxl

   impure elemental real(dp) function ppoisxl(q,mu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      integer::k
      if(mu<=0.0_dp)then;v=nan_dp();return;end if
      if(q<0.0_dp)then;p=0.0_dp
      else;k=int(floor(q));p=1.0_dp-(1.0_dp+mu*(real(k,dp)+4.0_dp+mu*(3.0_dp+mu)))/(1.0_dp+mu)**real(k+4,dp);end if
      v=prob_output(p,lower_tail,log_p)
   end function ppoisxl

   impure elemental integer function qpoisxl(p,mu,lower_tail,log_p) result(q)
      real(dp),intent(in)::p,mu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::t
      integer::lo,hi,mid
      t=prob_input(p,lower_tail,log_p)
      if(mu<=0.0_dp.or.t<0.0_dp.or.t>1.0_dp)then;q=-1;return;end if
      if(t<=0.0_dp)then;q=0;return;end if
      if(t>=1.0_dp)then;q=DD_QINF;return;end if
      lo=0;hi=1;do while(ppoisxl(real(hi,dp),mu)<t.and.hi<huge(1)/2);hi=2*hi+1;end do
      do while(lo<hi);mid=lo+(hi-lo)/2;if(ppoisxl(real(mid,dp),mu)>=t)then;hi=mid;else;lo=mid+1;end if;end do
      q=lo
   end function qpoisxl

   subroutine rpoisxl(n,mu,x)
      integer,intent(in)::n;real(dp),intent(in)::mu;integer,intent(out)::x(n)
      real(dp)::u;integer::i
      do i=1,n;call random_number(u);x(i)=qpoisxl(u,mu);end do
   end subroutine rpoisxl

   impure real(dp) function compo_mean_exact(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma
      v=ecmp(mu,sigma)
   end function compo_mean_exact

   impure real(dp) function compo_variance_exact(mu,sigma) result(v)
      real(dp),intent(in)::mu,sigma
      v=vcmp(mu,sigma)
   end function compo_variance_exact

end module discretedists_distributions
