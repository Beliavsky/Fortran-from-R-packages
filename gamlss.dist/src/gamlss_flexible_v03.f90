! Flexible gamma/negative-binomial families remaining from gamlss.dist.
! Upstream gamlss.dist GPL-2 | GPL-3. Translation SPDX-License-Identifier: GPL-3.0-only
module gamlss_flexible_v03
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp
   use gamlss_base, only : dgamma_v, pgamma_v, qgamma_v, rgamma_v, dpois_v, ppois_v, qpois_v, rpois_v, &
                           dnbinom_v, pnbinom_v, qnbinom_v, rnbinom_v
   implicit none
   private
   public :: dGAF,pGAF,qGAF,rGAF,dNBF,pNBF,qNBF,rNBF,dZINBF,pZINBF,qZINBF,rZINBF
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
            if(v==0.0_dp)then; v=-infv(); else; v=log(v); end if
         end if
      end if
   end function finish_prob
   elemental logical function integer_value(x) result(v)
      real(dp),intent(in)::x
      v=abs(x-real(nint(x),dp))<=16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x))
   end function integer_value

   elemental real(dp) function dGAF(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::s1
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      if(x<=0.0_dp)then
         v=merge(-infv(),0.0_dp,want_log(log_density)); return
      end if
      s1=sigma*mu**(0.5_dp*nu-1.0_dp)
      v=dgamma_v(x,1.0_dp/(s1*s1),mu*s1*s1,log_density)
   end function dGAF
   elemental real(dp) function pGAF(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::s1,p
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      s1=sigma*mu**(0.5_dp*nu-1.0_dp)
      p=pgamma_v(q,1.0_dp/(s1*s1),mu*s1*s1)
      v=finish_prob(p,lower_tail,log_p)
   end function pGAF
   real(dp) function qGAF(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::s1,pp
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      pp=input_prob(p,lower_tail,log_p)
      if(pp<0.0_dp.or.pp>1.0_dp)then; v=nanv(); return; end if
      s1=sigma*mu**(0.5_dp*nu-1.0_dp)
      v=qgamma_v(pp,1.0_dp/(s1*s1),mu*s1*s1)
   end function qGAF
   real(dp) function rGAF(mu,sigma,nu) result(v)
      real(dp),intent(in)::mu,sigma,nu
      real(dp)::s1
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      s1=sigma*mu**(0.5_dp*nu-1.0_dp)
      v=rgamma_v(1.0_dp/(s1*s1),mu*s1*s1)
   end function rGAF

   elemental real(dp) function dNBF(x,mu,sigma,nu,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu
      logical,intent(in),optional::log_density
      real(dp)::s1,pr
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      if(x<0.0_dp.or..not.integer_value(x))then
         v=merge(-infv(),0.0_dp,want_log(log_density)); return
      end if
      k=nint(x); s1=sigma*mu**(nu-2.0_dp)
      if(s1<1.0e-4_dp)then
         v=dpois_v(k,mu,log_density)
      else
         pr=1.0_dp/(1.0_dp+mu*s1)
         v=dnbinom_v(k,1.0_dp/s1,pr,log_density)
      end if
   end function dNBF
   elemental real(dp) function pNBF(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::s1,pr,p
      integer::k
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=nanv(); return; end if
      if(q<0.0_dp)then; v=finish_prob(0.0_dp,lower_tail,log_p); return; end if
      k=floor(q); s1=sigma*mu**(nu-2.0_dp)
      if(s1<1.0e-4_dp)then
         p=ppois_v(k,mu)
      else
         pr=1.0_dp/(1.0_dp+mu*s1)
         p=pnbinom_v(k,1.0_dp/s1,pr)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function pNBF
   integer function qNBF(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::s1,pr,pp
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=-1; return; end if
      pp=input_prob(p,lower_tail,log_p)
      if(pp<0.0_dp.or.pp>1.0_dp)then; v=-1; return; end if
      s1=sigma*mu**(nu-2.0_dp)
      if(s1<1.0e-4_dp)then
         v=qpois_v(pp,mu)
      else
         pr=1.0_dp/(1.0_dp+mu*s1)
         v=qnbinom_v(pp,1.0_dp/s1,pr)
      end if
   end function qNBF
   integer function rNBF(mu,sigma,nu) result(v)
      real(dp),intent(in)::mu,sigma,nu
      real(dp)::s1,pr
      if(mu<=0.0_dp.or.sigma<=0.0_dp)then; v=-1; return; end if
      s1=sigma*mu**(nu-2.0_dp)
      if(s1<1.0e-4_dp)then
         v=rpois_v(mu)
      else
         pr=1.0_dp/(1.0_dp+mu*s1)
         v=rnbinom_v(1.0_dp/s1,pr)
      end if
   end function rNBF

   elemental real(dp) function dZINBF(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp),intent(in)::x,mu,sigma,nu,tau
      logical,intent(in),optional::log_density
      real(dp)::f,ld
      if(tau<=0.0_dp.or.tau>=1.0_dp)then; v=nanv(); return; end if
      if(x<0.0_dp.or..not.integer_value(x))then
         v=merge(-infv(),0.0_dp,want_log(log_density)); return
      end if
      f=dNBF(x,mu,sigma,nu)
      if(nint(x)==0)then
         ld=log(tau+(1.0_dp-tau)*f)
      else
         if(f<=0.0_dp)then; ld=-infv(); else; ld=log(1.0_dp-tau)+log(f); end if
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function dZINBF
   elemental real(dp) function pZINBF(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::q,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::p
      if(tau<=0.0_dp.or.tau>=1.0_dp)then; v=nanv(); return; end if
      if(q<0.0_dp)then
         p=0.0_dp
      else
         p=tau+(1.0_dp-tau)*pNBF(q,mu,sigma,nu)
      end if
      v=finish_prob(p,lower_tail,log_p)
   end function pZINBF
   integer function qZINBF(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,mu,sigma,nu,tau
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp,pnew
      if(tau<=0.0_dp.or.tau>=1.0_dp)then; v=-1; return; end if
      pp=input_prob(p,lower_tail,log_p)
      if(pp<0.0_dp.or.pp>1.0_dp)then; v=-1; return; end if
      pnew=max(0.0_dp,(pp-tau)/(1.0_dp-tau)-1.0e-7_dp)
      v=qNBF(pnew,mu,sigma,nu)
   end function qZINBF
   integer function rZINBF(mu,sigma,nu,tau) result(v)
      real(dp),intent(in)::mu,sigma,nu,tau
      real(dp)::u
      if(tau<=0.0_dp.or.tau>=1.0_dp)then; v=-1; return; end if
      call random_number(u)
      v=qZINBF(u,mu,sigma,nu,tau)
   end function rZINBF
end module gamlss_flexible_v03
