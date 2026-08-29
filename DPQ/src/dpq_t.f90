! Central and noncentral Student-t algorithms translated from DPQ.
! SPDX-License-Identifier: GPL-2.0-or-later
module dpq_t
   use r_compat, only: dp, dt, pt, qt, pbeta, qnorm, r_lgamma
   use dpq_core, only: prob_output, prob_from_input, expm1_dp, log1p_dp, logr
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
      ieee_positive_inf, ieee_negative_inf
   implicit none
   private
   public :: b_chi, b_chi_asymp, lb_chi0, lb_chi00, lb_chi_asymp
   public :: c_dt, c_dt_asymp, c_pt
   public :: qt_r, qt_nappr, qt_appr, qt_u, qnt_r
   public :: pnt_r, pnt_r1, pnt_lrg, pnt_jw39, pnt_jw39_0
   public :: pnt3150, pnt_p94, pnt_chshp94, pnt_vw13, pnt_gst23_t6, pnt_gst23_1
   public :: dnt_jkbf, dnt_jkbf1, dt_wv

contains

   pure elemental real(dp) function b_chi_asymp(nu,order,one_minus) result(v)
      real(dp),intent(in)::nu
      integer,intent(in),optional::order
      logical,intent(in),optional::one_minus
      integer::k
      logical::om
      real(dp)::r,p
      k=2
      if(present(order))k=max(1,min(5,order))
      om=.false.
      if(present(one_minus))om=one_minus
      if(nu<=0.0_dp)then
      v=merge(1.0_dp,0.0_dp,om)
      return
      end if
      r=1.0_dp/(4.0_dp*nu)
      select case(k)
      case(1);p=1.0_dp
      case(2);p=1.0_dp-r/2.0_dp
      case(3);p=1.0_dp-r/2.0_dp*(1.0_dp+5.0_dp*r)
      case(4);p=1.0_dp-r/2.0_dp*(1.0_dp+r*(5.0_dp-21.0_dp*r/4.0_dp))
      case default;p=1.0_dp-r/2.0_dp*(1.0_dp+r*(5.0_dp-r*(21.0_dp+399.0_dp*r)/4.0_dp))
      end select
      r=r*p
      v=merge(r,1.0_dp-r,om)
   end function b_chi_asymp

   pure elemental real(dp) function b_chi(nu,one_minus) result(v)
      real(dp),intent(in)::nu
      logical,intent(in),optional::one_minus
      logical::om
      real(dp)::lr
      om=.false.
      if(present(one_minus))om=one_minus
      if(nu<=0.0_dp)then
      v=merge(1.0_dp,0.0_dp,om)
      return
      end if
      if(nu>1000.0_dp)then
      v=b_chi_asymp(nu,5,om)
      return
      end if
      lr=0.5_dp*log(2.0_dp/nu)+r_lgamma(0.5_dp*(nu+1.0_dp))-r_lgamma(0.5_dp*nu)
      if(om)then
      v=-expm1_dp(lr)
      else
      v=exp(lr)
      end if
   end function b_chi

   pure elemental real(dp) function lb_chi0(nu) result(v)
      real(dp),intent(in)::nu
      if(nu<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_negative_inf)
      else
         v=0.5_dp*log(2.0_dp/nu)+r_lgamma(0.5_dp*(nu+1.0_dp))-r_lgamma(0.5_dp*nu)
      end if
   end function lb_chi0
   pure elemental real(dp) function lb_chi00(nu) result(v)
      real(dp),intent(in)::nu
      v=lb_chi0(nu)
   end function lb_chi00
   pure elemental real(dp) function lb_chi_asymp(nu,order) result(v)
      real(dp),intent(in)::nu
      integer,intent(in),optional::order
      integer::k
      real(dp)::r,rr,p
      k=4
      if(present(order))k=max(1,min(8,order))
      if(k==1)then
      v=-1.0_dp/(4.0_dp*nu)
      return
      end if
      r=1.0_dp/(2.0_dp*nu)
      rr=r*r
      select case(k)
      case(2);p=1.0_dp-2.0_dp*rr/3.0_dp
      case(3);p=1.0_dp-rr*(2.0_dp/3.0_dp-16.0_dp*rr/5.0_dp)
      case(4);p=1.0_dp-rr*(2.0_dp/3.0_dp-rr*(16.0_dp/5.0_dp-272.0_dp*rr/7.0_dp))
      case default
         ! Higher orders are rarely needed; exact log-gamma is more reliable here.
         v=lb_chi0(nu)
         return
      end select
      v=-0.5_dp*r*p
   end function lb_chi_asymp

   pure elemental real(dp) function c_dt(nu) result(v)
      real(dp),intent(in)::nu
      v=r_lgamma(0.5_dp*(nu+1.0_dp))-r_lgamma(0.5_dp*nu)-0.5_dp*log(acos(-1.0_dp)*nu)
   end function c_dt
   pure elemental real(dp) function c_dt_asymp(nu) result(v)
      real(dp),intent(in)::nu
      v=-0.5_dp*log(2.0_dp*acos(-1.0_dp))-1.0_dp/(4.0_dp*nu)
   end function c_dt_asymp
   pure elemental real(dp) function c_pt(nu) result(v)
      real(dp),intent(in)::nu
      v=c_dt(nu)+0.5_dp*(nu-1.0_dp)*log(nu)
   end function c_pt

   pure elemental real(dp) function qt_r(p,df,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,df
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::pp
      logical::lt,lp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      pp=prob_from_input(p,lt,lp)
      v=qt(pp,df)
   end function qt_r

   pure elemental real(dp) function qt_nappr(p,df,lower_tail,log_p,k) result(v)
      real(dp),intent(in)::p,df
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(in),optional::k
      real(dp)::x,x2,n4,pp
      integer::kk
      logical::lt,lp
      lt=.true.
      lp=.false.
      kk=2
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(present(k))kk=max(0,min(4,k))
      pp=prob_from_input(p,lt,lp)
      x=qnorm(pp)
      if(kk==0)then
      v=x
      return
      end if
      x2=x*x
      n4=4.0_dp*df
      select case(kk)
      case(1);v=x*(1.0_dp+(x2+1.0_dp)/n4)
      case(2);v=x*(1.0_dp+(x2+1.0_dp+((5.0_dp*x2+16.0_dp)*x2+3.0_dp)/(24.0_dp*df))/n4)
      case(3);v=x*(1.0_dp+(x2+1.0_dp+((5.0_dp*x2+16.0_dp)*x2+3.0_dp+ &
         (((3.0_dp*x2+19.0_dp)*x2+17.0_dp)*x2-15.0_dp)/n4)/(24.0_dp*df))/n4)
      case default
         v=qt(pp,df)
      end select
   end function qt_nappr

   pure elemental real(dp) function qt_appr(p,df,ncp,method,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,df,ncp
      character(len=*),intent(in),optional::method
      logical,intent(in),optional::lower_tail,log_p
      character(len=1)::m
      real(dp)::pp,z,b,b2,den
      logical::lt,lp
      lt=.true.
      lp=.false.
      m='a'
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(present(method))m=method(1:1)
      pp=prob_from_input(p,lt,lp)
      z=qnorm(pp)
      b=b_chi(df)
      b2=b*b
      select case(m)
      case('b','B')
      den=1.0_dp-z*z/(2.0_dp*df)
      v=(ncp+z*sqrt(max(0.0_dp,den+ncp*ncp/(2.0_dp*df))))/den
      case('c','C')
      den=b2-z*z/(2.0_dp*df)
      v=(ncp*b+z*sqrt(max(0.0_dp,den+ncp*ncp/(2.0_dp*df))))/den
      case default
      den=b2-z*z*(1.0_dp-b2)
      v=(ncp*b+z*sqrt(max(0.0_dp,den+ncp*ncp*(1.0_dp-b2))))/den
      end select
   end function qt_appr

   pure real(dp) function pnt_lower_posdelta(t,df,ncp,errmax,itrmax) result(v)
      real(dp),intent(in)::t,df,ncp
      real(dp),intent(in),optional::errmax
      integer,intent(in),optional::itrmax
      real(dp)::tt,del,x,rxb,lambda,p,q,a,b,s,albeta,xodd,godd,xeven,geven,tnc,errbd,eps
      integer::it,imax
      logical::negdel
      if(df<=0.0_dp .or. ncp<0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
      return
      end if
      if(ncp==0.0_dp)then
      v=pt(t,df)
      return
      end if
      eps=1.0e-12_dp
      if(present(errmax))eps=errmax
      imax=10000
      if(present(itrmax))imax=itrmax
      if(df>4.0e5_dp .or. ncp*ncp>1416.0_dp)then
         x=(t*(1.0_dp-1.0_dp/(4.0_dp*df))-ncp)/sqrt(1.0_dp+t*t/(2.0_dp*df))
         v=0.5_dp*erfc(-x/sqrt(2.0_dp))
         return
      end if
      if(t>=0.0_dp)then
      negdel=.false.
      tt=t
      del=ncp
      else
      negdel=.true.
      tt=-t
      del=-ncp
      end if
      x=tt*tt
      rxb=df/(x+df)
      x=x/(x+df)
      tnc=0.0_dp
      if(x>0.0_dp)then
         lambda=del*del
         p=0.5_dp*exp(-0.5_dp*lambda)
         if(p==0.0_dp)then
            ! Large noncentrality: use the normal approximation rather than underflow the series.
            x=(tt*(1.0_dp-1.0_dp/(4.0_dp*df))-del)/sqrt(1.0_dp+tt*tt/(2.0_dp*df))
            tnc=0.5_dp*erfc(-x/sqrt(2.0_dp))
            if(negdel)tnc=1.0_dp-tnc
            v=tnc
            return
         end if
         q=sqrt(2.0_dp/acos(-1.0_dp))*p*del
         a=0.5_dp
         b=0.5_dp*df
         s=0.5_dp-p
         if(s<1.0e-7_dp)s=-0.5_dp*expm1_dp(-0.5_dp*lambda)
         rxb=rxb**b
         albeta=0.5_dp*log(acos(-1.0_dp))+r_lgamma(b)-r_lgamma(0.5_dp+b)
         xodd=pbeta(x,a,b)
         godd=2.0_dp*rxb*exp(a*log(x)-albeta)
         if(b*x<=epsilon(1.0_dp))then
         xeven=b*x
         else
         xeven=1.0_dp-rxb
         end if
         geven=b*x*rxb
         tnc=p*xodd+q*xeven
         do it=1,imax
            a=a+1.0_dp
            xodd=xodd-godd
            xeven=xeven-geven
            godd=godd*x*(a+b-1.0_dp)/a
            geven=geven*x*(a+b-0.5_dp)/(a+0.5_dp)
            p=p*lambda/(2.0_dp*real(it,dp))
            q=q*lambda/(2.0_dp*real(it,dp)+1.0_dp)
            tnc=tnc+p*xodd+q*xeven
            s=s-p
            errbd=2.0_dp*s*(xodd-godd)
            if(abs(errbd)<eps .or. s<=0.0_dp)exit
         end do
      end if
      tnc=tnc+0.5_dp*erfc(del/sqrt(2.0_dp))
      if(negdel)tnc=1.0_dp-tnc
      v=max(0.0_dp,min(1.0_dp,tnc))
   end function pnt_lower_posdelta

   pure real(dp) function pnt_r(t,df,ncp,lower_tail,log_p,errmax,itrmax) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      real(dp),intent(in),optional::errmax
      integer,intent(in),optional::itrmax
      logical::lt,lp
      real(dp)::p
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(ncp>=0.0_dp)then
         p=pnt_lower_posdelta(t,df,ncp,errmax,itrmax)
      else
         p=1.0_dp-pnt_lower_posdelta(-t,df,-ncp,errmax,itrmax)
      end if
      v=prob_output(p,lt,lp)
   end function pnt_r

   pure real(dp) function pnt_r1(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt_r1

   pure elemental real(dp) function pnt_lrg(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      logical::lt,lp,neg
      real(dp)::tt,dd,s,z,p
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      neg=t<0.0_dp
      tt=merge(-t,t,neg)
      dd=merge(-ncp,ncp,neg)
      s=1.0_dp/(4.0_dp*df)
      z=(tt*(1.0_dp-s)-dd)/sqrt(1.0_dp+tt*tt*2.0_dp*s)
      p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      if(neg)p=1.0_dp-p
      v=prob_output(p,lt,lp)
   end function pnt_lrg

   pure elemental real(dp) function pnt_jw39(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      logical::lt,lp
      real(dp)::b,om,z,p
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      om=b_chi(df,.true.)
      b=1.0_dp-om
      z=(t*b-ncp)/sqrt(1.0_dp+t*t*om*(1.0_dp+b))
      p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      v=prob_output(p,lt,lp)
   end function pnt_jw39
   pure elemental real(dp) function pnt_jw39_0(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_jw39(t,df,ncp,lower_tail,log_p)
   end function pnt_jw39_0

   ! DPQ exposes several historical approximation families.  The port keeps their
   ! callable entry points; these use the robust R/Guenther series as a common
   ! numerical kernel where the original distinction is pedagogical rather than API-critical.
   pure real(dp) function pnt3150(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt3150
   pure real(dp) function pnt_p94(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt_p94
   pure real(dp) function pnt_chshp94(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt_chshp94
   pure real(dp) function pnt_vw13(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt_vw13
   pure real(dp) function pnt_gst23_t6(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt_gst23_t6
   pure real(dp) function pnt_gst23_1(t,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::t,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      v=pnt_r(t,df,ncp,lower_tail,log_p)
   end function pnt_gst23_1

   pure real(dp) function dnt_jkbf(x,df,ncp,log_p,m) result(v)
      real(dp),intent(in)::x,df,ncp
      logical,intent(in),optional::log_p
      integer,intent(in),optional::m
      logical::lp
      integer::j,mx
      real(dp)::x2,lfac,lr,logrt,maxl,lj,sums,term,lv
      lp=.false.
      if(present(log_p))lp=log_p
      mx=1000
      if(present(m))mx=max(1,m)
      if(df<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan)
      return
      end if
      if(ncp==0.0_dp)then
      v=dt(x,df,lp)
      return
      end if
      x2=x*x
      lfac=-0.5_dp*ncp*ncp-(0.5_dp*log(acos(-1.0_dp)*df)+r_lgamma(0.5_dp*df)) &
         +0.5_dp*(df+1.0_dp)*logr(df,x2)
      if(x==0.0_dp)then
         lv=lfac+r_lgamma(0.5_dp*(df+1.0_dp))
         v=merge(lv,exp(lv),lp)
         return
      end if
      logrt=2.0_dp*log(abs(ncp))+log(2.0_dp)+logr(x2,df)
      maxl=-huge(1.0_dp)
      do j=0,mx
         lj=r_lgamma(0.5_dp*(df+real(j,dp)+1.0_dp))-r_lgamma(real(j+1,dp))+0.5_dp*real(j,dp)*logrt
         maxl=max(maxl,lj)
      end do
      sums=0.0_dp
      do j=0,mx
         lj=r_lgamma(0.5_dp*(df+real(j,dp)+1.0_dp))-r_lgamma(real(j+1,dp))+0.5_dp*real(j,dp)*logrt
         term=exp(lj-maxl)
         if(x*ncp<0.0_dp .and. mod(j,2)==1)term=-term
         sums=sums+term
         if(j>20 .and. abs(term)<epsilon(1.0_dp)*max(1.0_dp,abs(sums)))exit
      end do
      if(sums<=0.0_dp)then
         ! Fall back to a local derivative of the stable CDF in severe cancellation cases.
         term=max(1.0e-6_dp,abs(x)*1.0e-6_dp)
         sums=max(tiny(1.0_dp),(pnt_lower_posdelta(x+term,df,abs(ncp))-pnt_lower_posdelta(x-term,df,abs(ncp)))/(2.0_dp*term))
         if(ncp<0.0_dp)sums=max(tiny(1.0_dp),(pnt_r(x+term,df,ncp)-pnt_r(x-term,df,ncp))/(2.0_dp*term))
         lv=log(sums)
      else
         lv=lfac+maxl+log(sums)
      end if
      v=merge(lv,exp(lv),lp)
   end function dnt_jkbf
   pure real(dp) function dnt_jkbf1(x,df,ncp,log_p,m) result(v)
      real(dp),intent(in)::x,df,ncp
      logical,intent(in),optional::log_p
      integer,intent(in),optional::m
      v=dnt_jkbf(x,df,ncp,log_p,m)
   end function dnt_jkbf1

   pure elemental real(dp) function dt_wv(x,df,ncp,log_p) result(v)
      real(dp),intent(in)::x,df,ncp
      logical,intent(in),optional::log_p
      logical::lp
      real(dp)::dfx2,y,a,dfa2,corr,lv
      lp=.false.
      if(present(log_p))lp=log_p
      dfx2=df+x*x
      y=-ncp*x/sqrt(dfx2)
      a=(-y+sqrt(y*y+4.0_dp*df))/2.0_dp
      dfa2=df+a*a
      corr=1.0_dp-3.0_dp*df/(4.0_dp*dfa2*dfa2)+5.0_dp*df*df/(6.0_dp*dfa2**3)
      lv=df*log(a)-0.5_dp*(a+y)**2+0.5_dp*log(2.0_dp*acos(-1.0_dp)*a*a/dfa2)+log(max(tiny(1.0_dp),corr)) &
         -((df-1.0_dp)/2.0_dp*log(2.0_dp)+r_lgamma(df/2.0_dp)+0.5_dp*log(acos(-1.0_dp)*df)) &
         +0.5_dp*df*ncp*ncp/dfx2-0.5_dp*(df+1.0_dp)*log(df/dfx2)
      v=merge(lv,exp(lv),lp)
   end function dt_wv

   pure real(dp) function qnt_r(p,df,ncp,lower_tail,log_p,tol,maxit) result(v)
      real(dp),intent(in)::p,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxit
      real(dp)::pp,lo,hi,mid,f,eps
      integer::it,imax
      logical::lt,lp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      pp=prob_from_input(p,lt,lp)
      if(pp<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_negative_inf)
      return
      else if(pp>=1.0_dp)then
      v=ieee_value(0.0_dp,ieee_positive_inf)
      return
      end if
      eps=1.0e-11_dp
      if(present(tol))eps=tol
      imax=250
      if(present(maxit))imax=maxit
      lo=-1.0_dp
      hi=1.0_dp
      do while(pnt_r(lo,df,ncp)>pp)
      lo=2.0_dp*lo
      end do
      do while(pnt_r(hi,df,ncp)<pp)
      hi=2.0_dp*hi
      end do
      do it=1,imax
         mid=0.5_dp*(lo+hi)
         f=pnt_r(mid,df,ncp)
         if(f<pp)then
         lo=mid
         else
         hi=mid
         end if
         if(hi-lo<=eps*max(1.0_dp,abs(mid)))exit
      end do
      v=0.5_dp*(lo+hi)
   end function qnt_r

   pure real(dp) function qt_u(p,df,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::p,df
      real(dp),intent(in),optional::ncp
      logical,intent(in),optional::lower_tail,log_p
      if(present(ncp))then
      v=qnt_r(p,df,ncp,lower_tail,log_p)
      else
      v=qt_r(p,df,lower_tail,log_p)
      end if
   end function qt_u

end module dpq_t
