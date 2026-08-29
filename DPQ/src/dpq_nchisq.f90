! Noncentral chi-square algorithms and approximations translated from DPQ.
! SPDX-License-Identifier: GPL-2.0-or-later
module dpq_nchisq
   use r_compat, only: dp, pgamma, qgamma, dnorm, qnorm, r_lgamma, besselI, dpois, ppois
   use dpq_core, only: prob_output, prob_from_input, log1p_dp
   use dpq_gamma_discrete, only: dgamma_r
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
      ieee_positive_inf
   implicit none
   private
   public :: pchisq_central, qchisq_central, dchisq_central
   public :: dnchisq_r, dnchisq_bessel, dnoncentchisq, dchisq_asym
   public :: pnchisq, pnchisq_rc, pnchisq_it, pnchisq_v, pnchisq_ss
   public :: pnchi1sq, pnchi3sq, pnchisq_patnaik, pnchisq_pearson
   public :: pnchisq_abdel_aty, pnchisq_sankaran_d, pnchisq_bolkuz
   public :: pnchisq_t93, pnchisq_t93_a, pnchisq_t93_b
   public :: qnchisq, qnchisq_patnaik, qnchisq_pearson, qnchisq_abdel_aty
   public :: qnchisq_sankaran_d, qnchisq_bolkuz, qchisq_n
   public :: qchisq_appr0, qchisq_appr1, qchisq_appr2, qchisq_appr3
   public :: qchisq_appr_cf1, qchisq_appr_cf2, qchisq_cappr2, r_pois
   public :: pnchisq_terms

contains

   pure elemental real(dp) function pchisq_central(x, df) result(v)
      real(dp), intent(in) :: x, df
      if (df <= 0.0_dp) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
      else if (x <= 0.0_dp) then
         v = 0.0_dp
      else
         v = pgamma(x,0.5_dp*df,0.5_dp)
      end if
   end function pchisq_central

   pure elemental real(dp) function qchisq_central(p, df) result(v)
      real(dp), intent(in) :: p, df
      if (df <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
      else if (p == 0.0_dp) then
         v = 0.0_dp
      else if (p == 1.0_dp) then
         v = ieee_value(0.0_dp,ieee_positive_inf)
      else
         v = qgamma(p,0.5_dp*df,0.5_dp)
      end if
   end function qchisq_central

   pure elemental real(dp) function dchisq_central(x, df, log_p) result(v)
      real(dp), intent(in) :: x, df
      logical, intent(in), optional :: log_p
      logical :: lp
      lp = .false.
      if (present(log_p)) lp = log_p
      v = dgamma_r(x,0.5_dp*df,2.0_dp,lp)
   end function dchisq_central

   pure real(dp) function poisson_weight_at(k, lam) result(w)
      integer, intent(in) :: k
      real(dp), intent(in) :: lam
      if (lam == 0.0_dp) then
         w = merge(1.0_dp,0.0_dp,k == 0)
      else
         w = exp(-lam + real(k,dp)*log(lam) - r_lgamma(real(k+1,dp)))
      end if
   end function poisson_weight_at

   pure real(dp) function pnchisq_lower(x, df, ncp, tol, maxit) result(v)
      real(dp), intent(in) :: x, df, ncp
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      real(dp) :: lam, eps, w0, w, sumv, mass
      integer :: j0, j, imax
      if (df <= 0.0_dp .or. ncp < 0.0_dp) then
         v = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (x <= 0.0_dp) then
         v = 0.0_dp
         return
      end if
      if (ncp == 0.0_dp) then
         v = pchisq_central(x,df)
         return
      end if
      eps = 8.0_dp*epsilon(1.0_dp)
      if (present(tol)) eps = max(tol,epsilon(1.0_dp))
      imax = 200000
      if (present(maxit)) imax = maxit
      lam = 0.5_dp*ncp
      j0 = max(0,int(floor(lam)))
      w0 = poisson_weight_at(j0,lam)
      sumv = w0*pchisq_central(x,df+2.0_dp*real(j0,dp))
      mass = w0
      w = w0
      do j = j0-1, 0, -1
         w = w*real(j+1,dp)/lam
         sumv = sumv + w*pchisq_central(x,df+2.0_dp*real(j,dp))
         mass = mass + w
         if (j0-j > imax) exit
      end do
      w = w0
      do j = j0+1, j0+imax
         w = w*lam/real(j,dp)
         if (w == 0.0_dp) exit
         sumv = sumv + w*pchisq_central(x,df+2.0_dp*real(j,dp))
         mass = mass + w
         if (1.0_dp-mass <= eps .and. w <= eps) exit
      end do
      v = min(1.0_dp,max(0.0_dp,sumv))
   end function pnchisq_lower

   pure real(dp) function pnchisq(x, df, ncp, lower_tail, log_p, tol, maxit) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: lower_tail, log_p
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      logical :: lt, lp
      real(dp) :: p
      lt=.true.
      lp=.false.
      if (present(lower_tail)) lt=lower_tail
      if (present(log_p)) lp=log_p
      p = pnchisq_lower(x,df,ncp,tol,maxit)
      v = prob_output(p,lt,lp)
   end function pnchisq

   pure real(dp) function pnchisq_rc(x, df, ncp, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: lower_tail, log_p
      v = pnchisq(x,df,ncp,lower_tail,log_p)
   end function pnchisq_rc

   pure real(dp) function pnchisq_it(x, df, ncp, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: lower_tail, log_p
      v = pnchisq(x,df,ncp,lower_tail,log_p)
   end function pnchisq_it

   pure real(dp) function pnchisq_v(x, df, ncp, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: lower_tail, log_p
      v = pnchisq(x,df,ncp,lower_tail,log_p)
   end function pnchisq_v

   pure real(dp) function pnchisq_ss(x, df, ncp, lower_tail, log_p) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: lower_tail, log_p
      v = pnchisq(x,df,ncp,lower_tail,log_p)
   end function pnchisq_ss

   pure real(dp) function dnchisq_r(x, df, ncp, log_p, tol, maxit) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: log_p
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxit
      logical :: lp
      real(dp) :: lam, eps, w0, w, sumv, mass
      integer :: j0, j, imax
      lp=.false.
      if(present(log_p)) lp=log_p
      if (df <= 0.0_dp .or. ncp < 0.0_dp .or. x < 0.0_dp) then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if(ncp==0.0_dp) then
         v=dchisq_central(x,df,lp)
         return
      end if
      eps=8.0_dp*epsilon(1.0_dp)
      if(present(tol)) eps=max(tol,epsilon(1.0_dp))
      imax=200000
      if(present(maxit)) imax=maxit
      lam=0.5_dp*ncp
      j0=max(0,int(floor(lam)))
      w0=poisson_weight_at(j0,lam)
      sumv=w0*dchisq_central(x,df+2.0_dp*real(j0,dp))
      mass=w0
      w=w0
      do j=j0-1,0,-1
         w=w*real(j+1,dp)/lam
         sumv=sumv+w*dchisq_central(x,df+2.0_dp*real(j,dp))
         mass=mass+w
      end do
      w=w0
      do j=j0+1,j0+imax
         w=w*lam/real(j,dp)
         if(w==0.0_dp) exit
         sumv=sumv+w*dchisq_central(x,df+2.0_dp*real(j,dp))
         mass=mass+w
         if(1.0_dp-mass<=eps .and. w<=eps) exit
      end do
      if(lp) then
         if(sumv>0.0_dp) then
         v=log(sumv)
         else
         v=-huge(1.0_dp)
         end if
      else
         v=sumv
      end if
   end function dnchisq_r

   pure real(dp) function dnoncentchisq(x, df, ncp, log_p) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: log_p
      v=dnchisq_r(x,df,ncp,log_p)
   end function dnoncentchisq

   pure real(dp) function dnchisq_bessel(x, df, ncp, log_p) result(v)
      real(dp), intent(in) :: x, df, ncp
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: nu, y, iy, lv
      lp=.false.
      if(present(log_p)) lp=log_p
      if(ncp==0.0_dp) then
         v=dchisq_central(x,df,lp)
         return
      end if
      if(x<=0.0_dp .or. df<=0.0_dp .or. ncp<0.0_dp) then
         if(x==0.0_dp .and. df>2.0_dp) then
         v=merge(-huge(1.0_dp),0.0_dp,lp)
         else
         v=ieee_value(0.0_dp,ieee_quiet_nan)
         end if
         return
      end if
      nu=0.5_dp*(df-2.0_dp)
      y=sqrt(ncp*x)
      iy=besselI(y,nu,.true.)
      if(iy<=0.0_dp) then
         v=dnchisq_r(x,df,ncp,lp)
         return
      end if
      lv=(y-0.5_dp*(ncp+x))+0.5_dp*nu*log(x/ncp)+log(0.5_dp*iy)
      v=merge(lv,exp(lv),lp)
   end function dnchisq_bessel

   pure elemental real(dp) function dchisq_asym(x,df,ncp,log_p) result(v)
      real(dp),intent(in)::x,df,ncp
      logical,intent(in),optional::log_p
      logical::lp
      real(dp)::nl,n2l,ic,lv
      lp=.false.
      if(present(log_p))lp=log_p
      nl=df+ncp
      n2l=nl+ncp
      ic=nl/n2l
      lv=dchisq_central(x*ic,nl*ic,.true.)+log(ic)
      v=merge(lv,exp(lv),lp)
   end function dchisq_asym

   pure real(dp) function pnchi1sq(x,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::x,ncp
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::sq,sl,p
      logical::lt,lp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(x<=0.0_dp) then
      p=0.0_dp
      else
      sq=sqrt(x)
      sl=sqrt(max(0.0_dp,ncp))
      p=0.5_dp*erfc(-(sq-sl)/sqrt(2.0_dp))-0.5_dp*erfc((sq+sl)/sqrt(2.0_dp))
      end if
      v=prob_output(max(0.0_dp,min(1.0_dp,p)),lt,lp)
   end function pnchi1sq

   pure real(dp) function pnchi3sq(x,ncp,lower_tail,log_p) result(v)
      real(dp),intent(in)::x,ncp
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::sq,sl,d,pv,p
      logical::lt,lp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(ncp<=epsilon(1.0_dp)) then
      p=pchisq_central(x,3.0_dp)
      else if(x<=0.0_dp) then
      p=0.0_dp
      else
         sq=sqrt(x)
         sl=sqrt(ncp)
         d=sl-sq
         pv=sl+sq
         p=0.5_dp*erfc(d/sqrt(2.0_dp))-0.5_dp*erfc(pv/sqrt(2.0_dp)) &
            +(exp(-0.5_dp*pv*pv)-exp(-0.5_dp*d*d))/(sqrt(2.0_dp*acos(-1.0_dp))*sl)
      end if
      v=prob_output(max(0.0_dp,min(1.0_dp,p)),lt,lp)
   end function pnchi3sq

   pure elemental real(dp) function pnchisq_patnaik(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::e,d,ic
      e=df+ncp
      d=e+ncp
      ic=e/d
      v=pchisq_central(x*ic,e*ic)
   end function pnchisq_patnaik

   pure elemental real(dp) function pnchisq_pearson(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::n2,n3,r
      n2=df+2.0_dp*ncp
      n3=n2+ncp
      r=n2/n3
      v=pchisq_central((x+ncp*ncp/n3)*r,n2*r*r)
   end function pnchisq_pearson

   pure elemental real(dp) function pnchisq_abdel_aty(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::r,n2,var,z
      if(x<=0.0_dp) then
      v=0.0_dp
      return
      end if
      r=df+ncp
      n2=r+ncp
      var=(2.0_dp/9.0_dp)*n2/(r*r)
      z=((x/r)**(1.0_dp/3.0_dp)-(1.0_dp-var))/sqrt(var)
      v=0.5_dp*erfc(-z/sqrt(2.0_dp))
   end function pnchisq_abdel_aty

   pure elemental real(dp) function pnchisq_sankaran_d(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::r,n2,n3,q1,h1,h,mu,sv,z
      if(x<=0.0_dp) then
      v=0.0_dp
      return
      end if
      r=df+ncp
      n2=r+ncp
      n3=n2+ncp
      q1=n2/(r*r)
      h1=-(2.0_dp/3.0_dp)*r*n3/(n2*n2)
      h=1.0_dp+h1
      mu=1.0_dp+h*h1*q1*(1.0_dp+0.5_dp*(h-2.0_dp)*(1.0_dp-3.0_dp*h)*q1)
      sv=h*sqrt(2.0_dp*q1*(1.0_dp+h1*(1.0_dp-3.0_dp*h)*q1))
      z=((x/r)**h-mu)/sv
      v=0.5_dp*erfc(-z/sqrt(2.0_dp))
   end function pnchisq_sankaran_d

   pure elemental real(dp) function pnchisq_bolkuz(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::lnu,w
      lnu=ncp/df
      w=x*(1.0_dp+lnu*(-1.0_dp+0.5_dp*lnu*(1.0_dp+x/(df+2.0_dp))))
      v=pchisq_central(max(0.0_dp,w),df)
   end function pnchisq_bolkuz

   pure elemental real(dp) function pnchisq_t93(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      v=pnchisq_sankaran_d(x,df,ncp)
   end function pnchisq_t93
   pure elemental real(dp) function pnchisq_t93_a(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      v=pnchisq_patnaik(x,df,ncp)
   end function pnchisq_t93_a
   pure elemental real(dp) function pnchisq_t93_b(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      v=pnchisq_pearson(x,df,ncp)
   end function pnchisq_t93_b

   pure real(dp) function qnchisq(p,df,ncp,lower_tail,log_p,tol,maxit) result(v)
      real(dp),intent(in)::p,df,ncp
      logical,intent(in),optional::lower_tail,log_p
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxit
      real(dp)::pp,lo,hi,mid,f,eps,meanv,sdv
      integer::it,imax
      logical::lt,lp
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      pp=prob_from_input(p,lt,lp)
      if(pp<=0.0_dp)then
      v=0.0_dp
      return
      else if(pp>=1.0_dp)then
      v=ieee_value(0.0_dp,ieee_positive_inf)
      return
      end if
      eps=1.0e-12_dp
      if(present(tol))eps=tol
      imax=200
      if(present(maxit))imax=maxit
      meanv=df+ncp
      sdv=sqrt(2.0_dp*(df+2.0_dp*ncp))
      lo=0.0_dp
      hi=max(1.0_dp,meanv+8.0_dp*sdv)
      do while(pnchisq_lower(hi,df,ncp)<pp)
      hi=2.0_dp*hi
      end do
      do it=1,imax
         mid=0.5_dp*(lo+hi)
         f=pnchisq_lower(mid,df,ncp)
         if(f<pp)then
         lo=mid
         else
         hi=mid
         end if
         if(hi-lo<=eps*max(1.0_dp,mid))exit
      end do
      v=0.5_dp*(lo+hi)
   end function qnchisq

   pure elemental real(dp) function qnchisq_patnaik(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::e,d,ic
      e=df+ncp
      d=e+ncp
      ic=e/d
      v=qchisq_central(p,e*ic)/ic
   end function qnchisq_patnaik
   pure elemental real(dp) function qnchisq_pearson(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::n2,n3,r
      n2=df+2.0_dp*ncp
      n3=n2+ncp
      r=n2/n3
      v=qchisq_central(p,n2*r*r)/r-ncp*ncp/n3
   end function qnchisq_pearson
   pure elemental real(dp) function qnchisq_abdel_aty(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::r,n2,var,z
      r=df+ncp
      n2=r+ncp
      var=(2.0_dp/9.0_dp)*n2/(r*r)
      z=qnorm(p,1.0_dp-var,sqrt(var))
      v=r*z**3
   end function qnchisq_abdel_aty
   pure elemental real(dp) function qnchisq_sankaran_d(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::r,n2,n3,q1,h1,h,mu,sv,z
      r=df+ncp
      n2=r+ncp
      n3=n2+ncp
      q1=n2/(r*r)
      h1=-(2.0_dp/3.0_dp)*r*n3/(n2*n2)
      h=1.0_dp+h1
      mu=1.0_dp+h*h1*q1*(1.0_dp+0.5_dp*(h-2.0_dp)*(1.0_dp-3.0_dp*h)*q1)
      sv=h*sqrt(2.0_dp*q1*(1.0_dp+h1*(1.0_dp-3.0_dp*h)*q1))
      z=qnorm(p,mu,sv)
      v=r*max(0.0_dp,z)**(1.0_dp/h)
   end function qnchisq_sankaran_d

   pure real(dp) function qnchisq_bolkuz(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::lo,hi,mid
      integer::it
      lo=0.0_dp
      hi=max(1.0_dp,df+ncp+10.0_dp*sqrt(2.0_dp*(df+2.0_dp*ncp)))
      do it=1,160
         mid=0.5_dp*(lo+hi)
         if(pnchisq_bolkuz(mid,df,ncp)<p)then
         lo=mid
         else
         hi=mid
         end if
      end do
      v=0.5_dp*(lo+hi)
   end function qnchisq_bolkuz

   pure real(dp) function qchisq_n(p,df,ncp) result(v)
      real(dp),intent(in)::p,df
      real(dp),intent(in),optional::ncp
      if(present(ncp))then
      v=qnchisq(p,df,ncp)
      else
      v=qchisq_central(p,df)
      end if
   end function qchisq_n

   pure subroutine pnchisq_terms(x,df,ncp,terms,weights)
      real(dp),intent(in)::x,df,ncp
      real(dp),allocatable,intent(out)::terms(:),weights(:)
      real(dp)::lam,w
      integer::j,n
      lam=0.5_dp*ncp
      n=max(20,int(lam+12.0_dp*sqrt(max(1.0_dp,lam))+20.0_dp))
      allocate(terms(0:n),weights(0:n))
      w=exp(-lam)
      weights(0)=w
      terms(0)=w*pchisq_central(x,df)
      do j=1,n
         w=w*lam/real(j,dp)
         weights(j)=w
         terms(j)=w*pchisq_central(x,df+2.0_dp*real(j,dp))
      end do
   end subroutine pnchisq_terms

   pure elemental real(dp) function qchisq_appr0(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::z
      z=qnorm(p)
      v=df+ncp+z*sqrt(2.0_dp*df+4.0_dp*ncp)
   end function qchisq_appr0

   pure elemental real(dp) function qchisq_appr1(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::z,b1
      z=qnorm(p)
      b1=1.0_dp+ncp/(ncp+df)
      v=b1*df*(1.0_dp-2.0_dp/(9.0_dp*df)+z*sqrt(2.0_dp/(9.0_dp*df)))**3
   end function qchisq_appr1

   pure elemental real(dp) function qchisq_appr2(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::z,a,b1,dfs,alp
      z=qnorm(p)
      a=df+ncp
      b1=1.0_dp+ncp/a
      dfs=a/b1
      alp=2.0_dp/(9.0_dp*dfs)
      v=a*(1.0_dp-alp+z*sqrt(alp))**3
   end function qchisq_appr2

   pure elemental real(dp) function qchisq_appr3(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::z,a,dfs,alp
      z=qnorm(p)
      a=df+ncp
      dfs=a/(1.0_dp+ncp/a)
      alp=2.0_dp/(9.0_dp*dfs)
      v=dfs*(1.0_dp-alp+z*sqrt(alp))**3
   end function qchisq_appr3

   pure elemental real(dp) function qchisq_appr_cf1(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::z,mu,sig2,gam1
      z=qnorm(p)
      mu=df+ncp
      sig2=2.0_dp*(df+2.0_dp*ncp)
      gam1=8.0_dp*(df+3.0_dp*ncp)*sig2**(-1.5_dp)
      v=mu+sqrt(sig2)*(z+gam1*(z*z-1.0_dp)/6.0_dp)
   end function qchisq_appr_cf1

   pure elemental real(dp) function qchisq_appr_cf2(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::z,mu,sig2,gam1,gam2
      z=qnorm(p)
      mu=df+ncp
      sig2=2.0_dp*(df+2.0_dp*ncp)
      gam1=8.0_dp*(df+3.0_dp*ncp)*sig2**(-1.5_dp)
      gam2=48.0_dp*(df+4.0_dp*ncp)/(sig2*sig2)
      v=mu+sqrt(sig2)*(z+gam1*(z*z-1.0_dp)/6.0_dp+gam2*z*(z*z-3.0_dp)/24.0_dp &
         -gam1*gam1*z*(2.0_dp*z*z-5.0_dp)/36.0_dp)
   end function qchisq_appr_cf2

   pure elemental real(dp) function qchisq_cappr2(p,df,ncp) result(v)
      real(dp),intent(in)::p,df,ncp
      real(dp)::c2,lnu
      c2=qchisq_central(p,df)
      lnu=ncp/df
      v=c2*(1.0_dp+lnu*(1.0_dp+0.5_dp*lnu*(1.0_dp-c2/(df+2.0_dp))))
   end function qchisq_cappr2

   pure elemental real(dp) function r_pois(i,lambda) result(v)
      real(dp),intent(in)::i,lambda
      real(dp)::den
      if(i<0.0_dp .or. lambda<0.0_dp)then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
      else if(i==0.0_dp)then
         v=ieee_value(0.0_dp,ieee_positive_inf)
      else
         den=ppois(i-1.0_dp,lambda)
         if(den>0.0_dp)then
            v=dpois(i,lambda)/den
         else
            v=exp(i*log(lambda)-r_lgamma(i+1.0_dp)-lambda-log(max(tiny(1.0_dp),den)))
         end if
      end if
   end function r_pois


end module dpq_nchisq
