! Extended continuous families from gamlss.dist 6.1-1.
! Original package GPL-2 | GPL-3. Translation SPDX-License-Identifier: GPL-3.0-only
module gamlss_continuous_v02
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp, pi
   use gamlss_special, only : normal_cdf, normal_quantile, regularized_gamma_p, regularized_beta, gamma_quantile
   use gamlss_base, only : dbeta_v, qbeta_v
   use gamlss_student_t, only : student_t_pdf, student_t_cdf, student_t_quantile
   use gamlss_v02_numerics, only : adaptive_integral, bisection_root, log_bessel_k
   implicit none
   private
   public :: dGIG,pGIG,qGIG,rGIG
   public :: dSHASHo,pSHASHo,qSHASHo,rSHASHo,dSHASH,pSHASH,qSHASH,rSHASH
   public :: dSIMPLEX,pSIMPLEX,qSIMPLEX,rSIMPLEX
   public :: dSEP,pSEP,qSEP,rSEP,dSEP1,pSEP1,qSEP1,rSEP1
   public :: dSEP2,pSEP2,qSEP2,rSEP2
   public :: dST1,pST1,qST1,rST1,dST2,pST2,qST2,rST2
   public :: dST3,pST3,qST3,rST3,dST4,pST4,qST4,rST4,dST5,pST5,qST5,rST5
   public :: dSEP3,pSEP3,qSEP3,rSEP3,dSEP4,pSEP4,qSEP4,rSEP4,dNET,pNET,qNET,rNET

   integer,parameter :: CB_GIG=1,CB_SHASH=2,CB_SIMPLEX=3,CB_SEP=4,CB_SEP1=5,CB_ST1=6,CB_ST2=7,CB_NET=8
   integer,save :: pdf_callback_family=0,cdf_callback_family=0
   real(dp),save :: pdf_callback_par(4)=0.0_dp,cdf_callback_par(4)=0.0_dp
contains

   elemental real(dp) function nanv() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nanv

   elemental real(dp) function infv() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function infv

   elemental logical function want_log(flag) result(v)
      logical, intent(in), optional :: flag
      v=.false.
      if (present(flag)) v=flag
   end function want_log

   elemental logical function lower(flag) result(v)
      logical, intent(in), optional :: flag
      v=.true.
      if (present(flag)) v=flag
   end function lower

   elemental real(dp) function finish_prob(p,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail,log_p
      v=max(0.0_dp,min(1.0_dp,p))
      if (.not.lower(lower_tail)) v=1.0_dp-v
      if (present(log_p)) then
         if (log_p) then
            if (v<=0.0_dp) then
               v=-infv()
            else
               v=log(v)
            end if
         end if
      end if
   end function finish_prob

   elemental real(dp) function input_prob(p,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail,log_p
      v=p
      if (present(log_p)) then
         if (log_p) v=exp(v)
      end if
      if (.not.lower(lower_tail)) v=1.0_dp-v
   end function input_prob

   real(dp) function v02_pdf_callback(x) result(y)
      real(dp),intent(in)::x
      select case(pdf_callback_family)
      case(CB_GIG)
         if(x<=0.0_dp)then;y=0.0_dp;else;y=dGIG(x,pdf_callback_par(1),pdf_callback_par(2),pdf_callback_par(3));end if
      case(CB_SIMPLEX)
         y=dSIMPLEX(x,pdf_callback_par(1),pdf_callback_par(2))
      case(CB_SEP)
         y=dSEP(x,0.0_dp,1.0_dp,pdf_callback_par(3),pdf_callback_par(4))
      case(CB_SEP1)
         y=dSEP1(x,0.0_dp,1.0_dp,pdf_callback_par(3),pdf_callback_par(4))
      case(CB_ST1)
         y=dST1(x,0.0_dp,1.0_dp,pdf_callback_par(3),pdf_callback_par(4))
      case(CB_ST2)
         y=dST2(x,0.0_dp,1.0_dp,pdf_callback_par(3),pdf_callback_par(4))
      case default
         y=0.0_dp
      end select
   end function v02_pdf_callback

   real(dp) function v02_cdf_callback(x) result(y)
      real(dp),intent(in)::x
      select case(cdf_callback_family)
      case(CB_GIG)
         y=pGIG(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3))
      case(CB_SHASH)
         y=pSHASH(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3),cdf_callback_par(4))
      case(CB_SIMPLEX)
         y=pSIMPLEX(x,cdf_callback_par(1),cdf_callback_par(2))
      case(CB_SEP)
         y=pSEP(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3),cdf_callback_par(4))
      case(CB_SEP1)
         y=pSEP1(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3),cdf_callback_par(4))
      case(CB_ST1)
         y=pST1(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3),cdf_callback_par(4))
      case(CB_ST2)
         y=pST2(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3),cdf_callback_par(4))
      case(CB_NET)
         y=pNET(x,cdf_callback_par(1),cdf_callback_par(2),cdf_callback_par(3),cdf_callback_par(4))
      case default
         y=0.0_dp
      end select
   end function v02_cdf_callback

   real(dp) function dGIG(x,mu,sigma,nu,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu
      logical, intent(in), optional :: log_density
      real(dp) :: c,lk0,lk1,ld
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
      else if (x<=0.0_dp) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         lk0=log_bessel_k(1.0_dp/(sigma*sigma),nu)
         lk1=log_bessel_k(1.0_dp/(sigma*sigma),nu+1.0_dp)
         c=exp(lk1-lk0)
         ld=nu*log(c)-nu*log(mu)+(nu-1.0_dp)*log(x)-log(2.0_dp)-lk0
         ld=ld-(c*x/mu+mu/(c*x))/(2.0_dp*sigma*sigma)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dGIG

   real(dp) function pGIG(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: raw
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
         return
      end if
      if (q<=0.0_dp) then
         raw=0.0_dp
      else
         pdf_callback_family=CB_GIG;pdf_callback_par=[mu,sigma,nu,0.0_dp]
         raw=adaptive_integral(v02_pdf_callback,0.0_dp,q,2.0e-9_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pGIG

   real(dp) function qGIG(p,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,hi
      pp=input_prob(p,lower_tail,log_p)
      if (pp<0.0_dp .or. pp>1.0_dp .or. mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
         return
      end if
      if (pp==0.0_dp) then
         v=0.0_dp
         return
      else if (pp==1.0_dp) then
         v=infv()
         return
      end if
      cdf_callback_family=CB_GIG;cdf_callback_par=[mu,sigma,nu,0.0_dp]
      hi=max(mu,1.0e-6_dp)
      do while (v02_cdf_callback(hi)<pp .and. hi<1.0e12_dp)
         hi=2.0_dp*hi
      end do
      v=bisection_root(v02_cdf_callback,pp,0.0_dp,hi,2.0e-8_dp,100)
   end function qGIG

   real(dp) function rGIG(mu,sigma,nu) result(v)
      real(dp), intent(in) :: mu,sigma,nu
      real(dp) :: u
      call random_number(u)
      v=qGIG(u,mu,sigma,nu)
   end function rGIG

   elemental real(dp) function dSHASHo(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,a,r,c,ld
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      a=tau*asinh(z)-nu
      r=sinh(a)
      c=cosh(a)
      ld=-log(sigma)+log(tau)-0.5_dp*log(2.0_dp*pi)-0.5_dp*log(1.0_dp+z*z)
      ld=ld+log(c)-0.5_dp*r*r
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSHASHo

   elemental real(dp) function pSHASHo(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,r
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      r=sinh(tau*asinh(z)-nu)
      v=finish_prob(normal_cdf(r),lower_tail,log_p)
   end function pSHASHo

   elemental real(dp) function qSHASHo(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,z
      pp=input_prob(p,lower_tail,log_p)
      if (sigma<=0.0_dp .or. tau<=0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
         v=nanv()
         return
      end if
      z=normal_quantile(pp)
      v=mu+sigma*sinh((asinh(z)+nu)/tau)
   end function qSHASHo

   real(dp) function rSHASHo(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qSHASHo(u,mu,sigma,nu,tau)
   end function rSHASHo

   elemental real(dp) function dSHASH(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,a,r,c,ld
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      a=asinh(z)
      r=0.5_dp*(exp(tau*a)-exp(-nu*a))
      c=0.5_dp*(tau*exp(tau*a)+nu*exp(-nu*a))
      ld=-log(sigma)-0.5_dp*log(2.0_dp*pi)-0.5_dp*log(1.0_dp+z*z)
      ld=ld+log(c)-0.5_dp*r*r
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSHASH

   elemental real(dp) function pSHASH(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,a,r
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      a=asinh(z)
      r=0.5_dp*(exp(tau*a)-exp(-nu*a))
      v=finish_prob(normal_cdf(r),lower_tail,log_p)
   end function pSHASH

   real(dp) function qSHASH(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,lo,hi,width
      pp=input_prob(p,lower_tail,log_p)
      if (pp<=0.0_dp .or. pp>=1.0_dp .or. sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      cdf_callback_family=CB_SHASH;cdf_callback_par=[mu,sigma,nu,tau]
      width=sigma
      lo=mu-width
      hi=mu+width
      do while (v02_cdf_callback(lo)>pp)
         width=2.0_dp*width
         lo=mu-width
      end do
      width=sigma
      do while (v02_cdf_callback(hi)<pp)
         width=2.0_dp*width
         hi=mu+width
      end do
      v=bisection_root(v02_cdf_callback,pp,lo,hi,1.0e-9_dp,100)
   end function qSHASH

   real(dp) function rSHASH(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qSHASH(u,mu,sigma,nu,tau)
   end function rSHASH

   elemental real(dp) function dSIMPLEX(x,mu,sigma,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma
      logical, intent(in), optional :: log_density
      real(dp) :: ld
      if (mu<=0.0_dp .or. mu>=1.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
      else if (x<=0.0_dp .or. x>=1.0_dp) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         ld=-0.5_dp*((x-mu)/(mu*(1.0_dp-mu)))**2/(x*(1.0_dp-x)*sigma*sigma)
         ld=ld-0.5_dp*(log(2.0_dp*pi*sigma*sigma)+3.0_dp*(log(x)+log(1.0_dp-x)))
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dSIMPLEX

   real(dp) function pSIMPLEX(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: raw
      if (mu<=0.0_dp .or. mu>=1.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
         return
      end if
      if (q<=0.0_dp) then
         raw=0.0_dp
      else if (q>=1.0_dp) then
         raw=1.0_dp
      else
         pdf_callback_family=CB_SIMPLEX;pdf_callback_par=[mu,sigma,0.0_dp,0.0_dp]
         raw=adaptive_integral(v02_pdf_callback,0.0_dp,q,2.0e-9_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pSIMPLEX

   real(dp) function qSIMPLEX(p,mu,sigma,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp
      pp=input_prob(p,lower_tail,log_p)
      if (pp<0.0_dp .or. pp>1.0_dp) then
         v=nanv()
      else if (pp==0.0_dp) then
         v=0.0_dp
      else if (pp==1.0_dp) then
         v=1.0_dp
      else
         cdf_callback_family=CB_SIMPLEX;cdf_callback_par=[mu,sigma,0.0_dp,0.0_dp]
         v=bisection_root(v02_cdf_callback,pp,1.0e-10_dp,1.0_dp-1.0e-10_dp,1.0e-9_dp,100)
      end if
   end function qSIMPLEX

   real(dp) function rSIMPLEX(mu,sigma) result(v)
      real(dp), intent(in) :: mu,sigma
      real(dp) :: u
      call random_number(u)
      v=qSIMPLEX(u,mu,sigma)
   end function rSIMPLEX

   elemental real(dp) function powerexp_pdf(z,tau) result(v)
      real(dp), intent(in) :: z,tau
      real(dp) :: ld
      if (tau<=0.0_dp) then
         v=nanv()
      else
         ld=(1.0_dp-1.0_dp/tau)*log(tau)-abs(z)**tau/tau-log_gamma(1.0_dp/tau)-log(2.0_dp)
         v=exp(ld)
      end if
   end function powerexp_pdf

   elemental real(dp) function powerexp_cdf(z,tau) result(v)
      real(dp), intent(in) :: z,tau
      real(dp) :: s
      if (tau<=0.0_dp) then
         v=nanv()
         return
      end if
      if (z==0.0_dp) then
         v=0.5_dp
      else
         s=abs(z)**tau/tau
         v=0.5_dp*(1.0_dp+sign(1.0_dp,z)*regularized_gamma_p(1.0_dp/tau,s))
      end if
   end function powerexp_cdf

   elemental real(dp) function dSEP(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,w,ld,pn
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      w=sign(1.0_dp,z)*abs(z)**(0.5_dp*tau)*nu*sqrt(2.0_dp/tau)
      pn=max(normal_cdf(w),tiny(1.0_dp))
      ld=log(pn)-abs(z)**tau/tau-log(sigma)-log_gamma(1.0_dp/tau)
      ld=ld-(1.0_dp/tau-1.0_dp)*log(tau)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSEP

   real(dp) function pSEP(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,lo,raw
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      lo=-max(40.0_dp,25.0_dp**(1.0_dp/tau))
      if (z<=lo) then
         raw=0.0_dp
      else
         pdf_callback_family=CB_SEP;pdf_callback_par=[0.0_dp,1.0_dp,nu,tau]
         raw=adaptive_integral(v02_pdf_callback,lo,z,2.0e-8_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pSEP

   real(dp) function qSEP(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,lo,hi,width
      pp=input_prob(p,lower_tail,log_p)
      if (pp<=0.0_dp .or. pp>=1.0_dp .or. sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      cdf_callback_family=CB_SEP;cdf_callback_par=[mu,sigma,nu,tau]
      width=sigma
      lo=mu-width
      hi=mu+width
      do while (v02_cdf_callback(lo)>pp)
         width=2.0_dp*width
         lo=mu-width
      end do
      width=sigma
      do while (v02_cdf_callback(hi)<pp)
         width=2.0_dp*width
         hi=mu+width
      end do
      v=bisection_root(v02_cdf_callback,pp,lo,hi,2.0e-8_dp,100)
   end function qSEP

   real(dp) function rSEP(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qSEP(u,mu,sigma,nu,tau)
   end function rSEP

   elemental real(dp) function dSEP1(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,base,cdf,ld
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      base=powerexp_pdf(z,tau)
      cdf=max(powerexp_cdf(nu*z,tau),tiny(1.0_dp))
      ld=log(2.0_dp)+log(base)+log(cdf)-log(sigma)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSEP1

   real(dp) function pSEP1(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,lo,raw
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      lo=-max(40.0_dp,25.0_dp**(1.0_dp/tau))
      if (z<=lo) then
         raw=0.0_dp
      else
         pdf_callback_family=CB_SEP1;pdf_callback_par=[0.0_dp,1.0_dp,nu,tau]
         raw=adaptive_integral(v02_pdf_callback,lo,z,2.0e-8_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pSEP1

   real(dp) function qSEP1(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,lo,hi,width
      pp=input_prob(p,lower_tail,log_p)
      if (pp<=0.0_dp .or. pp>=1.0_dp .or. sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      cdf_callback_family=CB_SEP1;cdf_callback_par=[mu,sigma,nu,tau]
      width=sigma
      lo=mu-width
      hi=mu+width
      do while (v02_cdf_callback(lo)>pp)
         width=2.0_dp*width
         lo=mu-width
      end do
      width=sigma
      do while (v02_cdf_callback(hi)<pp)
         width=2.0_dp*width
         hi=mu+width
      end do
      v=bisection_root(v02_cdf_callback,pp,lo,hi,2.0e-8_dp,100)
   end function qSEP1

   real(dp) function rSEP1(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qSEP1(u,mu,sigma,nu,tau)
   end function rSEP1

   elemental real(dp) function dSEP2(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      v=dSEP(x,mu,sigma,nu,tau,log_density)
   end function dSEP2

   real(dp) function pSEP2(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      v=pSEP(q,mu,sigma,nu,tau,lower_tail,log_p)
   end function pSEP2

   real(dp) function qSEP2(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      v=qSEP(p,mu,sigma,nu,tau,lower_tail,log_p)
   end function qSEP2

   real(dp) function rSEP2(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      v=rSEP(mu,sigma,nu,tau)
   end function rSEP2

   elemental real(dp) function dST1(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,ld,cdf
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      if (tau<1.0e6_dp) then
         cdf=max(student_t_cdf(nu*z,tau),tiny(1.0_dp))
         ld=log(2.0_dp)+log(max(student_t_pdf(z,tau),tiny(1.0_dp)))+log(cdf)-log(sigma)
      else
         cdf=max(normal_cdf(nu*z),tiny(1.0_dp))
         ld=log(2.0_dp)-0.5_dp*z*z-0.5_dp*log(2.0_dp*pi)+log(cdf)-log(sigma)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function dST1

   real(dp) function pST1(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,lo,raw
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      lo=-max(50.0_dp,20.0_dp*sqrt(max(tau,1.0_dp)))
      if (z<=lo) then
         raw=0.0_dp
      else
         pdf_callback_family=CB_ST1;pdf_callback_par=[0.0_dp,1.0_dp,nu,tau]
         raw=adaptive_integral(v02_pdf_callback,lo,z,4.0e-8_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pST1

   real(dp) function qST1(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,lo,hi,width
      pp=input_prob(p,lower_tail,log_p)
      if (pp<=0.0_dp .or. pp>=1.0_dp .or. sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      cdf_callback_family=CB_ST1;cdf_callback_par=[mu,sigma,nu,tau]
      width=sigma
      lo=mu-width
      hi=mu+width
      do while (v02_cdf_callback(lo)>pp)
         width=2.0_dp*width
         lo=mu-width
      end do
      width=sigma
      do while (v02_cdf_callback(hi)<pp)
         width=2.0_dp*width
         hi=mu+width
      end do
      v=bisection_root(v02_cdf_callback,pp,lo,hi,3.0e-8_dp,100)
   end function qST1

   real(dp) function rST1(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qST1(u,mu,sigma,nu,tau)
   end function rST1

   elemental real(dp) function dST2(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,lam,w,ld,cdf
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      if (tau<1.0e6_dp) then
         lam=(tau+1.0_dp)/(tau+z*z)
         w=nu*sqrt(lam)*z
         cdf=max(student_t_cdf(w,tau+1.0_dp),tiny(1.0_dp))
         ld=log(2.0_dp)+log(max(student_t_pdf(z,tau),tiny(1.0_dp)))+log(cdf)-log(sigma)
      else
         w=nu*z
         cdf=max(normal_cdf(w),tiny(1.0_dp))
         ld=log(2.0_dp)-0.5_dp*z*z-0.5_dp*log(2.0_dp*pi)+log(cdf)-log(sigma)
      end if
      v=merge(ld,exp(ld),want_log(log_density))
   end function dST2

   real(dp) function pST2(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,lo,raw
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      lo=-max(50.0_dp,20.0_dp*sqrt(max(tau,1.0_dp)))
      if (z<=lo) then
         raw=0.0_dp
      else
         pdf_callback_family=CB_ST2;pdf_callback_par=[0.0_dp,1.0_dp,nu,tau]
         raw=adaptive_integral(v02_pdf_callback,lo,z,4.0e-8_dp,22)
      end if
      v=finish_prob(raw,lower_tail,log_p)
   end function pST2

   real(dp) function qST2(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,lo,hi,width
      pp=input_prob(p,lower_tail,log_p)
      if (pp<=0.0_dp .or. pp>=1.0_dp .or. sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      cdf_callback_family=CB_ST2;cdf_callback_par=[mu,sigma,nu,tau]
      width=sigma
      lo=mu-width
      hi=mu+width
      do while (v02_cdf_callback(lo)>pp)
         width=2.0_dp*width
         lo=mu-width
      end do
      width=sigma
      do while (v02_cdf_callback(hi)<pp)
         width=2.0_dp*width
         hi=mu+width
      end do
      v=bisection_root(v02_cdf_callback,pp,lo,hi,3.0e-8_dp,100)
   end function qST2

   real(dp) function rST2(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qST2(u,mu,sigma,nu,tau)
   end function rST2

   elemental real(dp) function dST5(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,ve,lam,a,b,root,ld
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      ve=2.0_dp/tau
      lam=2.0_dp*nu/(tau*sqrt(2.0_dp*tau+nu*nu))
      a=0.5_dp*(ve+lam)
      b=0.5_dp*(ve-lam)
      if (a<=0.0_dp .or. b<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      root=sqrt(a+b+z*z)
      ld=(a+0.5_dp)*log(1.0_dp+z/root)+(b+0.5_dp)*log(1.0_dp-z/root)
      ld=ld-(a+b-1.0_dp)*log(2.0_dp)-0.5_dp*log(a+b)
      ld=ld-(log_gamma(a)+log_gamma(b)-log_gamma(a+b))-log(sigma)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dST5

   elemental real(dp) function pST5(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,ve,lam,a,b,alpha
      if (sigma<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      ve=2.0_dp/tau
      lam=2.0_dp*nu/(tau*sqrt(2.0_dp*tau+nu*nu))
      a=0.5_dp*(ve+lam)
      b=0.5_dp*(ve-lam)
      if (a<=0.0_dp .or. b<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      alpha=0.5_dp*(1.0_dp+z/sqrt(a+b+z*z))
      v=finish_prob(regularized_beta(alpha,a,b),lower_tail,log_p)
   end function pST5

   real(dp) function qST5(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,ve,lam,a,b,alpha,z
      pp=input_prob(p,lower_tail,log_p)
      if (sigma<=0.0_dp .or. tau<=0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
         v=nanv()
         return
      end if
      ve=2.0_dp/tau
      lam=2.0_dp*nu/(tau*sqrt(2.0_dp*tau+nu*nu))
      a=0.5_dp*(ve+lam)
      b=0.5_dp*(ve-lam)
      if (a<=0.0_dp .or. b<=0.0_dp) then
         v=nanv()
         return
      end if
      alpha=qbeta_v(pp,a,b)
      z=sqrt(a+b)*(2.0_dp*alpha-1.0_dp)/(2.0_dp*sqrt(alpha*(1.0_dp-alpha)))
      v=mu+sigma*z
   end function qST5

   real(dp) function rST5(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qST5(u,mu,sigma,nu,tau)
   end function rST5

   elemental real(dp) function dNET(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: k1,k2,c1,c2,c3,ct,t,d,ld
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<nu .or. nu*tau<=1.0_dp) then
         v=nanv()
         return
      end if
      k1=nu
      k2=tau
      c1=(1.0_dp-2.0_dp*normal_cdf(-k1))*sqrt(2.0_dp*pi)
      c2=2.0_dp*exp(-0.5_dp*k1*k1)/k1
      c3=2.0_dp*exp(-k1*k2+0.5_dp*k1*k1)/((k1*k2-1.0_dp)*k1)
      ct=1.0_dp/(c1+c2+c3)
      t=(x-mu)/sigma
      if (abs(t)<=k1) then
         d=-0.5_dp*t*t
      else if (abs(t)<=k2) then
         d=-k1*abs(t)+0.5_dp*k1*k1
      else
         d=-k1*k2*log(abs(t)/k2)-k1*k2+0.5_dp*k1*k1
      end if
      ld=log(ct)-log(sigma)+d
      v=merge(ld,exp(ld),want_log(log_density))
   end function dNET

   elemental real(dp) function pNET(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: k1,k2,c1,c2,c3,ct,t,b,pneg,abst
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<nu .or. nu*tau<=1.0_dp) then
         v=nanv()
         return
      end if
      k1=nu
      k2=tau
      c1=(1.0_dp-2.0_dp*normal_cdf(-k1))*sqrt(2.0_dp*pi)
      c2=2.0_dp*exp(-0.5_dp*k1*k1)/k1
      c3=2.0_dp*exp(-k1*k2+0.5_dp*k1*k1)/((k1*k2-1.0_dp)*k1)
      ct=1.0_dp/(c1+c2+c3)
      t=(q-mu)/sigma
      abst=abs(t)
      b=(ct/(k1*(k1*k2-1.0_dp)))*exp(-k1*k2+0.5_dp*k1*k1)
      if (abst>=k2) then
         pneg=ct*k2**(k1*k2)*exp(-k1*k2+0.5_dp*k1*k1)
         pneg=pneg*abst**(-k1*k2+1.0_dp)/(k1*k2-1.0_dp)
      else if (abst>=k1) then
         pneg=b+(ct/k1)*exp(-k1*abst+0.5_dp*k1*k1)
      else
         pneg=b+(ct/k1)*exp(-0.5_dp*k1*k1)
         pneg=pneg+ct*sqrt(2.0_dp*pi)*(normal_cdf(-abst)-normal_cdf(-k1))
      end if
      if (t<=0.0_dp) then
         v=finish_prob(pneg,lower_tail,log_p)
      else
         v=finish_prob(1.0_dp-pneg,lower_tail,log_p)
      end if
   end function pNET

   real(dp) function qNET(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,lo,hi,width
      pp=input_prob(p,lower_tail,log_p)
      if (pp<=0.0_dp .or. pp>=1.0_dp .or. sigma<=0.0_dp .or. nu<=0.0_dp) then
         v=nanv()
         return
      end if
      cdf_callback_family=CB_NET;cdf_callback_par=[mu,sigma,nu,tau]
      width=sigma
      lo=mu-width
      hi=mu+width
      do while (v02_cdf_callback(lo)>pp)
         width=2.0_dp*width
         lo=mu-width
      end do
      width=sigma
      do while (v02_cdf_callback(hi)<pp)
         width=2.0_dp*width
         hi=mu+width
      end do
      v=bisection_root(v02_cdf_callback,pp,lo,hi,1.0e-10_dp,100)
   end function qNET

   real(dp) function rNET(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qNET(u,mu,sigma,nu,tau)
   end function rNET


   elemental real(dp) function dST3(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,ld
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      if (x<mu) then
         ld=log(max(student_t_pdf(nu*z,tau),tiny(1.0_dp)))
      else
         ld=log(max(student_t_pdf(z/nu,tau),tiny(1.0_dp)))
      end if
      ld=ld+log(2.0_dp*nu/(1.0_dp+nu*nu))-log(sigma)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dST3

   elemental real(dp) function pST3(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: cdf
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      if (q<mu) then
         cdf=2.0_dp*student_t_cdf(nu*(q-mu)/sigma,tau)/(1.0_dp+nu*nu)
      else
         cdf=(1.0_dp+2.0_dp*nu*nu*(student_t_cdf((q-mu)/(sigma*nu),tau)-0.5_dp)) &
            /(1.0_dp+nu*nu)
      end if
      v=finish_prob(cdf,lower_tail,log_p)
   end function pST3

   real(dp) function qST3(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,cut,arg
      pp=input_prob(p,lower_tail,log_p)
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
         v=nanv()
         return
      end if
      cut=1.0_dp/(1.0_dp+nu*nu)
      if (pp<cut) then
         arg=pp*(1.0_dp+nu*nu)/2.0_dp
         v=mu+(sigma/nu)*student_t_quantile(arg,tau)
      else
         arg=(pp*(1.0_dp+nu*nu)-1.0_dp)/(2.0_dp*nu*nu)+0.5_dp
         v=mu+sigma*nu*student_t_quantile(arg,tau)
      end if
   end function qST3

   real(dp) function rST3(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qST3(u,mu,sigma,nu,tau)
   end function rST3

   elemental real(dp) function st4_logk(df) result(v)
      real(dp), intent(in) :: df
      if (df>=1.0e6_dp) then
         v=0.5_dp*log(2.0_dp*pi)
      else
         v=log_gamma(0.5_dp)+log_gamma(0.5_dp*df)-log_gamma(0.5_dp*(df+1.0_dp)) &
            +0.5_dp*log(df)
      end if
   end function st4_logk

   elemental real(dp) function dST4(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,k,ld
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      k=exp(st4_logk(tau)-st4_logk(nu))
      if (x<mu) then
         ld=log(max(student_t_pdf(z,nu),tiny(1.0_dp)))
      else
         ld=log(max(student_t_pdf(z,tau),tiny(1.0_dp)))+log(k)
      end if
      ld=ld+log(2.0_dp/(1.0_dp+k))-log(sigma)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dST4

   elemental real(dp) function pST4(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z,k,cdf
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(q-mu)/sigma
      k=exp(st4_logk(tau)-st4_logk(nu))
      if (q<mu) then
         cdf=2.0_dp*student_t_cdf(z,nu)/(1.0_dp+k)
      else
         cdf=(1.0_dp+2.0_dp*k*(student_t_cdf(z,tau)-0.5_dp))/(1.0_dp+k)
      end if
      v=finish_prob(cdf,lower_tail,log_p)
   end function pST4

   real(dp) function qST4(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,k,cut,arg
      pp=input_prob(p,lower_tail,log_p)
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
         v=nanv()
         return
      end if
      k=exp(st4_logk(tau)-st4_logk(nu))
      cut=1.0_dp/(1.0_dp+k)
      if (pp<cut) then
         arg=pp*(1.0_dp+k)/2.0_dp
         v=mu+sigma*student_t_quantile(arg,nu)
      else
         arg=(pp*(1.0_dp+k)-1.0_dp)/(2.0_dp*k)+0.5_dp
         v=mu+sigma*student_t_quantile(arg,tau)
      end if
   end function qST4

   real(dp) function rST4(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qST4(u,mu,sigma,nu,tau)
   end function rST4

   elemental real(dp) function dSEP3(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,ld
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      z=(x-mu)/sigma
      if (x<mu) then
         ld=-0.5_dp*(nu*abs(z))**tau
      else
         ld=-0.5_dp*(abs(z)/nu)**tau
      end if
      ld=ld-log(sigma)+log(nu)-log(1.0_dp+nu*nu)-(1.0_dp/tau)*log(2.0_dp)
      ld=ld-log_gamma(1.0_dp+1.0_dp/tau)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSEP3

   elemental real(dp) function pSEP3(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: z1,z2,s1,s2,k,cdf
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      k=nu*nu
      z1=nu*(q-mu)/(sigma*2.0_dp**(1.0_dp/tau))
      z2=(q-mu)/(sigma*nu*2.0_dp**(1.0_dp/tau))
      s1=abs(z1)**tau
      s2=abs(z2)**tau
      if (q<mu) then
         cdf=(1.0_dp-regularized_gamma_p(1.0_dp/tau,s1))/(1.0_dp+k)
      else
         cdf=(1.0_dp+k*regularized_gamma_p(1.0_dp/tau,s2))/(1.0_dp+k)
      end if
      v=finish_prob(cdf,lower_tail,log_p)
   end function pSEP3

   real(dp) function qSEP3(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,k,cut,g
      pp=input_prob(p,lower_tail,log_p)
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
         v=nanv()
         return
      end if
      k=nu*nu
      cut=1.0_dp/(1.0_dp+k)
      if (pp<cut) then
         g=gamma_quantile(1.0_dp-pp*(1.0_dp+k),1.0_dp/tau,1.0_dp)
         v=mu-(sigma*2.0_dp**(1.0_dp/tau)/nu)*g**(1.0_dp/tau)
      else
         g=gamma_quantile((pp*(1.0_dp+k)-1.0_dp)/k,1.0_dp/tau,1.0_dp)
         v=mu+sigma*nu*2.0_dp**(1.0_dp/tau)*g**(1.0_dp/tau)
      end if
   end function qSEP3

   real(dp) function rSEP3(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qSEP3(u,mu,sigma,nu,tau)
   end function rSEP3

   elemental real(dp) function dSEP4(x,mu,sigma,nu,tau,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu,tau
      logical, intent(in), optional :: log_density
      real(dp) :: z,k1,k2,ld
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      k1=exp(log_gamma(1.0_dp+1.0_dp/nu))
      k2=exp(log_gamma(1.0_dp+1.0_dp/tau))
      z=(x-mu)/sigma
      if (x<mu) then
         ld=-abs(z)**nu
      else
         ld=-abs(z)**tau
      end if
      ld=ld-log(k1+k2)-log(sigma)
      v=merge(ld,exp(ld),want_log(log_density))
   end function dSEP4

   elemental real(dp) function pSEP4(q,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: lk1,lk2,k,z,s,cdf
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp) then
         v=nanv()
         return
      end if
      lk1=log_gamma(1.0_dp+1.0_dp/nu)
      lk2=log_gamma(1.0_dp+1.0_dp/tau)
      k=exp(lk2-lk1)
      z=(q-mu)/sigma
      if (q<mu) then
         s=abs(z)**nu
         cdf=(1.0_dp-regularized_gamma_p(1.0_dp/nu,s))/(1.0_dp+k)
      else
         s=abs(z)**tau
         cdf=(1.0_dp+k*regularized_gamma_p(1.0_dp/tau,s))/(1.0_dp+k)
      end if
      v=finish_prob(cdf,lower_tail,log_p)
   end function pSEP4

   real(dp) function qSEP4(p,mu,sigma,nu,tau,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu,tau
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: pp,k,cut,g,lk1,lk2
      pp=input_prob(p,lower_tail,log_p)
      if (sigma<=0.0_dp .or. nu<=0.0_dp .or. tau<=0.0_dp .or. pp<=0.0_dp .or. pp>=1.0_dp) then
         v=nanv()
         return
      end if
      lk1=log_gamma(1.0_dp+1.0_dp/nu)
      lk2=log_gamma(1.0_dp+1.0_dp/tau)
      k=exp(lk2-lk1)
      cut=1.0_dp/(1.0_dp+k)
      if (pp<cut) then
         g=gamma_quantile(1.0_dp-pp*(1.0_dp+k),1.0_dp/nu,1.0_dp)
         v=mu-sigma*g**(1.0_dp/nu)
      else
         g=gamma_quantile((pp*(1.0_dp+k)-1.0_dp)/k,1.0_dp/tau,1.0_dp)
         v=mu+sigma*g**(1.0_dp/tau)
      end if
   end function qSEP4

   real(dp) function rSEP4(mu,sigma,nu,tau) result(v)
      real(dp), intent(in) :: mu,sigma,nu,tau
      real(dp) :: u
      call random_number(u)
      v=qSEP4(u,mu,sigma,nu,tau)
   end function rSEP4

end module gamlss_continuous_v02
