! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_actuarial
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
      ieee_positive_inf
   use vgam_kinds, only : dp, pi
   use vgam_special, only : log1p_v, expm1_v, normal_quantile, &
      regularized_gamma_p, gamma_quantile, lambert_w0
   use vgam_random, only : random_gamma, random_poisson
   implicit none
   private
   public :: dgompertz, pgompertz, qgompertz, rgompertz
   public :: dmakeham, pmakeham, qmakeham, rmakeham
   public :: dperks, pperks, qperks, rperks
   public :: dlindley, plindley, qlindley, rlindley
   public :: dnakagami, pnakagami, qnakagami, rnakagami
   public :: dmaxwell, pmaxwell, qmaxwell, rmaxwell
   public :: dbenini, pbenini, qbenini, rbenini
   public :: dlevy, plevy, qlevy, rlevy, dskellam, rskellam
contains
   elemental real(dp) function nanv() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function
   elemental real(dp) function infv() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function

   elemental real(dp) function dgompertz(x,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in),optional::log_density
      real(dp)::ld
      logical::lg
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
         v=nanv()
      return
      end if
      if(x<0.0_dp)then
         ld=-infv()
      else
         ld=log(shape)+x*scale-(shape/scale)*expm1_v(x*scale)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pgompertz(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
         v=nanv()
      else if(q<=0.0_dp)then
         v=0.0_dp
      else
         v=-expm1_v((-shape/scale)*expm1_v(scale*q))
      end if
   end function
   elemental real(dp) function qgompertz(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
         v=nanv()
      else if(p==1.0_dp)then
         v=infv()
      else
         v=log1p_v((-scale/shape)*log1p_v(-p))/scale
      end if
   end function
   real(dp) function rgompertz(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      real(dp)::u
      call random_number(u)
      v=qgompertz(u,scale,shape)
   end function

   elemental real(dp) function dmakeham(x,scale,shape,epsilon,log_density) result(v)
      real(dp),intent(in)::x,scale,shape,epsilon
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.epsilon<0.0_dp)then
         v=nanv()
      return
      end if
      if(x<0.0_dp)then
         ld=-infv()
      else
         ld=log(epsilon*exp(-x*scale)+shape)+x*(scale-epsilon)- &
            (shape/scale)*expm1_v(x*scale)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pmakeham(q,scale,shape,epsilon) result(v)
      real(dp),intent(in)::q,scale,shape,epsilon
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.epsilon<0.0_dp)then
         v=nanv()
      else if(q<=0.0_dp)then
         v=0.0_dp
      else
         v=-expm1_v(-q*epsilon-(shape/scale)*expm1_v(scale*q))
      end if
   end function
   elemental real(dp) function qmakeham(p,scale,shape,epsilon) result(v)
      real(dp),intent(in)::p,scale,shape,epsilon
      real(dp)::lp,arg
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.epsilon<0.0_dp.or. &
         p<0.0_dp.or.p>1.0_dp)then
         v=nanv()
      return
      end if
      if(epsilon==0.0_dp)then
         v=qgompertz(p,scale,shape)
      return
      end if
      if(p==1.0_dp)then
      v=infv()
      return
      end if
      lp=log1p_v(-p)
      arg=(shape/epsilon)*exp(shape/epsilon)* &
          exp((-scale/epsilon)*lp)
      v=shape/(scale*epsilon)-lp/epsilon-lambert_w0(arg)/scale
   end function
   real(dp) function rmakeham(scale,shape,epsilon) result(v)
      real(dp),intent(in)::scale,shape,epsilon
      real(dp)::u
      call random_number(u)
      v=qmakeham(u,scale,shape,epsilon)
   end function

   elemental real(dp) function dperks(x,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=nanv()
      return
      end if
      if(x<0.0_dp)then
         ld=-infv()
      else
         ld=log(shape)-x+log1p_v(shape)/scale- &
            (1.0_dp+1.0_dp/scale)*log(shape+exp(-x*scale))
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pperks(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      real(dp)::logs
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
         v=nanv()
      else if(q<=0.0_dp)then
         v=0.0_dp
      else
         logs=-q+(log1p_v(shape)-log(shape+exp(-q*scale)))/scale
         v=-expm1_v(logs)
      end if
   end function
   elemental real(dp) function qperks(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      real(dp)::tmp,onem
      if(scale<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
         v=nanv()
      else if(p==1.0_dp)then
         v=infv()
      else
         tmp=scale*log1p_v(-p)
      onem=exp(tmp)
         v=(log1p_v(shape-onem)-log(shape)-tmp)/scale
      end if
   end function
   real(dp) function rperks(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      real(dp)::u
      call random_number(u)
      v=qperks(u,scale,shape)
   end function

   elemental real(dp) function dlindley(x,theta,log_density) result(v)
      real(dp),intent(in)::x,theta
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(theta<=0.0_dp)then
      v=nanv()
      return
      end if
      if(x<0.0_dp)then
      ld=-infv()
      else
      ld=2.0_dp*log(theta)+log1p_v(x)-theta*x-log1p_v(theta)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function plindley(q,theta) result(v)
      real(dp),intent(in)::q,theta
      if(theta<=0.0_dp)then
      v=nanv()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=-expm1_v(-theta*q+log1p_v(q/(1.0_dp+1.0_dp/theta)))
      end if
   end function
   real(dp) function qlindley(p,theta) result(v)
      real(dp),intent(in)::p,theta
      real(dp)::lo,hi,mid
      integer::i
      if(theta<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=nanv()
      return
      end if
      if(p==1.0_dp)then
      v=infv()
      return
      end if
      lo=0.0_dp
      hi=max(1.0_dp,2.0_dp/theta)
      do while(plindley(hi,theta)<p)
      hi=2.0_dp*hi
      end do
      do i=1,100
         mid=(lo+hi)/2
         if(plindley(mid,theta)<p)then
      lo=mid
      else
      hi=mid
      end if
      end do
      v=(lo+hi)/2
   end function
   real(dp) function rlindley(theta) result(v)
      real(dp),intent(in)::theta
      real(dp)::u
      call random_number(u)
      v=qlindley(u,theta)
   end function

   elemental real(dp) function dnakagami(x,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,scale,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=nanv()
      return
      end if
      if(x<=0.0_dp)then
      ld=-infv()
      else
         ld=log(2.0_dp)+log(x)+(shape-1.0_dp)*log(x*x)- &
            x*x/(scale/shape)-log_gamma(shape)-shape*log(scale/shape)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pnakagami(q,scale,shape) result(v)
      real(dp),intent(in)::q,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=nanv()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
      v=regularized_gamma_p(shape,shape*q*q/scale)
      end if
   end function
   real(dp) function qnakagami(p,scale,shape) result(v)
      real(dp),intent(in)::p,scale,shape
      if(scale<=0.0_dp.or.shape<=0.0_dp)then
      v=nanv()
      else
      v=sqrt(scale*gamma_quantile(p,shape,1.0_dp)/shape)
      end if
   end function
   real(dp) function rnakagami(scale,shape) result(v)
      real(dp),intent(in)::scale,shape
      v=sqrt(random_gamma(shape,scale/shape))
   end function

   elemental real(dp) function dmaxwell(x,rate,log_density) result(v)
      real(dp),intent(in)::x,rate
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(rate<=0.0_dp)then
      v=nanv()
      return
      end if
      if(x<=0.0_dp)then
      ld=-infv()
      else
      ld=0.5_dp*log(2.0_dp/pi)+1.5_dp*log(rate)+ &
           2.0_dp*log(x)-0.5_dp*rate*x*x
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pmaxwell(q,rate) result(v)
      real(dp),intent(in)::q,rate
      if(rate<=0.0_dp)then
      v=nanv()
      else if(q<=0.0_dp)then
      v=0.0_dp
      else
         v=erf(q*sqrt(rate/2.0_dp))-q*exp(-0.5_dp*rate*q*q)* &
           sqrt(2.0_dp*rate/pi)
      end if
   end function
   real(dp) function qmaxwell(p,rate) result(v)
      real(dp),intent(in)::p,rate
      if(rate<=0.0_dp)then
      v=nanv()
      else
      v=sqrt(2.0_dp*gamma_quantile(p,1.5_dp,1.0_dp)/rate)
      end if
   end function
   real(dp) function rmaxwell(rate) result(v)
      real(dp),intent(in)::rate
      v=sqrt(2.0_dp*random_gamma(1.5_dp)/rate)
   end function

   elemental real(dp) function dbenini(x,y0,shape,log_density) result(v)
      real(dp),intent(in)::x,y0,shape
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::t,ld
      lg=.false.
      if(present(log_density))lg=log_density
      if(y0<=0.0_dp.or.shape<=0.0_dp)then
      v=nanv()
      return
      end if
      if(x<=y0)then
      ld=-infv()
      else
      t=log(x/y0)
      ld=log(2.0_dp*shape)-shape*t*t+log(t)-log(x)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function pbenini(q,y0,shape) result(v)
      real(dp),intent(in)::q,y0,shape
      if(y0<=0.0_dp.or.shape<=0.0_dp)then
      v=nanv()
      else if(q<=y0)then
      v=0.0_dp
      else
      v=-expm1_v(-shape*log(q/y0)**2)
      end if
   end function
   elemental real(dp) function qbenini(p,y0,shape) result(v)
      real(dp),intent(in)::p,y0,shape
      if(y0<=0.0_dp.or.shape<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
         v=nanv()
      else
      v=y0*exp(sqrt(-log1p_v(-p)/shape))
      end if
   end function
   real(dp) function rbenini(y0,shape) result(v)
      real(dp),intent(in)::y0,shape
      real(dp)::u
      call random_number(u)
      v=qbenini(u,y0,shape)
   end function

   elemental real(dp) function dlevy(x,location,scale,log_density) result(v)
      real(dp),intent(in)::x,location,scale
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld,z
      lg=.false.
      if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then
      v=nanv()
      return
      end if
      z=x-location
      if(z<=0.0_dp)then
      ld=-infv()
      else
      ld=0.5_dp*log(scale/(2.0_dp*pi))-1.5_dp*log(z)-0.5_dp*scale/z
      end if
      v=merge(ld,exp(ld),lg)
   end function
   elemental real(dp) function plevy(q,location,scale) result(v)
      real(dp),intent(in)::q,location,scale
      if(scale<=0.0_dp)then
      v=nanv()
      else if(q<=location)then
      v=0.0_dp
      else
      v=erfc(sqrt(0.5_dp*scale/(q-location)))
      end if
   end function
   elemental real(dp) function qlevy(p,location,scale) result(v)
      real(dp),intent(in)::p,location,scale
      real(dp)::z
      if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then
      v=nanv()
      else if(p==0.0_dp)then
      v=location
      else if(p==1.0_dp)then
      v=infv()
      else
      z=normal_quantile(0.5_dp*p)
      v=location+scale/(z*z)
      end if
   end function
   real(dp) function rlevy(location,scale) result(v)
      real(dp),intent(in)::location,scale
      real(dp)::z
      do
      z=random_gamma(0.5_dp,2.0_dp)
      if(z>0.0_dp)exit
      end do
      v=location+scale/z
   end function

   pure real(dp) function log_bessel_i_integer(n,x) result(v)
      integer,intent(in)::n
      real(dp),intent(in)::x
      real(dp)::logterm,mx,s
      integer::k
      if(x==0.0_dp)then
         if(n==0)then
      v=0.0_dp
      else
      v=-infv()
      end if
         return
      end if
      mx=-huge(1.0_dp)
      s=0.0_dp
      do k=0,10000
         logterm=(real(2*k+n,dp))*log(x/2.0_dp)- &
            log_gamma(real(k+1,dp))-log_gamma(real(k+n+1,dp))
         if(logterm>mx)then
            if(mx>-huge(1.0_dp)/2)s=s*exp(mx-logterm)
            mx=logterm
      s=s+1.0_dp
         else
            s=s+exp(logterm-mx)
         end if
         if(k>10.and.exp(logterm-mx)<1.0e-15_dp*s)exit
      end do
      v=mx+log(s)
   end function
   real(dp) function dskellam(x,mu1,mu2,log_density) result(v)
      integer,intent(in)::x
      real(dp),intent(in)::mu1,mu2
      logical,intent(in),optional::log_density
      logical::lg
      real(dp)::ld,arg
      lg=.false.
      if(present(log_density))lg=log_density
      if(mu1<0.0_dp.or.mu2<0.0_dp)then
      v=nanv()
      return
      end if
      if(mu1==0.0_dp.and.mu2==0.0_dp)then
         ld=merge(0.0_dp,-infv(),x==0)
      else if(mu1==0.0_dp)then
         if(x<=0)then
            ld=-mu2+real(-x,dp)*log(mu2)-log_gamma(real(-x+1,dp))
         else
      ld=-infv()
      end if
      else if(mu2==0.0_dp)then
         if(x>=0)then
            ld=-mu1+real(x,dp)*log(mu1)-log_gamma(real(x+1,dp))
         else
      ld=-infv()
      end if
      else
         arg=2.0_dp*sqrt(mu1*mu2)
         ld=-mu1-mu2+0.5_dp*real(x,dp)*log(mu1/mu2)+ &
            log_bessel_i_integer(abs(x),arg)
      end if
      v=merge(ld,exp(ld),lg)
   end function
   integer function rskellam(mu1,mu2) result(x)
      real(dp),intent(in)::mu1,mu2
      x=random_poisson(mu1)-random_poisson(mu2)
   end function
end module vgam_actuarial
