! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_extremes
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use vgam_kinds, only : dp
   use vgam_distributions, only : dgumbel, pgumbel, qgumbel, rgumbel
   implicit none
   private
   public :: dgev, pgev, qgev, rgev
   public :: dgpd, pgpd, qgpd, rgpd
   public :: dlaplace, plaplace, qlaplace, rlaplace
   public :: dkumar, pkumar, qkumar, rkumar
contains
   elemental real(dp) function nanv() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nanv

   elemental real(dp) function dgev(x,location,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,location,scale,shape
      logical,intent(in),optional::log_density
      real(dp)::z,t,ld
      logical::lg
      lg=.false.; if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then; v=nanv(); return; end if
      z=(x-location)/scale
      if(abs(shape)<1.0e-10_dp)then
         v=dgumbel(x,location,scale,lg); return
      end if
      t=1.0_dp+shape*z
      if(t<=0.0_dp)then
         ld=-huge(1.0_dp)
      else
         ld=-log(scale)+(-1.0_dp/shape-1.0_dp)*log(t)-t**(-1.0_dp/shape)
      end if
      v=merge(ld,exp(ld),lg)
   end function dgev

   elemental real(dp) function pgev(q,location,scale,shape) result(v)
      real(dp),intent(in)::q,location,scale,shape
      real(dp)::z,t
      if(scale<=0.0_dp)then; v=nanv(); return; end if
      if(abs(shape)<1.0e-10_dp)then; v=pgumbel(q,location,scale); return; end if
      z=(q-location)/scale; t=1.0_dp+shape*z
      if(t<=0.0_dp)then
         v=merge(0.0_dp,1.0_dp,shape>0.0_dp)
      else
         v=exp(-t**(-1.0_dp/shape))
      end if
   end function pgev

   elemental real(dp) function qgev(p,location,scale,shape) result(v)
      real(dp),intent(in)::p,location,scale,shape
      if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then; v=nanv(); return; end if
      if(abs(shape)<1.0e-10_dp)then
         v=qgumbel(p,location,scale)
      else
         v=location+scale*((-log(p))**(-shape)-1.0_dp)/shape
      end if
   end function qgev

   real(dp) function rgev(location,scale,shape) result(v)
      real(dp),intent(in)::location,scale,shape
      real(dp)::u
      call random_number(u); v=qgev(u,location,scale,shape)
   end function rgev

   elemental real(dp) function dgpd(x,threshold,scale,shape,log_density) result(v)
      real(dp),intent(in)::x,threshold,scale,shape
      logical,intent(in),optional::log_density
      real(dp)::z,t,ld
      logical::lg
      lg=.false.; if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then; v=nanv(); return; end if
      z=(x-threshold)/scale
      if(z<0.0_dp)then
         ld=-huge(1.0_dp)
      else if(abs(shape)<1.0e-10_dp)then
         ld=-log(scale)-z
      else
         t=1.0_dp+shape*z
         if(t<=0.0_dp)then
            ld=-huge(1.0_dp)
         else
            ld=-log(scale)+(-1.0_dp/shape-1.0_dp)*log(t)
         end if
      end if
      v=merge(ld,exp(ld),lg)
   end function dgpd

   elemental real(dp) function pgpd(q,threshold,scale,shape) result(v)
      real(dp),intent(in)::q,threshold,scale,shape
      real(dp)::z,t
      if(scale<=0.0_dp)then;v=nanv();return;end if
      z=(q-threshold)/scale
      if(z<=0.0_dp)then;v=0.0_dp
      else if(abs(shape)<1.0e-10_dp)then;v=1.0_dp-exp(-z)
      else
         t=1.0_dp+shape*z
         if(t<=0.0_dp)then;v=1.0_dp
         else;v=1.0_dp-t**(-1.0_dp/shape)
         end if
      end if
   end function pgpd

   elemental real(dp) function qgpd(p,threshold,scale,shape) result(v)
      real(dp),intent(in)::p,threshold,scale,shape
      if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;v=nanv();return;end if
      if(abs(shape)<1.0e-10_dp)then
         v=threshold-scale*log(1.0_dp-p)
      else
         v=threshold+scale*((1.0_dp-p)**(-shape)-1.0_dp)/shape
      end if
   end function qgpd

   real(dp) function rgpd(threshold,scale,shape) result(v)
      real(dp),intent(in)::threshold,scale,shape
      real(dp)::u
      call random_number(u);v=qgpd(u,threshold,scale,shape)
   end function rgpd

   elemental real(dp) function dlaplace(x,location,scale,log_density) result(v)
      real(dp),intent(in)::x,location,scale
      logical,intent(in),optional::log_density
      real(dp)::ld;logical::lg
      lg=.false.;if(present(log_density))lg=log_density
      if(scale<=0.0_dp)then;v=nanv();return;end if
      ld=-log(2.0_dp*scale)-abs(x-location)/scale
      v=merge(ld,exp(ld),lg)
   end function dlaplace
   elemental real(dp) function plaplace(q,location,scale) result(v)
      real(dp),intent(in)::q,location,scale
      if(scale<=0.0_dp)then;v=nanv()
      else if(q<location)then;v=0.5_dp*exp((q-location)/scale)
      else;v=1.0_dp-0.5_dp*exp(-(q-location)/scale)
      end if
   end function plaplace
   elemental real(dp) function qlaplace(p,location,scale) result(v)
      real(dp),intent(in)::p,location,scale
      if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;v=nanv()
      else if(p<0.5_dp)then;v=location+scale*log(2.0_dp*p)
      else;v=location-scale*log(2.0_dp*(1.0_dp-p))
      end if
   end function qlaplace
   real(dp) function rlaplace(location,scale) result(v)
      real(dp),intent(in)::location,scale
      real(dp)::u
      call random_number(u);v=qlaplace(u,location,scale)
   end function rlaplace

   elemental real(dp) function dkumar(x,a,b,log_density) result(v)
      real(dp),intent(in)::x,a,b
      logical,intent(in),optional::log_density
      real(dp)::ld;logical::lg
      lg=.false.;if(present(log_density))lg=log_density
      if(a<=0.0_dp.or.b<=0.0_dp)then;v=nanv();return;end if
      if(x<=0.0_dp.or.x>=1.0_dp)then;ld=-huge(1.0_dp)
      else;ld=log(a)+log(b)+(a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x**a)
      end if
      v=merge(ld,exp(ld),lg)
   end function dkumar
   elemental real(dp) function pkumar(q,a,b) result(v)
      real(dp),intent(in)::q,a,b
      if(a<=0.0_dp.or.b<=0.0_dp)then;v=nanv()
      else if(q<=0.0_dp)then;v=0.0_dp
      else if(q>=1.0_dp)then;v=1.0_dp
      else;v=1.0_dp-(1.0_dp-q**a)**b
      end if
   end function pkumar
   elemental real(dp) function qkumar(p,a,b) result(v)
      real(dp),intent(in)::p,a,b
      if(a<=0.0_dp.or.b<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;v=nanv()
      else;v=(1.0_dp-(1.0_dp-p)**(1.0_dp/b))**(1.0_dp/a)
      end if
   end function qkumar
   real(dp) function rkumar(a,b) result(v)
      real(dp),intent(in)::a,b
      real(dp)::u
      call random_number(u);v=qkumar(u,a,b)
   end function rkumar
end module vgam_extremes
