! Wiener-germ approximations for the noncentral chi-square distribution.
! Translated from DPQ R code and the package's original Fortran noncechi routine.
! SPDX-License-Identifier: GPL-2.0-or-later
module dpq_wiener
   use r_compat, only: dp
   use dpq_core, only: prob_output, log1p_dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: h0, h1, hnt, h2, h, gnt, g2, sw, qs_w
   public :: z0_w, z_f, z_s, pchisq_w

contains

   pure elemental real(dp) function h0(y) result(v)
      real(dp), intent(in) :: y
      real(dp) :: y2
      y2=y*y
      if(y==0.0_dp) then
         v=0.0_dp
      else
         v=((1.0_dp-y)*log(1.0_dp-y)+y-0.5_dp*y2)/y2
      end if
   end function h0

   pure elemental real(dp) function h1(y) result(v)
      real(dp), intent(in) :: y
      real(dp) :: y2
      y2=y*y
      if(y==0.0_dp) then
         v=0.0_dp
      else
         v=((1.0_dp-y)*log1p_dp(-y)+y*(1.0_dp-0.5_dp*y))/y2
      end if
   end function h1

   pure real(dp) function hnt(y,nterms) result(v)
      real(dp), intent(in) :: y
      integer, intent(in) :: nterms
      real(dp) :: s
      integer :: k
      if(nterms<=0) then
         if(y==0.0_dp) then
            v=0.0_dp
         else
            v=((1.0_dp-y)*log1p_dp(-y)+y*(1.0_dp-y/2.0_dp))/(y*y)
         end if
      else if(nterms<=1) then
         v=y/6.0_dp*(1.0_dp+y/2.0_dp)
      else
         s=0.0_dp
         do k=nterms,1,-1
            s=s*y+1.0_dp/(real(k+1,dp)*real(k+2,dp))
         end do
         v=y*s
      end if
   end function hnt

   pure elemental real(dp) function h2(y) result(v)
      real(dp), intent(in) :: y
      real(dp) :: ay,eps,c1,c2
      eps=epsilon(1.0_dp)
      ay=abs(y)
      c1=sqrt((10.0_dp/3.0_dp)*eps)
      c2=(5.0_dp*eps)**(1.0_dp/3.0_dp)
      if(ay<c1) then
         v=y/6.0_dp*(1.0_dp+y/2.0_dp)
      else if(ay<c2) then
         v=y*(1.0_dp/6.0_dp+y*(1.0_dp/12.0_dp+y/20.0_dp))
      else if(y==1.0_dp) then
         v=0.5_dp
      else
         v=((1.0_dp-y)*log1p_dp(-y)+y*(1.0_dp-y/2.0_dp))/(y*y)
      end if
   end function h2

   pure elemental real(dp) function h(y) result(v)
      real(dp), intent(in) :: y
      v=h2(y)
   end function h

   pure real(dp) function gnt(u,nterms,times_u) result(v)
      real(dp), intent(in) :: u
      integer, intent(in) :: nterms
      logical, intent(in), optional :: times_u
      logical :: tu
      real(dp) :: s
      integer :: k
      tu=.false.
      if(present(times_u))tu=times_u
      if(nterms<=0) then
         if(u==0.0_dp) then
            s=0.5_dp
            v=merge(0.0_dp,s,tu)
         else
            s=6.0_dp*h(u)/u-1.0_dp
            if(tu) then
            v=s
            else
            v=s/u
            end if
         end if
      else if(nterms<=1) then
         s=(1.0_dp+0.6_dp*u)/2.0_dp
         if(tu) then
         v=u*s
         else
         v=s
         end if
      else
         s=0.0_dp
         do k=nterms,0,-1
            s=s*u+1.0_dp/(real(k+3,dp)*real(k+4,dp))
         end do
         if(tu) then
         v=6.0_dp*u*s
         else
         v=6.0_dp*s
         end if
      end if
   end function gnt

   pure elemental real(dp) function g2(u) result(v)
      real(dp),intent(in)::u
      real(dp)::au,eps
      eps=epsilon(1.0_dp)
      au=abs(u)
      if(au<sqrt(2.5_dp*eps)) then
         v=(1.0_dp+0.6_dp*u)/2.0_dp
      else if(au<(3.5_dp*eps)**(1.0_dp/3.0_dp)) then
         v=(1.0_dp+u/5.0_dp*(3.0_dp+2.0_dp*u))/2.0_dp
      else if(au<((14.0_dp/3.0_dp)*eps)**0.25_dp) then
         v=(1.0_dp+u*(0.6_dp+u*(0.4_dp+u*2.0_dp/7.0_dp)))/2.0_dp
      else
         v=(6.0_dp*h(u)/u-1.0_dp)/u
      end if
   end function g2

   pure subroutine sw(x,df,ncp,ff,s,ifault)
      real(dp),intent(in)::x,df,ncp
      real(dp),intent(out)::ff,s
      integer,intent(out),optional::ifault
      real(dp)::mu2,e
      if(present(ifault))ifault=0
      if(ncp<0.0_dp) then
         ff=ieee_value(0.0_dp,ieee_quiet_nan)
         s=ff
         if(present(ifault))ifault=1
         return
      end if
      if(df<=0.0_dp) then
         ff=ieee_value(0.0_dp,ieee_quiet_nan)
         s=ff
         if(present(ifault))ifault=2
         return
      end if
      mu2=ncp/df
      if(mu2==0.0_dp) then
         ff=1.0_dp
         s=x/df
         return
      end if
      ff=sqrt(1.0_dp+4.0_dp*x*mu2/df)
      e=2.0_dp*mu2*x/df
      if(e<(8.0_dp/7.0_dp)*epsilon(1.0_dp)**0.25_dp) then
         s=x/df*(1.0_dp+e*(-0.5_dp+e*(0.5_dp-5.0_dp*e/8.0_dp)))
      else
         s=(ff-1.0_dp)/(2.0_dp*mu2)
      end if
   end subroutine sw

   pure real(dp) function qs_w(x,df,ncp,eps1,smax) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp),intent(in),optional::eps1,smax
      real(dp)::ff,s,mu2,e1,sm,is
      e1=0.5_dp
      if(present(eps1))e1=eps1
      sm=1.0e100_dp
      if(present(smax))sm=smax
      call sw(x,df,ncp,ff,s)
      mu2=ncp/df
      if(s<e1) then
         is=1.0_dp/(1.0_dp-s)
         v=2.0_dp*(mu2-is*(log(s)*is+1.0_dp))
      else if(s>sm) then
         v=2.0_dp*mu2+(1.0_dp+2.0_dp*h(1.0_dp-1.0_dp/s))/s
      else
         v=(ff+2.0_dp*h(1.0_dp-1.0_dp/s))/s
      end if
   end function qs_w

   pure real(dp) function z0_w(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::ff,s,s1,qv
      call sw(x,df,ncp,ff,s)
      s1=s-1.0_dp
      qv=qs_w(x,df,ncp)
      v=0.5_dp*df*s1*s1*qv-log(qv/ff)
   end function z0_w

   pure real(dp) function z_f(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::mu2,d2
      mu2=ncp/df
      d2=1.0_dp+3.0_dp*mu2
      v=z0_w(x,df,ncp)+(2.0_dp/9.0_dp)*d2*d2/(df*(1.0_dp+2.0_dp*mu2)**3)
   end function z_f

   pure real(dp) function z_s(x,df,ncp) result(v)
      real(dp),intent(in)::x,df,ncp
      real(dp)::ff,s,s1,qv,z0,m2s,f2,d2,dh,h1s,eta,gg,t1,t2
      call sw(x,df,ncp,ff,s)
      s1=s-1.0_dp
      qv=qs_w(x,df,ncp)
      z0=0.5_dp*df*s1*s1*qv-log(qv/ff)
      m2s=(ncp/df)*s
      f2=ff*ff
      d2=1.0_dp+3.0_dp*m2s
      dh=-1.5_dp*(1.0_dp+4.0_dp*m2s)/f2+(5.0_dp/3.0_dp)*d2*d2/(ff*f2)
      h1s=-h(1.0_dp-1.0_dp/s)
      eta=1.0_dp-ff/qv
      gg=eta/(s1*s1*ff)
      t2=gg*(3.0_dp-(0.5_dp+h(eta))*eta)
      t1=2.0_dp*d2/(s1*f2)
      v=z0+2.0_dp*(dh+t1+t2)/df
      if(.not.(v==v)) then
         ! Equivalent algebra from the original noncechi Fortran is often
         ! more stable away from s=1.
         eta=(ff-2.0_dp*h1s-s*ff)/(ff-2.0_dp*h1s)
         gg=eta/(s1*s1*ff)
         t2=gg*(3.0_dp-(0.5_dp+h(eta))*eta)
         v=z0+2.0_dp*(dh+t1+t2)/df
      end if
   end function z_s

   real(dp) function pchisq_w(x,df,ncp,variant,lower_tail,log_p,ifault) result(v)
      real(dp),intent(in)::x,df,ncp
      character(len=*),intent(in),optional::variant
      logical,intent(in),optional::lower_tail,log_p
      integer,intent(out),optional::ifault
      character(len=1)::vv
      real(dp)::ff,s,z,p
      logical::lt,lp
      lt=.true.
      lp=.false.
      vv='s'
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      if(present(variant))vv=variant(1:1)
      if(present(ifault))ifault=0
      if(ncp<0.0_dp)then
      if(present(ifault))ifault=1
      v=ieee_value(0.0_dp,ieee_quiet_nan)
      return
      end if
      if(df<=0.0_dp)then
      if(present(ifault))ifault=2
      v=ieee_value(0.0_dp,ieee_quiet_nan)
      return
      end if
      if(x<=0.0_dp)then
      v=prob_output(0.0_dp,lt,lp)
      return
      end if
      call sw(x,df,ncp,ff,s)
      if(abs(s-1.0_dp)<=4.0_dp*epsilon(1.0_dp)) then
         p=0.5_dp
      else
         if(vv=='f' .or. vv=='F')then
         z=z_f(x,df,ncp)
         else
         z=z_s(x,df,ncp)
         end if
         z=sqrt(abs(z))
         if(s<1.0_dp)z=-z
         p=0.5_dp*erfc(-z/sqrt(2.0_dp))
      end if
      v=prob_output(p,lt,lp)
   end function pchisq_w

end module dpq_wiener
