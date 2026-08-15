! SPDX-License-Identifier: GPL-3.0-only
module gamlss_random
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use gamlss_kinds, only : dp, pi
   implicit none
   private
   public :: random_normal, random_exponential, random_gamma, random_beta
   public :: random_poisson, random_binomial
contains
   real(dp) function random_normal() result(z)
      real(dp) :: u1,u2
      call random_number(u1)
       call random_number(u2)
      u1=max(u1,tiny(1.0_dp))
      z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   real(dp) function random_exponential(rate) result(x)
      real(dp),intent(in),optional::rate
      real(dp)::u,r
      r=1.0_dp
      if(present(rate))r=rate
      call random_number(u)
       x=-log(max(u,tiny(1.0_dp)))/r
   end function random_exponential

   recursive real(dp) function random_gamma(shape,scale) result(x)
      real(dp),intent(in)::shape
      real(dp),intent(in),optional::scale
      real(dp)::sc,d,c,z,u,v
      sc=1.0_dp
      if(present(scale))sc=scale
      if(shape<=0.0_dp .or. sc<=0.0_dp)then
         x=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if(shape<1.0_dp)then
         call random_number(u)
         x=random_gamma(shape+1.0_dp,sc)*u**(1.0_dp/shape)
         return
      end if
      d=shape-1.0_dp/3.0_dp
      c=1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z=random_normal()
      v=(1.0_dp+c*z)**3
            if(v>0.0_dp)exit
         end do
         call random_number(u)
         if(u<1.0_dp-0.0331_dp*z**4)exit
         if(log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v)))exit
      end do
      x=sc*d*v
   end function random_gamma

   real(dp) function random_beta(a,b) result(x)
      real(dp),intent(in)::a,b
      real(dp)::g1,g2
      g1=random_gamma(a)
      g2=random_gamma(b)
      x=g1/(g1+g2)
   end function random_beta

   integer function random_poisson(lambda) result(k)
      real(dp),intent(in)::lambda
      real(dp)::l,p,u,z
      integer::n
      if(lambda<0.0_dp)then
      k=-1
      return
      else if(lambda==0.0_dp)then
      k=0
      return
      else if(lambda<30.0_dp)then
         l=exp(-lambda)
      p=1.0_dp
      n=0
         do
            n=n+1
      call random_number(u)
      p=p*u
            if(p<=l)exit
         end do
         k=n-1
      else
         do
            z=random_normal()
      n=nint(lambda+sqrt(lambda)*z)
            if(n>=0)exit
         end do
         k=n
      end if
   end function random_poisson

   integer function random_binomial(n,p) result(k)
      integer,intent(in)::n
      real(dp),intent(in)::p
      integer::i
      real(dp)::u
      if(n<0 .or. p<0.0_dp .or. p>1.0_dp)then
      k=-1
      return
      end if
      k=0
      if(n<100)then
         do i=1,n
      call random_number(u)
      if(u<p)k=k+1
      end do
      else
         do
            k=nint(real(n,dp)*p+sqrt(real(n,dp)*p*(1.0_dp-p))*random_normal())
            if(k>=0 .and. k<=n)exit
         end do
      end if
   end function random_binomial
end module gamlss_random
