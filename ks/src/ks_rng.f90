! SPDX-License-Identifier: GPL-2.0-only
module ks_rng
   use ks_kinds, only: dp, i8, pi
   implicit none
   private
   public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_gamma, rng_chisq

   type :: rng_state
      integer(i8) :: state = 123456789_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state
contains
   subroutine rng_seed(rng, seed)
      type(rng_state), intent(inout) :: rng
      integer(i8), intent(in) :: seed
      integer(i8), parameter :: m=2147483647_i8
      rng%state=modulo(seed,m)
      if (rng%state==0_i8) rng%state=1_i8
      rng%has_spare=.false.
      rng%spare=0.0_dp
   end subroutine rng_seed

   function rng_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(i8), parameter :: a=16807_i8, m=2147483647_i8
      rng%state=mod(a*rng%state,m)
      if (rng%state<=0_i8) rng%state=rng%state+m
      u=real(rng%state,dp)/real(m,dp)
      u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u))
   end function rng_uniform

   function rng_normal(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z,u1,u2,r
      if (rng%has_spare) then
         z=rng%spare
         rng%has_spare=.false.
         return
      end if
      u1=rng_uniform(rng); u2=rng_uniform(rng)
      r=sqrt(-2.0_dp*log(u1))
      z=r*cos(2.0_dp*pi*u2)
      rng%spare=r*sin(2.0_dp*pi*u2)
      rng%has_spare=.true.
   end function rng_normal

   recursive function rng_gamma(rng,shape) result(x)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: shape
      real(dp) :: x,d,c,z,u,v
      if (shape<=0.0_dp) then
         x=0.0_dp; return
      end if
      if (shape<1.0_dp) then
         u=rng_uniform(rng)
         x=rng_gamma(rng,shape+1.0_dp)*u**(1.0_dp/shape)
         return
      end if
      d=shape-1.0_dp/3.0_dp
      c=1.0_dp/sqrt(9.0_dp*d)
      do
         z=rng_normal(rng)
         v=1.0_dp+c*z
         if (v<=0.0_dp) cycle
         v=v*v*v
         u=rng_uniform(rng)
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
      end do
      x=d*v
   end function rng_gamma

   function rng_chisq(rng,df) result(x)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: df
      real(dp) :: x
      x=2.0_dp*rng_gamma(rng,0.5_dp*df)
   end function rng_chisq
end module ks_rng
