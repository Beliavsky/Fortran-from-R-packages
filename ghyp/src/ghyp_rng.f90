! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_rng
   use ghyp_kinds, only : dp, i8, pi
   implicit none
   private
   public :: rng_state, seed_rng, uniform_rng, normal_rng, gamma_rng

   type :: rng_state
      integer(i8) :: state = 88172645463325252_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state

contains

   subroutine seed_rng(rng, seed)
      type(rng_state), intent(inout) :: rng
      integer(i8), intent(in) :: seed
      rng%state = seed
      if (rng%state == 0_i8) rng%state = 88172645463325252_i8
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine seed_rng

   function next_u64(rng) result(x)
      type(rng_state), intent(inout) :: rng
      integer(i8) :: x
      x = rng%state
      x = ieor(x,shiftl(x,13))
      x = ieor(x,shiftr(x,7))
      x = ieor(x,shiftl(x,17))
      rng%state = x
   end function next_u64

   function uniform_rng(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(i8) :: x
      x = next_u64(rng)
      u = real(iand(x,int(z'001FFFFFFFFFFFFF',i8)),dp)/9007199254740992.0_dp
      u = max(u,epsilon(1.0_dp))
      u = min(u,1.0_dp-epsilon(1.0_dp))
   end function uniform_rng

   function normal_rng(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z, u1, u2, r
      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
         return
      end if
      u1 = uniform_rng(rng)
      u2 = uniform_rng(rng)
      r = sqrt(-2.0_dp*log(u1))
      z = r*cos(2.0_dp*pi*u2)
      rng%spare = r*sin(2.0_dp*pi*u2)
      rng%has_spare = .true.
   end function normal_rng

   recursive function gamma_rng(shape, scale, rng) result(x)
      real(dp), intent(in) :: shape, scale
      type(rng_state), intent(inout) :: rng
      real(dp) :: x, d, c, z, v, u
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = 0.0_dp
         return
      end if
      if (shape < 1.0_dp) then
         x = gamma_rng(shape+1.0_dp,scale,rng)*uniform_rng(rng)**(1.0_dp/shape)
         return
      end if
      d = shape-1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z = normal_rng(rng)
            v = 1.0_dp+c*z
            if (v > 0.0_dp) exit
         end do
         v = v*v*v
         u = uniform_rng(rng)
         if (u < 1.0_dp-0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
      end do
      x = scale*d*v
   end function gamma_rng

end module ghyp_rng
