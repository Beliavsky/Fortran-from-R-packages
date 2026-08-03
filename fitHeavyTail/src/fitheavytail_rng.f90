! SPDX-License-Identifier: GPL-3.0-only
module fitheavytail_rng
   use fitheavytail_kinds, only: dp
   implicit none
   private
   public :: rng_state, seed_rng, uniform_random, normal_random, gamma_random, random_mvt_identity

   type :: rng_state
      integer(kind=8) :: state = 104729_8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state

contains

   subroutine seed_rng(rng, seed)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: seed
      rng%state = max(1_8, int(abs(seed),8))
      rng%has_spare = .false.
   end subroutine seed_rng

   function uniform_random(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(kind=8), parameter :: a = 48271_8, m = 2147483647_8
      rng%state = modulo(a*rng%state, m)
      u = real(rng%state,dp)/real(m,dp)
      u = max(tiny(1.0_dp), min(1.0_dp-epsilon(1.0_dp), u))
   end function uniform_random

   function normal_random(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z, r, theta
      real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp
      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
      else
         r = sqrt(-2.0_dp*log(uniform_random(rng)))
         theta = two_pi*uniform_random(rng)
         z = r*cos(theta)
         rng%spare = r*sin(theta)
         rng%has_spare = .true.
      end if
   end function normal_random

   recursive function gamma_random(rng, shape) result(g)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: shape
      real(dp) :: g, d, c, x, v, u
      if (shape <= 0.0_dp) then
         g = 0.0_dp
      else if (shape < 1.0_dp) then
         g = gamma_random(rng,shape+1.0_dp)*uniform_random(rng)**(1.0_dp/shape)
      else
         d = shape - 1.0_dp/3.0_dp
         c = 1.0_dp/sqrt(9.0_dp*d)
         do
            x = normal_random(rng)
            v = 1.0_dp + c*x
            if (v <= 0.0_dp) cycle
            v = v*v*v
            u = uniform_random(rng)
            if (u < 1.0_dp-0.0331_dp*x**4) exit
            if (log(u) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
         end do
         g = d*v
      end if
   end function gamma_random

   subroutine random_mvt_identity(nobs, nvar, nu, x, seed)
      integer, intent(in) :: nobs, nvar
      real(dp), intent(in) :: nu
      real(dp), intent(out) :: x(nobs,nvar)
      integer, intent(in), optional :: seed
      type(rng_state) :: rng
      real(dp) :: chi
      integer :: i, j
      call seed_rng(rng, 12345)
      if (present(seed)) call seed_rng(rng,seed)
      do i = 1, nobs
         chi = 2.0_dp*gamma_random(rng,0.5_dp*nu)
         do j = 1, nvar
            x(i,j) = normal_random(rng)*sqrt(nu/max(chi,tiny(1.0_dp)))
         end do
      end do
   end subroutine random_mvt_identity

end module fitheavytail_rng
