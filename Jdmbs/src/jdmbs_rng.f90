! SPDX-License-Identifier: GPL-2.0-or-later
module jdmbs_rng
   use jdmbs_kinds, only : dp, int64
   implicit none
   private

   integer(int64), parameter :: pm_modulus = 2147483647_int64
   integer(int64), parameter :: pm_multiplier = 16807_int64
   integer(int64), parameter :: pm_q = 127773_int64
   integer(int64), parameter :: pm_r = 2836_int64
   real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp

   type, public :: rng_state
      integer(int64) :: state = 123456789_int64
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state

   public :: rng_seed, rng_uniform, rng_normal, rng_integer

contains

   subroutine rng_seed(rng, seed)
      type(rng_state), intent(out) :: rng
      integer(int64), intent(in) :: seed
      integer(int64) :: s

      s = modulo(abs(seed), pm_modulus - 1_int64) + 1_int64
      rng%state = s
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine rng_seed

   function rng_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(int64) :: hi, lo, test

      hi = rng%state / pm_q
      lo = modulo(rng%state, pm_q)
      test = pm_multiplier * lo - pm_r * hi
      if (test > 0_int64) then
         rng%state = test
      else
         rng%state = test + pm_modulus
      end if
      u = real(rng%state, dp) / real(pm_modulus, dp)
   end function rng_uniform

   function rng_normal(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z
      real(dp) :: radius, angle, u1, u2

      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
         return
      end if
      u1 = max(rng_uniform(rng), tiny(1.0_dp))
      u2 = rng_uniform(rng)
      radius = sqrt(-2.0_dp * log(u1))
      angle = two_pi * u2
      z = radius * cos(angle)
      rng%spare = radius * sin(angle)
      rng%has_spare = .true.
   end function rng_normal

   function rng_integer(rng, lower, upper) result(value)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: lower, upper
      integer :: value

      value = lower + int(rng_uniform(rng) * real(upper - lower + 1, dp))
      value = min(value, upper)
   end function rng_integer

end module jdmbs_rng
