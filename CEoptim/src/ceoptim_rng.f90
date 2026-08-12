! SPDX-License-Identifier: GPL-2.0-or-later
module ceoptim_rng
   use ceoptim_kinds, only : dp, i64
   implicit none
   private

   integer(i64), parameter :: pm_mod = 2147483647_i64
   integer(i64), parameter :: pm_mul = 16807_i64
   real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

   type, public :: rng_state
      integer(i64) :: state = 123456789_i64
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state

   public :: rng_seed, rng_uniform, rng_normal, rng_gamma, rng_categorical

contains

   subroutine rng_seed(rng, seed)
      type(rng_state), intent(inout) :: rng
      integer(i64), intent(in) :: seed
      integer(i64) :: s

      s = modulo(abs(seed), pm_mod - 1_i64) + 1_i64
      rng%state = s
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine rng_seed

   function rng_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u

      rng%state = modulo(pm_mul * rng%state, pm_mod)
      if (rng%state <= 0_i64) rng%state = 1_i64
      u = real(rng%state, dp) / real(pm_mod, dp)
   end function rng_uniform

   function rng_normal(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z
      real(dp) :: u1, u2, r

      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
         return
      end if

      u1 = max(rng_uniform(rng), tiny(1.0_dp))
      u2 = rng_uniform(rng)
      r = sqrt(-2.0_dp * log(u1))
      z = r * cos(two_pi * u2)
      rng%spare = r * sin(two_pi * u2)
      rng%has_spare = .true.
   end function rng_normal

   recursive function rng_gamma(rng, shape) result(x)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: shape
      real(dp) :: x
      real(dp) :: d, c, z, v, u

      if (shape <= 0.0_dp) then
         x = 0.0_dp
         return
      end if

      if (shape < 1.0_dp) then
         u = max(rng_uniform(rng), tiny(1.0_dp))
         x = rng_gamma(rng, shape + 1.0_dp) * u**(1.0_dp / shape)
         return
      end if

      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         do
            z = rng_normal(rng)
            v = 1.0_dp + c * z
            if (v > 0.0_dp) exit
         end do
         v = v**3
         u = rng_uniform(rng)
         if (u < 1.0_dp - 0.0331_dp * z**4) exit
         if (log(u) < 0.5_dp * z*z + d * (1.0_dp - v + log(v))) exit
      end do
      x = d * v
   end function rng_gamma

   integer function rng_categorical(rng, probs, ncat) result(k)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: probs(:)
      integer, intent(in) :: ncat
      real(dp) :: u, csum, total
      integer :: j

      total = sum(probs(1:ncat))
      if (total <= 0.0_dp) then
         k = 0
         return
      end if
      u = rng_uniform(rng) * total
      csum = 0.0_dp
      do j = 1, ncat
         csum = csum + probs(j)
         if (u <= csum) then
            k = j - 1
            return
         end if
      end do
      k = ncat - 1
   end function rng_categorical

end module ceoptim_rng
