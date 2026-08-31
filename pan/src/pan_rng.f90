! SPDX-License-Identifier: GPL-3.0-only
! Derived from computational code in R package pan 2.0.
! Upstream authorship/maintenance: Joseph L. Schafer and Jing Hua Zhao.
! Random-number kernels derived from the computational role of pan/src/pan.f.
module pan_rng
   use pan_kinds, only : dp, i8
   implicit none
   private

   real(dp), parameter :: pi_dp = acos(-1.0_dp)
   integer(i8), parameter :: pm_modulus = 2147483647_i8
   integer(i8), parameter :: pm_multiplier = 16807_i8

   type, public :: rng_state
      integer(i8) :: state = 1_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state

   public :: rng_seed
   public :: rng_uniform
   public :: rng_normal
   public :: rng_gamma

contains

   subroutine rng_seed(rng, seed)
      type(rng_state), intent(out) :: rng !! Random-number generator state initialized by this call.
      integer, intent(in) :: seed !! Positive deterministic seed; nonpositive values are mapped to one.

      integer(i8) :: s

      s = int(seed, i8)
      if (s <= 0_i8) s = 1_i8
      s = modulo(s, pm_modulus)
      if (s == 0_i8) s = 1_i8

      rng%state = s
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine rng_seed

   function rng_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng !! Mutable Park-Miller generator state.
      real(dp) :: u

      rng%state = modulo(pm_multiplier * rng%state, pm_modulus)
      if (rng%state <= 0_i8) rng%state = 1_i8
      u = real(rng%state, dp) / real(pm_modulus, dp)

      if (u <= 0.0_dp) u = epsilon(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
   end function rng_uniform

   function rng_normal(rng) result(z)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for a standard-normal draw.
      real(dp) :: z

      real(dp) :: radius
      real(dp) :: theta
      real(dp) :: u1
      real(dp) :: u2

      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
         return
      end if

      u1 = rng_uniform(rng)
      u2 = rng_uniform(rng)
      radius = sqrt(-2.0_dp * log(u1))
      theta = 2.0_dp * pi_dp * u2

      z = radius * cos(theta)
      rng%spare = radius * sin(theta)
      rng%has_spare = .true.
   end function rng_normal

   recursive function rng_gamma(rng, shape) result(x)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the gamma draw.
      real(dp), intent(in) :: shape !! Positive gamma shape parameter; the scale is one.
      real(dp) :: x

      real(dp) :: c
      real(dp) :: d
      real(dp) :: u
      real(dp) :: v
      real(dp) :: z

      if (shape <= 0.0_dp) then
         x = 0.0_dp
         return
      end if

      if (shape < 1.0_dp) then
         u = rng_uniform(rng)
         x = rng_gamma(rng, shape + 1.0_dp) * u**(1.0_dp / shape)
         return
      end if

      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)

      do
         z = rng_normal(rng)
         v = 1.0_dp + c * z
         if (v <= 0.0_dp) cycle
         v = v * v * v
         u = rng_uniform(rng)

         if (u < 1.0_dp - 0.0331_dp * z**4) exit
         if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
      end do

      x = d * v
   end function rng_gamma

end module pan_rng
