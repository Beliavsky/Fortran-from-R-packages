! Deterministic RNG and basic random deviates used by jomo MCMC kernels.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! The upstream C code also contains LGPL random-deviate routines derived from
! Barry Brown, James Lovato, and John Burkardt. This module is an independent
! modern Fortran implementation of equivalent distributional operations.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_rng
   use jomo_kinds, only : dp, i8
   implicit none
   private

   integer(i8), parameter :: pm_modulus = 2147483647_i8
   integer(i8), parameter :: pm_multiplier = 48271_i8
   real(dp), parameter :: pi_dp = acos(-1.0_dp)

   type, public :: rng_state
      integer(i8) :: state = 123456789_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state

   public :: rng_seed
   public :: rng_uniform
   public :: rng_normal
   public :: rng_gamma
   public :: rng_chisq
   public :: rng_student_t

contains

   subroutine rng_seed(rng, seed)
      type(rng_state), intent(inout) :: rng !! Mutable pseudorandom-number generator state to initialize.
      integer(i8), intent(in) :: seed !! Deterministic integer seed; zero is remapped to a nonzero valid state.

      rng%state = modulo(abs(seed), pm_modulus - 1_i8) + 1_i8
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine rng_seed

   function rng_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by one uniform draw.
      real(dp) :: u

      rng%state = modulo(pm_multiplier * rng%state, pm_modulus)
      if (rng%state <= 0_i8) rng%state = rng%state + pm_modulus - 1_i8
      u = real(rng%state, dp) / real(pm_modulus, dp)
      u = max(tiny(1.0_dp), min(1.0_dp - epsilon(1.0_dp), u))
   end function rng_uniform

   function rng_normal(rng, mean, sd) result(x)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the Gaussian draw.
      real(dp), intent(in) :: mean !! Mean of the requested normal distribution.
      real(dp), intent(in) :: sd !! Nonnegative standard deviation of the requested normal distribution.
      real(dp) :: x
      real(dp) :: r
      real(dp) :: theta

      if (sd < 0.0_dp) error stop "rng_normal: sd must be nonnegative"
      if (sd <= 0.0_dp) then
         x = mean
         return
      end if

      if (rng%has_spare) then
         x = mean + sd * rng%spare
         rng%has_spare = .false.
         return
      end if

      r = sqrt(-2.0_dp * log(rng_uniform(rng)))
      theta = 2.0_dp * pi_dp * rng_uniform(rng)
      x = mean + sd * r * cos(theta)
      rng%spare = r * sin(theta)
      rng%has_spare = .true.
   end function rng_normal

   recursive function rng_gamma(rng, shape, scale) result(x)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the gamma draw.
      real(dp), intent(in) :: shape !! Positive gamma shape parameter.
      real(dp), intent(in) :: scale !! Positive gamma scale parameter multiplying a unit-scale draw.
      real(dp) :: x
      real(dp) :: d
      real(dp) :: c
      real(dp) :: z
      real(dp) :: v
      real(dp) :: u

      if (shape <= 0.0_dp) error stop "rng_gamma: shape must be positive"
      if (scale <= 0.0_dp) error stop "rng_gamma: scale must be positive"

      if (shape < 1.0_dp) then
         u = rng_uniform(rng)
         x = rng_gamma(rng, shape + 1.0_dp, scale) * u ** (1.0_dp / shape)
         return
      end if

      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         z = rng_normal(rng, 0.0_dp, 1.0_dp)
         v = 1.0_dp + c * z
         if (v <= 0.0_dp) cycle
         v = v * v * v
         u = rng_uniform(rng)
         if (u < 1.0_dp - 0.0331_dp * z ** 4) exit
         if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
      end do
      x = scale * d * v
   end function rng_gamma

   function rng_chisq(rng, df) result(x)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the chi-square draw.
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
      real(dp) :: x

      if (df <= 0.0_dp) error stop "rng_chisq: df must be positive"
      x = rng_gamma(rng, 0.5_dp * df, 2.0_dp)
   end function rng_chisq

   function rng_student_t(rng, df) result(x)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the Student-t draw.
      real(dp), intent(in) :: df !! Positive Student-t degrees of freedom.
      real(dp) :: x
      real(dp) :: z
      real(dp) :: v

      if (df <= 0.0_dp) error stop "rng_student_t: df must be positive"
      z = rng_normal(rng, 0.0_dp, 1.0_dp)
      v = rng_chisq(rng, df)
      x = z / sqrt(v / df)
   end function rng_student_t

end module jomo_rng
