! SPDX-License-Identifier: GPL-2.0-or-later
module gnorm_rng
   use gnorm_kinds, only : dp, i8
   implicit none
   private

   integer(i8), parameter :: modulus = 2147483647_i8
   integer(i8), parameter :: multiplier = 16807_i8

   type, public :: gnorm_rng_state
      integer(i8) :: state = 104729_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: normal => rng_normal
      procedure :: gamma => rng_gamma
   end type gnorm_rng_state

contains

   subroutine rng_seed(self, seed_value)
      class(gnorm_rng_state), intent(inout) :: self
      integer(i8), intent(in) :: seed_value

      self%state = modulo(abs(seed_value), modulus - 1_i8) + 1_i8
      self%has_spare = .false.
      self%spare = 0.0_dp
   end subroutine rng_seed

   function rng_uniform(self) result(value)
      class(gnorm_rng_state), intent(inout) :: self
      real(dp) :: value

      self%state = modulo(multiplier * self%state, modulus)
      value = real(self%state, dp) / real(modulus, dp)
   end function rng_uniform

   function rng_normal(self) result(value)
      class(gnorm_rng_state), intent(inout) :: self
      real(dp) :: value
      real(dp) :: radius, angle, u1, u2
      real(dp), parameter :: two_pi = 2.0_dp * acos(-1.0_dp)

      if (self%has_spare) then
         value = self%spare
         self%has_spare = .false.
         return
      end if
      u1 = max(self%uniform(), tiny(1.0_dp))
      u2 = self%uniform()
      radius = sqrt(-2.0_dp * log(u1))
      angle = two_pi * u2
      value = radius * cos(angle)
      self%spare = radius * sin(angle)
      self%has_spare = .true.
   end function rng_normal

   recursive function rng_gamma(self, shape) result(value)
      class(gnorm_rng_state), intent(inout) :: self
      real(dp), intent(in) :: shape
      real(dp) :: value
      real(dp) :: d, c, x, v, u

      if (shape <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (shape < 1.0_dp) then
         u = max(self%uniform(), tiny(1.0_dp))
         value = self%gamma(shape + 1.0_dp) * u**(1.0_dp / shape)
         return
      end if

      d = shape - 1.0_dp / 3.0_dp
      c = 1.0_dp / sqrt(9.0_dp * d)
      do
         x = self%normal()
         v = 1.0_dp + c * x
         if (v <= 0.0_dp) cycle
         v = v * v * v
         u = self%uniform()
         if (u < 1.0_dp - 0.0331_dp * x**4) exit
         if (log(max(u, tiny(1.0_dp))) < 0.5_dp * x*x + d*(1.0_dp - v + log(v))) exit
      end do
      value = d * v
   end function rng_gamma

end module gnorm_rng
