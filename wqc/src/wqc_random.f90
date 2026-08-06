! SPDX-License-Identifier: GPL-3.0-only
module wqc_random
   use wqc_kinds, only : dp
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   private

   integer(int64), parameter :: modulus = 2147483647_int64
   integer(int64), parameter :: multiplier = 16807_int64
   integer(int64), parameter :: quotient = 127773_int64
   integer(int64), parameter :: remainder = 2836_int64
   real(dp), parameter :: two_pi = 6.28318530717958647692528676655900577_dp

   type, public :: normal_rng
      private
      integer(int64) :: state = 1_int64
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   contains
      procedure :: seed => seed_rng
      procedure :: uniform => next_uniform
      procedure :: normal => next_normal
      procedure :: fill_normal
   end type normal_rng

contains

   subroutine seed_rng(self, seed)
      class(normal_rng), intent(inout) :: self
      integer, intent(in), optional :: seed

      integer :: values(8)
      integer(int64) :: raw

      if (present(seed)) then
         raw = int(seed, int64)
      else
         call date_and_time(values=values)
         raw = int(values(1), int64)
         raw = 37_int64 * raw + int(values(2), int64)
         raw = 37_int64 * raw + int(values(3), int64)
         raw = 37_int64 * raw + int(values(5), int64)
         raw = 37_int64 * raw + int(values(6), int64)
         raw = 37_int64 * raw + int(values(7), int64)
         raw = 1000_int64 * raw + int(values(8), int64)
      end if

      self%state = modulo(abs(raw), modulus - 1_int64) + 1_int64
      self%has_spare = .false.
      self%spare = 0.0_dp
   end subroutine seed_rng


   function next_uniform(self) result(u)
      class(normal_rng), intent(inout) :: self
      real(dp) :: u

      integer(int64) :: hi, lo, test

      hi = self%state / quotient
      lo = modulo(self%state, quotient)
      test = multiplier * lo - remainder * hi
      if (test > 0_int64) then
         self%state = test
      else
         self%state = test + modulus
      end if
      u = real(self%state, dp) / real(modulus, dp)
   end function next_uniform


   function next_normal(self) result(z)
      class(normal_rng), intent(inout) :: self
      real(dp) :: z

      real(dp) :: radius, theta, u1, u2

      if (self%has_spare) then
         z = self%spare
         self%has_spare = .false.
         return
      end if

      u1 = max(self%uniform(), tiny(1.0_dp))
      u2 = self%uniform()
      radius = sqrt(-2.0_dp * log(u1))
      theta = two_pi * u2
      z = radius * cos(theta)
      self%spare = radius * sin(theta)
      self%has_spare = .true.
   end function next_normal


   subroutine fill_normal(self, x, location, scale)
      class(normal_rng), intent(inout) :: self
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: location, scale

      integer :: i

      do i = 1, size(x)
         x(i) = location + scale * self%normal()
      end do
   end subroutine fill_normal

end module wqc_random
