! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_rng
   use, intrinsic :: iso_fortran_env, only : int64
   use dfoptim_kinds, only : dp
   implicit none
   private

   type, public :: rng_t
      integer(int64) :: state = 1138_int64
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: permutation => rng_permutation
   end type rng_t

contains

   subroutine rng_seed(self, seed)
      class(rng_t), intent(inout) :: self
      integer, intent(in) :: seed
      integer(int64), parameter :: modulus = 2147483647_int64

      self%state = modulo(abs(int(seed, int64)), modulus - 1_int64) + 1_int64
   end subroutine rng_seed

   function rng_uniform(self) result(u)
      class(rng_t), intent(inout) :: self
      real(dp) :: u
      integer(int64), parameter :: modulus = 2147483647_int64
      integer(int64), parameter :: multiplier = 16807_int64
      integer(int64), parameter :: quotient = 127773_int64
      integer(int64), parameter :: remainder = 2836_int64
      integer(int64) :: high, low, test

      high = self%state / quotient
      low = modulo(self%state, quotient)
      test = multiplier * low - remainder * high
      if (test <= 0_int64) test = test + modulus
      self%state = test
      u = real(self%state, dp) / real(modulus, dp)
   end function rng_uniform

   subroutine rng_permutation(self, p)
      class(rng_t), intent(inout) :: self
      integer, intent(out) :: p(:)
      integer :: i, j, tmp

      p = [(i, i = 1, size(p))]
      do i = size(p), 2, -1
         j = 1 + int(self%uniform() * real(i, dp))
         j = min(i, max(1, j))
         tmp = p(i)
         p(i) = p(j)
         p(j) = tmp
      end do
   end subroutine rng_permutation

end module dfoptim_rng
