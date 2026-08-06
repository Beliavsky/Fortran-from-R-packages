! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_random
   use iso_fortran_env, only : int64
   use spillover_kinds, only : dp
   implicit none
   private

   integer(int64), parameter :: modulus = 2147483647_int64
   integer(int64), parameter :: multiplier = 16807_int64

   type, public :: spillover_rng
      integer(int64) :: state = 1_int64
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: permutation => rng_permutation
   end type spillover_rng

contains

   subroutine rng_seed(self, seed)
      class(spillover_rng), intent(inout) :: self
      integer, intent(in) :: seed

      self%state = modulo(int(seed, int64), modulus - 1_int64) + 1_int64
   end subroutine rng_seed

   function rng_uniform(self) result(u)
      class(spillover_rng), intent(inout) :: self
      real(dp) :: u

      self%state = modulo(multiplier * self%state, modulus)
      u = real(self%state, dp) / real(modulus, dp)
   end function rng_uniform

   subroutine rng_permutation(self, perm)
      class(spillover_rng), intent(inout) :: self
      integer, intent(out) :: perm(:)

      integer :: i, j, tmp

      perm = [(i, i = 1, size(perm))]
      do i = size(perm), 2, -1
         j = 1 + int(self%uniform() * real(i, dp))
         if (j > i) j = i
         tmp = perm(i)
         perm(i) = perm(j)
         perm(j) = tmp
      end do
   end subroutine rng_permutation

end module spillover_random
