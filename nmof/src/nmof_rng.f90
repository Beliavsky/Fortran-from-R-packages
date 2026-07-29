! SPDX-License-Identifier: GPL-3.0-only
module nmof_rng
   use nmof_kinds, only: dp, i8, pi
   implicit none
   private
   public :: rng_state, rng_seed, rng_uniform, rng_normal, rng_integer
   public :: rng_shuffle, rng_logical

   type :: rng_state
      integer(i8) :: state = 88172645463393265_i8
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type rng_state
contains
   subroutine rng_seed(rng, seed)
      type(rng_state), intent(inout) :: rng
      integer(i8), intent(in) :: seed
      integer(i8) :: z
      z = seed
      if (z == 0_i8) z = 88172645463393265_i8
      z = ieor(z, shiftr(z, 30)) * int(z'BF58476D1CE4E5B9', i8)
      z = ieor(z, shiftr(z, 27)) * int(z'94D049BB133111EB', i8)
      rng%state = ieor(z, shiftr(z, 31))
      if (rng%state == 0_i8) rng%state = 88172645463393265_i8
      rng%has_spare = .false.
   end subroutine rng_seed

   function next_uint64(rng) result(x)
      type(rng_state), intent(inout) :: rng
      integer(i8) :: x
      x = rng%state
      x = ieor(x, shiftl(x, 13))
      x = ieor(x, shiftr(x, 7))
      x = ieor(x, shiftl(x, 17))
      rng%state = x
   end function next_uint64

   function rng_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(i8) :: x
      x = iand(shiftr(next_uint64(rng), 11), int(z'001FFFFFFFFFFFFF', i8))
      u = real(x, dp) * 1.11022302462515654042e-16_dp
      if (u <= 0.0_dp) u = epsilon(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
   end function rng_uniform

   function rng_normal(rng) result(z)
      type(rng_state), intent(inout) :: rng
      real(dp) :: z, r, theta
      if (rng%has_spare) then
         z = rng%spare
         rng%has_spare = .false.
      else
         r = sqrt(-2.0_dp * log(rng_uniform(rng)))
         theta = 2.0_dp * pi * rng_uniform(rng)
         z = r * cos(theta)
         rng%spare = r * sin(theta)
         rng%has_spare = .true.
      end if
   end function rng_normal

   function rng_integer(rng, n) result(k)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: n
      integer :: k
      if (n <= 1) then
         k = 1
      else
         k = 1 + int(rng_uniform(rng) * real(n, dp))
         if (k > n) k = n
      end if
   end function rng_integer

   function rng_logical(rng, probability) result(value)
      type(rng_state), intent(inout) :: rng
      real(dp), intent(in) :: probability
      logical :: value
      value = rng_uniform(rng) < probability
   end function rng_logical

   subroutine rng_shuffle(rng, indices)
      type(rng_state), intent(inout) :: rng
      integer, intent(inout) :: indices(:)
      integer :: i, j, tmp
      do i = size(indices), 2, -1
         j = rng_integer(rng, i)
         tmp = indices(i)
         indices(i) = indices(j)
         indices(j) = tmp
      end do
   end subroutine rng_shuffle
end module nmof_rng
