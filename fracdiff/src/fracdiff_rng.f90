! SPDX-License-Identifier: GPL-2.0-or-later
module fracdiff_rng
   use, intrinsic :: iso_fortran_env, only : int64
   use fracdiff_kinds, only : dp, pi_dp
   implicit none
   private

   type, public :: fracdiff_rng_state
      integer(int64) :: state = 1_int64
      logical :: has_spare = .false.
      real(dp) :: spare = 0.0_dp
   end type fracdiff_rng_state

   public :: seed_rng, random_uniform, random_normal, fill_normal

contains

   subroutine seed_rng(rng, seed)
      type(fracdiff_rng_state), intent(inout) :: rng
      integer(int64), intent(in) :: seed
      integer(int64), parameter :: modulus = 2147483647_int64

      rng%state = modulo(abs(seed), modulus - 1_int64) + 1_int64
      rng%has_spare = .false.
      rng%spare = 0.0_dp
   end subroutine seed_rng

   function random_uniform(rng) result(value)
      type(fracdiff_rng_state), intent(inout) :: rng
      real(dp) :: value
      integer(int64), parameter :: multiplier = 16807_int64
      integer(int64), parameter :: modulus = 2147483647_int64
      integer(int64), parameter :: quotient = 127773_int64
      integer(int64), parameter :: remainder = 2836_int64
      integer(int64) :: high, low, test

      high = rng%state/quotient
      low = modulo(rng%state, quotient)
      test = multiplier*low - remainder*high
      if (test > 0_int64) then
         rng%state = test
      else
         rng%state = test + modulus
      end if
      value = real(rng%state, dp)/real(modulus, dp)
   end function random_uniform

   function random_normal(rng) result(value)
      type(fracdiff_rng_state), intent(inout) :: rng
      real(dp) :: value
      real(dp) :: u1, u2, radius, angle

      if (rng%has_spare) then
         value = rng%spare
         rng%has_spare = .false.
         return
      end if

      u1 = max(random_uniform(rng), tiny(1.0_dp))
      u2 = random_uniform(rng)
      radius = sqrt(-2.0_dp*log(u1))
      angle = 2.0_dp*pi_dp*u2
      value = radius*cos(angle)
      rng%spare = radius*sin(angle)
      rng%has_spare = .true.
   end function random_normal

   subroutine fill_normal(rng, values, standard_deviation)
      type(fracdiff_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: values(:)
      real(dp), intent(in), optional :: standard_deviation
      real(dp) :: sd
      integer :: i

      sd = 1.0_dp
      if (present(standard_deviation)) sd = standard_deviation
      do i = 1, size(values)
         values(i) = sd*random_normal(rng)
      end do
   end subroutine fill_normal

end module fracdiff_rng
