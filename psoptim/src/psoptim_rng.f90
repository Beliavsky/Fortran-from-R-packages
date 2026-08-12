! psoptim-fortran: modern Fortran translation of psoptim 1.0
! Original package copyright: Krzysztof Ciupke
! License: GPL-2.0-or-later. See COPYING and original/DESCRIPTION.
module psoptim_rng
   use, intrinsic :: iso_fortran_env, only : int64, real64
   implicit none
   private

   integer, parameter, public :: dp = real64

   type, public :: ps_rng
      integer(int64) :: state = 1_int64
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
   end type ps_rng

contains

   subroutine rng_seed(self, seed)
      class(ps_rng), intent(inout) :: self
      integer(int64), intent(in) :: seed

      if (seed == 0_int64) then
         self%state = 1_int64
      else
         self%state = seed
      end if
   end subroutine rng_seed

   function rng_uniform(self) result(u)
      class(ps_rng), intent(inout) :: self
      real(dp) :: u
      integer(int64) :: z

      ! xorshift64*; deterministic and independent of compiler runtime RNG.
      z = self%state
      z = ieor(z, shiftl(z, 13))
      z = ieor(z, shiftr(z, 7))
      z = ieor(z, shiftl(z, 17))
      self%state = z

      ! Use the low 53 bits after masking off the sign bit.
      z = iand(z, int(z'001fffffffffffff', int64))
      u = real(z, dp) / real(int(z'0020000000000000', int64), dp)
   end function rng_uniform

end module psoptim_rng
