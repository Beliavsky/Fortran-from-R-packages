! Bit operations used by the WELL recurrences.
! See LICENSES/WELL-SOURCE-NOTICE.txt for the notice carried by the upstream
! WELL C sources from Francois Panneton, Pierre L'Ecuyer and Makoto Matsumoto.
module rngwell_bitops
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   implicit none
   private

   integer(int64), parameter :: two32 = 4294967296_int64
   real(real64), parameter :: inv_two32 = 2.3283064365386962890625e-10_real64

   public :: lxor, rxor, lshift32, rshift32, mat2, mat4neg, mat5
   public :: u32_to_real, u32_to_i64, i64_to_u32, inv_two32

contains

   pure elemental integer(int32) function lxor(v, k) result(y)
      integer(int32), intent(in) :: v
      integer, intent(in) :: k
      y = ieor(v, shiftl(v, k))
   end function lxor

   pure elemental integer(int32) function rxor(v, k) result(y)
      integer(int32), intent(in) :: v
      integer, intent(in) :: k
      y = ieor(v, shiftr(v, k))
   end function rxor

   pure elemental integer(int32) function lshift32(v, k) result(y)
      integer(int32), intent(in) :: v
      integer, intent(in) :: k
      y = shiftl(v, k)
   end function lshift32

   pure elemental integer(int32) function rshift32(v, k) result(y)
      integer(int32), intent(in) :: v
      integer, intent(in) :: k
      y = shiftr(v, k)
   end function rshift32

   pure elemental integer(int32) function mat2(a, v) result(y)
      integer(int32), intent(in) :: a, v
      if (btest(v, 0)) then
         y = ieor(shiftr(v, 1), a)
      else
         y = shiftr(v, 1)
      end if
   end function mat2

   pure elemental integer(int32) function mat4neg(k, b, v) result(y)
      integer, intent(in) :: k
      integer(int32), intent(in) :: b, v
      y = ieor(v, iand(shiftl(v, k), b))
   end function mat4neg

   pure elemental integer(int32) function mat5(r, a, ds, dt, v) result(y)
      integer, intent(in) :: r
      integer(int32), intent(in) :: a, ds, dt, v
      integer(int32) :: t

      t = iand(ieor(shiftl(v, r), shiftr(v, 32-r)), ds)
      if (iand(v, dt) /= 0_int32) then
         y = ieor(t, a)
      else
         y = t
      end if
   end function mat5

   pure elemental integer(int64) function u32_to_i64(x) result(y)
      integer(int32), intent(in) :: x
      y = int(x, int64)
      if (y < 0_int64) y = y + two32
   end function u32_to_i64

   pure elemental integer(int32) function i64_to_u32(x) result(y)
      integer(int64), intent(in) :: x
      integer(int64) :: z

      z = modulo(x, two32)
      if (z <= int(huge(0_int32), int64)) then
         y = int(z, int32)
      else
         y = int(z - two32, int32)
      end if
   end function i64_to_u32

   pure elemental real(real64) function u32_to_real(x) result(y)
      integer(int32), intent(in) :: x
      y = real(u32_to_i64(x), real64) * inv_two32
   end function u32_to_real

end module rngwell_bitops
