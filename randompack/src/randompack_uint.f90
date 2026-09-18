module randompack_uint
   use iso_fortran_env, only : int64
   implicit none
   private
   integer, parameter :: i128 = selected_int_kind(38)
   integer(int64), parameter :: mask32 = int(z'00000000FFFFFFFF', int64)
   integer(i128), parameter :: base32 = 4294967296_i128
   public :: uadd, usub, umul, umul_wide, urotl, add128, mul128_mod
   public :: u64_from_u32, lo32, hi32, unsigned_lt
contains
   pure elemental function lo32(x) result(y)
      integer(int64), intent(in) :: x !! Packed 64-bit word whose low 32 bits are requested.
      integer(int64) :: y
      y = iand(x, mask32)
   end function lo32

   pure elemental function hi32(x) result(y)
      integer(int64), intent(in) :: x !! Packed 64-bit word whose high 32 bits are requested.
      integer(int64) :: y
      y = shiftr(x, 32)
   end function hi32

   pure elemental function u64_from_u32(low, high) result(x)
      integer(int64), intent(in) :: low !! Unsigned low 32-bit half stored in an int64.
      integer(int64), intent(in) :: high !! Unsigned high 32-bit half stored in an int64.
      integer(int64) :: x
      x = ior(iand(low, mask32), shiftl(iand(high, mask32), 32))
   end function u64_from_u32

   pure elemental function uadd(a, b) result(c)
      integer(int64), intent(in) :: a !! First unsigned 64-bit addend represented by its bits.
      integer(int64), intent(in) :: b !! Second unsigned 64-bit addend represented by its bits.
      integer(int64) :: c
      integer(int64) :: alo, ahi, blo, bhi, low, high, carry
      alo = lo32(a)
      ahi = hi32(a)
      blo = lo32(b)
      bhi = hi32(b)
      low = alo + blo
      carry = shiftr(low, 32)
      high = iand(ahi + bhi + carry, mask32)
      c = u64_from_u32(iand(low, mask32), high)
   end function uadd

   pure elemental function usub(a, b) result(c)
      integer(int64), intent(in) :: a !! Unsigned 64-bit minuend represented by its bits.
      integer(int64), intent(in) :: b !! Unsigned 64-bit subtrahend represented by its bits.
      integer(int64) :: c
      c = uadd(a, uadd(not(b), 1_int64))
   end function usub

   pure elemental function umul(a, b) result(c)
      integer(int64), intent(in) :: a !! First unsigned 64-bit factor represented by its bits.
      integer(int64), intent(in) :: b !! Second unsigned 64-bit factor represented by its bits.
      integer(int64) :: c
      integer(i128) :: a0, a1, b0, b1, p00, cross, low128
      integer(int64) :: low, high
      a0 = int(lo32(a), i128)
      a1 = int(hi32(a), i128)
      b0 = int(lo32(b), i128)
      b1 = int(hi32(b), i128)
      p00 = a0 * b0
      cross = a0 * b1 + a1 * b0
      low128 = p00 + modulo(cross, base32) * base32
      low = int(modulo(low128, base32), int64)
      high = int(modulo(low128 / base32, base32), int64)
      c = u64_from_u32(low, high)
   end function umul

   pure subroutine umul_wide(a, b, high, low)
      integer(int64), intent(in) :: a !! First unsigned 64-bit factor represented by its bits.
      integer(int64), intent(in) :: b !! Second unsigned 64-bit factor represented by its bits.
      integer(int64), intent(out) :: high !! High 64 bits of the exact 128-bit unsigned product.
      integer(int64), intent(out) :: low !! Low 64 bits of the exact 128-bit unsigned product.
      integer(i128) :: a0, a1, b0, b1, p00, p01, p10, p11, mid, hi128
      integer(int64) :: l0, l1, h0, h1
      a0 = int(lo32(a), i128)
      a1 = int(hi32(a), i128)
      b0 = int(lo32(b), i128)
      b1 = int(hi32(b), i128)
      p00 = a0 * b0
      p01 = a0 * b1
      p10 = a1 * b0
      p11 = a1 * b1
      l0 = int(modulo(p00, base32), int64)
      mid = p00 / base32 + modulo(p01, base32) + modulo(p10, base32)
      l1 = int(modulo(mid, base32), int64)
      hi128 = p11 + p01 / base32 + p10 / base32 + mid / base32
      h0 = int(modulo(hi128, base32), int64)
      h1 = int(modulo(hi128 / base32, base32), int64)
      low = u64_from_u32(l0, l1)
      high = u64_from_u32(h0, h1)
   end subroutine umul_wide

   pure elemental function urotl(x, k) result(y)
      integer(int64), intent(in) :: x !! Unsigned 64-bit word represented by its bits.
      integer, intent(in) :: k !! Left rotation distance in bits, normally in 0..63.
      integer(int64) :: y
      integer :: r
      r = modulo(k, 64)
      if (r == 0) then
         y = x
      else
         y = ior(shiftl(x, r), shiftr(x, 64 - r))
      end if
   end function urotl

   pure subroutine add128(a_lo, a_hi, b_lo, b_hi, c_lo, c_hi)
      integer(int64), intent(in) :: a_lo !! Low word of the first unsigned 128-bit addend.
      integer(int64), intent(in) :: a_hi !! High word of the first unsigned 128-bit addend.
      integer(int64), intent(in) :: b_lo !! Low word of the second unsigned 128-bit addend.
      integer(int64), intent(in) :: b_hi !! High word of the second unsigned 128-bit addend.
      integer(int64), intent(out) :: c_lo !! Low word of the modulo-2^128 result.
      integer(int64), intent(out) :: c_hi !! High word of the modulo-2^128 result.
      logical :: carry
      c_lo = uadd(a_lo, b_lo)
      carry = unsigned_lt(c_lo, a_lo)
      c_hi = uadd(a_hi, b_hi)
      if (carry) c_hi = uadd(c_hi, 1_int64)
   end subroutine add128

   pure subroutine mul128_mod(a_lo, a_hi, b_lo, b_hi, c_lo, c_hi)
      integer(int64), intent(in) :: a_lo !! Low word of the first unsigned 128-bit factor.
      integer(int64), intent(in) :: a_hi !! High word of the first unsigned 128-bit factor.
      integer(int64), intent(in) :: b_lo !! Low word of the second unsigned 128-bit factor.
      integer(int64), intent(in) :: b_hi !! High word of the second unsigned 128-bit factor.
      integer(int64), intent(out) :: c_lo !! Low word of the product modulo 2^128.
      integer(int64), intent(out) :: c_hi !! High word of the product modulo 2^128.
      integer(int64) :: hh, ll
      call umul_wide(a_lo, b_lo, hh, ll)
      c_lo = ll
      c_hi = uadd(hh, uadd(umul(a_lo, b_hi), umul(a_hi, b_lo)))
   end subroutine mul128_mod

   pure elemental logical function unsigned_lt(a, b)
      integer(int64), intent(in) :: a !! First unsigned 64-bit operand represented by its bits.
      integer(int64), intent(in) :: b !! Second unsigned 64-bit operand represented by its bits.
      integer(int64), parameter :: signbit = int(z'8000000000000000', int64)
      unsigned_lt = (ieor(a, signbit) < ieor(b, signbit))
   end function unsigned_lt
end module randompack_uint
