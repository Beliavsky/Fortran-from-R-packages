module gmp_bigq
   use, intrinsic :: iso_fortran_env, only : int64
   use gmp_kinds, only : dp
   use gmp_bigz, only : bigz, bigz_from_int64, bigz_from_string, bigz_to_string, bigz_to_dp, &
      bigz_is_zero, bigz_is_one, bigz_compare, bigz_equal, bigz_add, bigz_sub, bigz_neg, bigz_abs, &
      bigz_mul, bigz_gcd, bigz_quotient, bigz_pow, bigz_modulo
   implicit none
   private

   type, public :: bigq
      type(bigz) :: num
      type(bigz) :: den
   end type bigq

   public :: bigq_make, bigq_from_int64, bigq_from_string, bigq_to_string, bigq_to_dp
   public :: bigq_add, bigq_sub, bigq_neg, bigq_abs, bigq_mul, bigq_div, bigq_pow
   public :: bigq_compare, bigq_equal, bigq_is_zero, bigq_is_integer
   public :: bigq_floor, bigq_ceiling, bigq_trunc, bigq_round0, bigq_round_digits
   public :: bigq_numerator, bigq_denominator

contains

   pure elemental function bigq_make(numerator, denominator) result(q)
      type(bigz), intent(in) :: numerator !! Exact rational numerator before reduction.
      type(bigz), intent(in) :: denominator !! Exact nonzero rational denominator before reduction.
      type(bigq) :: q
      type(bigz) :: g, n, d

      if (bigz_is_zero(denominator)) then
         q%num = bigz_from_int64(0_int64)
         q%den = bigz_from_int64(0_int64)
         return
      end if
      if (bigz_is_zero(numerator)) then
         q%num = bigz_from_int64(0_int64)
         q%den = bigz_from_int64(1_int64)
         return
      end if
      n = numerator
      d = denominator
      if (d%sign < 0) then
         n = bigz_neg(n)
         d = bigz_neg(d)
      end if
      g = bigz_gcd(bigz_abs(n), d)
      q%num = bigz_quotient(n, g)
      q%den = bigz_quotient(d, g)
   end function bigq_make

   pure elemental function bigq_from_int64(value) result(q)
      integer(int64), intent(in) :: value !! Signed machine integer converted exactly to a rational with denominator one.
      type(bigq) :: q
      q = bigq_make(bigz_from_int64(value), bigz_from_int64(1_int64))
   end function bigq_from_int64

   pure function bigq_from_string(text) result(q)
      character(len=*), intent(in) :: text !! Integer, decimal, or numerator/denominator character representation.
      type(bigq) :: q
      character(len=:), allocatable :: s, left, right, digits
      integer :: slash, point, i, nfrac
      logical :: valid
      type(bigz) :: n, d

      s = trim(adjustl(text))
      valid = len(s) > 0
      slash = index(s, '/')
      point = index(s, '.')
      if (.not. valid) then
         q = bigq_from_int64(0_int64)
      else if (slash > 0) then
         if (index(s(slash + 1:), '/') > 0) then
            valid = .false.
            q = bigq_from_int64(0_int64)
         else
            left = s(:slash - 1)
            right = s(slash + 1:)
            n = bigz_from_string(left)
            d = bigz_from_string(right)
            if (bigz_is_zero(d)) valid = .false.
            q = bigq_make(n, d)
         end if
      else if (point > 0) then
         if (index(s(point + 1:), '.') > 0) then
            valid = .false.
            q = bigq_from_int64(0_int64)
         else
            left = s(:point - 1)
            right = s(point + 1:)
            nfrac = len(right)
            digits = left // right
            if (len_trim(digits) == 0 .or. digits == '+' .or. digits == '-') digits = digits // '0'
            n = bigz_from_string(digits)
            d = bigz_from_int64(1_int64)
            do i = 1, nfrac
               d = bigz_mul(d, bigz_from_int64(10_int64))
            end do
            q = bigq_make(n, d)
         end if
      else
         q = bigq_make(bigz_from_string(s), bigz_from_int64(1_int64))
      end if
   end function bigq_from_string

   pure function bigq_to_string(q, base) result(text)
      type(bigq), intent(in) :: q !! Exact rational value to format.
      integer, intent(in), optional :: base !! Integer radix from 2 through 36 used for numerator and denominator.
      character(len=:), allocatable :: text
      integer :: radix
      character(len=:), allocatable :: ns, ds

      radix = 10
      if (present(base)) radix = base
      if (bigz_is_zero(q%den)) then
         text = 'NaN'
         return
      end if
      ns = bigz_to_string(q%num, radix)
      if (bigz_is_one(q%den)) then
         text = ns
      else
         ds = bigz_to_string(q%den, radix)
         text = ns // '/' // ds
      end if
   end function bigq_to_string

   pure elemental real(dp) function bigq_to_dp(q) result(value)
      type(bigq), intent(in) :: q !! Exact rational converted approximately to the package real kind.
      if (bigz_is_zero(q%den)) then
         value = huge(1.0_dp)
      else
         value = bigz_to_dp(q%num) / bigz_to_dp(q%den)
      end if
   end function bigq_to_dp

   pure elemental function bigq_add(a, b) result(c)
      type(bigq), intent(in) :: a !! Left rational addend.
      type(bigq), intent(in) :: b !! Right rational addend.
      type(bigq) :: c
      type(bigz) :: n, d, g, ad, bd

      if (bigz_is_zero(a%den) .or. bigz_is_zero(b%den)) then
         c%num = bigz_from_int64(0_int64)
         c%den = bigz_from_int64(0_int64)
         return
      end if
      g = bigz_gcd(a%den, b%den)
      ad = bigz_quotient(a%den, g)
      bd = bigz_quotient(b%den, g)
      n = bigz_add(bigz_mul(a%num, bd), bigz_mul(b%num, ad))
      d = bigz_mul(ad, b%den)
      c = bigq_make(n, d)
   end function bigq_add

   pure elemental function bigq_sub(a, b) result(c)
      type(bigq), intent(in) :: a !! Rational minuend.
      type(bigq), intent(in) :: b !! Rational subtrahend.
      type(bigq) :: c
      c = bigq_add(a, bigq_neg(b))
   end function bigq_sub

   pure elemental function bigq_neg(a) result(c)
      type(bigq), intent(in) :: a !! Rational whose sign is reversed.
      type(bigq) :: c
      c = a
      c%num = bigz_neg(c%num)
   end function bigq_neg

   pure elemental function bigq_abs(a) result(c)
      type(bigq), intent(in) :: a !! Rational whose nonnegative magnitude is returned.
      type(bigq) :: c
      c = a
      c%num = bigz_abs(c%num)
   end function bigq_abs

   pure elemental function bigq_mul(a, b) result(c)
      type(bigq), intent(in) :: a !! Left rational multiplicand.
      type(bigq), intent(in) :: b !! Right rational multiplicand.
      type(bigq) :: c
      type(bigz) :: g1, g2, an, bn, ad, bd

      if (bigz_is_zero(a%den) .or. bigz_is_zero(b%den)) then
         c%num = bigz_from_int64(0_int64)
         c%den = bigz_from_int64(0_int64)
         return
      end if
      g1 = bigz_gcd(bigz_abs(a%num), b%den)
      g2 = bigz_gcd(bigz_abs(b%num), a%den)
      an = bigz_quotient(a%num, g1)
      bd = bigz_quotient(b%den, g1)
      bn = bigz_quotient(b%num, g2)
      ad = bigz_quotient(a%den, g2)
      c = bigq_make(bigz_mul(an, bn), bigz_mul(ad, bd))
   end function bigq_mul

   pure elemental function bigq_div(a, b) result(c)
      type(bigq), intent(in) :: a !! Rational dividend.
      type(bigq), intent(in) :: b !! Nonzero rational divisor.
      type(bigq) :: c

      if (bigz_is_zero(a%den) .or. bigz_is_zero(b%den) .or. bigz_is_zero(b%num)) then
         c%num = bigz_from_int64(0_int64)
         c%den = bigz_from_int64(0_int64)
         return
      end if
      c = bigq_make(bigz_mul(a%num, b%den), bigz_mul(a%den, b%num))
   end function bigq_div

   pure elemental function bigq_pow(a, exponent) result(c)
      type(bigq), intent(in) :: a !! Rational base.
      integer(int64), intent(in) :: exponent !! Signed machine-integer exponent.
      type(bigq) :: c
      integer(int64) :: e

      if (exponent == 0_int64) then
         c = bigq_from_int64(1_int64)
         return
      end if
      if (exponent == int(z'8000000000000000', int64)) then
         c%num = bigz_from_int64(0_int64)
         c%den = bigz_from_int64(0_int64)
         return
      end if
      e = abs(exponent)
      if (exponent > 0_int64) then
         c = bigq_make(bigz_pow(a%num, e), bigz_pow(a%den, e))
      else
         c = bigq_make(bigz_pow(a%den, e), bigz_pow(a%num, e))
      end if
   end function bigq_pow

   pure elemental integer function bigq_compare(a, b) result(cmp)
      type(bigq), intent(in) :: a !! Left rational in the comparison.
      type(bigq), intent(in) :: b !! Right rational in the comparison.
      type(bigz) :: left, right

      if (bigz_is_zero(a%den) .or. bigz_is_zero(b%den)) then
         cmp = 0
         return
      end if
      left = bigz_mul(a%num, b%den)
      right = bigz_mul(b%num, a%den)
      cmp = bigz_compare(left, right)
   end function bigq_compare

   pure elemental logical function bigq_equal(a, b) result(equal)
      type(bigq), intent(in) :: a !! Left rational value.
      type(bigq), intent(in) :: b !! Right rational value.
      equal = bigq_compare(a, b) == 0 .and. .not. bigz_is_zero(a%den) .and. .not. bigz_is_zero(b%den)
   end function bigq_equal

   pure elemental logical function bigq_is_zero(a) result(is_zero)
      type(bigq), intent(in) :: a !! Rational tested for exact equality to zero.
      is_zero = .not. bigz_is_zero(a%den) .and. bigz_is_zero(a%num)
   end function bigq_is_zero

   pure elemental logical function bigq_is_integer(a) result(is_integer)
      type(bigq), intent(in) :: a !! Rational tested for an exact integer value.
      is_integer = .not. bigz_is_zero(a%den) .and. bigz_is_one(a%den)
   end function bigq_is_integer

   pure elemental function bigq_floor(a) result(z)
      type(bigq), intent(in) :: a !! Rational whose mathematical floor is returned.
      type(bigz) :: z
      if (bigz_is_zero(a%den)) then
         z = bigz_from_int64(0_int64)
      else
         z = bigz_quotient(a%num, a%den)
      end if
   end function bigq_floor

   pure elemental function bigq_ceiling(a) result(z)
      type(bigq), intent(in) :: a !! Rational whose mathematical ceiling is returned.
      type(bigz) :: z
      z = bigz_neg(bigq_floor(bigq_neg(a)))
   end function bigq_ceiling

   pure elemental function bigq_trunc(a) result(z)
      type(bigq), intent(in) :: a !! Rational truncated toward zero.
      type(bigz) :: z, q

      if (bigz_is_zero(a%den)) then
         z = bigz_from_int64(0_int64)
         return
      end if
      q = bigz_quotient(bigz_abs(a%num), a%den)
      if (a%num%sign < 0) q = bigz_neg(q)
      z = q
   end function bigq_trunc

   pure elemental function bigq_round0(a) result(z)
      type(bigq), intent(in) :: a !! Rational rounded to nearest integer with ties to even, matching R round().
      type(bigz) :: z, floor_z, rem, twice_rem
      type(bigq) :: abs_a
      integer :: cmp

      if (bigz_is_zero(a%den)) then
         z = bigz_from_int64(0_int64)
         return
      end if
      abs_a = bigq_abs(a)
      floor_z = bigq_floor(abs_a)
      rem = bigz_sub(abs_a%num, bigz_mul(floor_z, abs_a%den))
      twice_rem = bigz_add(rem, rem)
      cmp = bigz_compare(twice_rem, abs_a%den)
      z = floor_z
      if (cmp > 0) then
         z = bigz_add(z, bigz_from_int64(1_int64))
      else if (cmp == 0) then
         if (.not. bigz_is_zero(bigz_modulo(z, bigz_from_int64(2_int64)))) then
            z = bigz_add(z, bigz_from_int64(1_int64))
         end if
      end if
      if (a%num%sign < 0) z = bigz_neg(z)
   end function bigq_round0

   pure elemental function bigq_round_digits(a, digits) result(c)
      type(bigq), intent(in) :: a !! Rational value to round.
      integer, intent(in) :: digits !! Decimal digits; positive rounds after the point and negative before it.
      type(bigq) :: c, scaled
      type(bigz) :: p10, rounded
      integer :: i

      p10 = bigz_from_int64(1_int64)
      do i = 1, abs(digits)
         p10 = bigz_mul(p10, bigz_from_int64(10_int64))
      end do
      if (digits >= 0) then
         scaled = bigq_mul(a, bigq_make(p10, bigz_from_int64(1_int64)))
         rounded = bigq_round0(scaled)
         c = bigq_make(rounded, p10)
      else
         scaled = bigq_div(a, bigq_make(p10, bigz_from_int64(1_int64)))
         rounded = bigq_round0(scaled)
         c = bigq_make(bigz_mul(rounded, p10), bigz_from_int64(1_int64))
      end if
   end function bigq_round_digits

   pure elemental function bigq_numerator(a) result(z)
      type(bigq), intent(in) :: a !! Rational whose reduced numerator is returned.
      type(bigz) :: z
      z = a%num
   end function bigq_numerator

   pure elemental function bigq_denominator(a) result(z)
      type(bigq), intent(in) :: a !! Rational whose positive reduced denominator is returned.
      type(bigz) :: z
      z = a%den
   end function bigq_denominator

end module gmp_bigq
