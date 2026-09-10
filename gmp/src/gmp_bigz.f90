module gmp_bigz
   use, intrinsic :: iso_fortran_env, only : int32, int64
   use gmp_kinds, only : dp
   implicit none
   private

   integer(int64), parameter :: limb_base = 10000_int64
   integer, parameter :: limb_width = 4

   type, public :: bigz
      integer :: sign = 0
      integer(int32), allocatable :: limb(:)
   end type bigz

   public :: bigz_from_int64, bigz_from_string, bigz_to_string, bigz_to_dp
   public :: bigz_compare, bigz_equal, bigz_is_zero, bigz_is_one, bigz_is_even
   public :: bigz_add, bigz_sub, bigz_neg, bigz_abs, bigz_mul, bigz_mul_small
   public :: bigz_divmod, bigz_quotient, bigz_modulo, bigz_div_small
   public :: bigz_gcd, bigz_lcm, bigz_gcdex, bigz_pow, bigz_powm, bigz_invmod
   public :: bigz_sizeinbase, bigz_frexp, bigz_log, bigz_log10, bigz_log2
   public :: bigz_to_int64, bigz_fits_int64

contains

   pure recursive function bigz_from_int64(value) result(z)
      integer(int64), intent(in) :: value !! Signed machine integer to convert exactly to arbitrary precision.
      type(bigz) :: z
      integer(int64) :: mag
      integer :: n

      if (value == 0_int64) return
      if (value < 0_int64) then
         z%sign = -1
         if (value == int(z'8000000000000000', int64)) then
            ! Avoid overflow in ABS for the most-negative int64.
            z = bigz_from_string('-9223372036854775808')
            return
         end if
         mag = -value
      else
         z%sign = 1
         mag = value
      end if

      n = 0
      do while (mag > 0_int64)
         n = n + 1
         mag = mag / limb_base
      end do
      allocate(z%limb(n))
      mag = merge(-value, value, value < 0_int64)
      do n = 1, size(z%limb)
         z%limb(n) = int(mod(mag, limb_base), int32)
         mag = mag / limb_base
      end do
   end function bigz_from_int64

   pure recursive function bigz_from_string(text, base) result(z)
      character(len=*), intent(in) :: text !! Character representation of an integer, optionally preceded by a sign.
      integer, intent(in), optional :: base !! Input radix from 2 through 36; defaults to decimal.
      type(bigz) :: z
      integer :: radix, i, first, digit, sgn
      character :: ch
      logical :: valid

      radix = 10
      if (present(base)) radix = base
      valid = radix >= 2 .and. radix <= 36
      sgn = 1
      first = 1
      do while (first <= len(text) .and. text(first:first) == ' ')
         first = first + 1
      end do
      if (first > len(text)) valid = .false.
      if (valid) then
         if (text(first:first) == '-') then
            sgn = -1
            first = first + 1
         else if (text(first:first) == '+') then
            first = first + 1
         end if
      end if
      if (first > len(text)) valid = .false.

      z = bigz_from_int64(0_int64)
      if (valid) then
         do i = first, len_trim(text)
            ch = text(i:i)
            if (ch >= '0' .and. ch <= '9') then
               digit = iachar(ch) - iachar('0')
            else if (ch >= 'A' .and. ch <= 'Z') then
               digit = 10 + iachar(ch) - iachar('A')
            else if (ch >= 'a' .and. ch <= 'z') then
               digit = 10 + iachar(ch) - iachar('a')
            else
               valid = .false.
               exit
            end if
            if (digit >= radix) then
               valid = .false.
               exit
            end if
            z = bigz_mul_small(z, int(radix, int64))
            z = bigz_add_small_nonnegative(z, int(digit, int64))
         end do
      end if
      if (valid .and. .not. bigz_is_zero(z)) z%sign = sgn
      if (.not. valid) z = bigz_from_int64(0_int64)
   end function bigz_from_string

   pure function bigz_to_string(z, base) result(text)
      type(bigz), intent(in) :: z !! Arbitrary-precision integer to format.
      integer, intent(in), optional :: base !! Output radix from 2 through 36; defaults to decimal.
      character(len=:), allocatable :: text
      type(bigz) :: work, q
      integer(int64) :: rem
      integer :: radix, n, i
      character(len=1), allocatable :: chars(:)
      character(len=36), parameter :: alphabet = '0123456789abcdefghijklmnopqrstuvwxyz'

      radix = 10
      if (present(base)) radix = base
      if (radix < 2 .or. radix > 36) then
         text = ''
         return
      end if
      if (bigz_is_zero(z)) then
         text = '0'
         return
      end if

      work = bigz_abs(z)
      n = max(1, bigz_sizeinbase(work, radix) + 1)
      allocate(chars(n))
      i = 0
      do while (.not. bigz_is_zero(work))
         call bigz_div_small(work, int(radix, int64), q, rem)
         i = i + 1
         chars(i) = alphabet(int(rem) + 1:int(rem) + 1)
         work = q
      end do
      if (z%sign < 0) then
         allocate(character(len=i + 1) :: text)
         text(1:1) = '-'
         do n = 1, i
            text(n + 1:n + 1) = chars(i - n + 1)
         end do
      else
         allocate(character(len=i) :: text)
         do n = 1, i
            text(n:n) = chars(i - n + 1)
         end do
      end if
   end function bigz_to_string

   pure elemental real(dp) function bigz_to_dp(z) result(value)
      type(bigz), intent(in) :: z !! Arbitrary-precision integer to approximate in the package real kind.
      integer :: i

      value = 0.0_dp
      if (.not. allocated(z%limb)) return
      do i = size(z%limb), 1, -1
         value = value * real(limb_base, dp) + real(z%limb(i), dp)
      end do
      value = real(z%sign, dp) * value
   end function bigz_to_dp

   pure elemental integer function bigz_compare(a, b) result(cmp)
      type(bigz), intent(in) :: a !! Left integer in the comparison.
      type(bigz), intent(in) :: b !! Right integer in the comparison.
      integer :: magcmp

      if (a%sign < b%sign) then
         cmp = -1
      else if (a%sign > b%sign) then
         cmp = 1
      else if (a%sign == 0) then
         cmp = 0
      else
         magcmp = compare_magnitude(a, b)
         cmp = a%sign * magcmp
      end if
   end function bigz_compare

   pure elemental logical function bigz_equal(a, b) result(equal)
      type(bigz), intent(in) :: a !! Left arbitrary-precision integer.
      type(bigz), intent(in) :: b !! Right arbitrary-precision integer.
      equal = bigz_compare(a, b) == 0
   end function bigz_equal

   pure elemental logical function bigz_is_zero(a) result(is_zero)
      type(bigz), intent(in) :: a !! Integer tested for exact equality to zero.
      is_zero = a%sign == 0
   end function bigz_is_zero

   pure elemental logical function bigz_is_one(a) result(is_one)
      type(bigz), intent(in) :: a !! Integer tested for exact equality to one.
      is_one = .false.
      if (a%sign /= 1) return
      if (.not. allocated(a%limb)) return
      if (size(a%limb) /= 1) return
      is_one = a%limb(1) == 1_int32
   end function bigz_is_one

   pure elemental logical function bigz_is_even(a) result(is_even)
      type(bigz), intent(in) :: a !! Integer tested for divisibility by two.
      if (a%sign == 0) then
         is_even = .true.
      else
         is_even = mod(a%limb(1), 2_int32) == 0_int32
      end if
   end function bigz_is_even

   pure elemental function bigz_add(a, b) result(c)
      type(bigz), intent(in) :: a !! Left addend.
      type(bigz), intent(in) :: b !! Right addend.
      type(bigz) :: c
      integer :: cmp

      if (a%sign == 0) then
         c = b
      else if (b%sign == 0) then
         c = a
      else if (a%sign == b%sign) then
         c = add_magnitude(a, b)
         c%sign = a%sign
      else
         cmp = compare_magnitude(a, b)
         if (cmp == 0) then
            c = bigz_from_int64(0_int64)
         else if (cmp > 0) then
            c = sub_magnitude(a, b)
            c%sign = a%sign
         else
            c = sub_magnitude(b, a)
            c%sign = b%sign
         end if
      end if
   end function bigz_add

   pure elemental function bigz_sub(a, b) result(c)
      type(bigz), intent(in) :: a !! Minuend.
      type(bigz), intent(in) :: b !! Subtrahend.
      type(bigz) :: c
      c = bigz_add(a, bigz_neg(b))
   end function bigz_sub

   pure elemental function bigz_neg(a) result(c)
      type(bigz), intent(in) :: a !! Integer whose sign is reversed.
      type(bigz) :: c
      c = a
      c%sign = -c%sign
   end function bigz_neg

   pure elemental function bigz_abs(a) result(c)
      type(bigz), intent(in) :: a !! Integer whose nonnegative magnitude is returned.
      type(bigz) :: c
      c = a
      if (c%sign < 0) c%sign = 1
   end function bigz_abs

   pure elemental function bigz_mul(a, b) result(c)
      type(bigz), intent(in) :: a !! Left multiplicand.
      type(bigz), intent(in) :: b !! Right multiplicand.
      type(bigz) :: c
      integer(int64), allocatable :: work(:)
      integer(int64) :: carry, value
      integer :: i, j, n

      if (a%sign == 0 .or. b%sign == 0) return
      n = size(a%limb) + size(b%limb)
      allocate(work(n))
      work = 0_int64
      do i = 1, size(a%limb)
         do j = 1, size(b%limb)
            work(i + j - 1) = work(i + j - 1) + int(a%limb(i), int64) * int(b%limb(j), int64)
         end do
      end do
      carry = 0_int64
      do i = 1, n
         value = work(i) + carry
         work(i) = mod(value, limb_base)
         carry = value / limb_base
      end do
      allocate(c%limb(n))
      c%limb = int(work, int32)
      c%sign = a%sign * b%sign
      call normalize_bigz(c)
   end function bigz_mul

   pure elemental function bigz_mul_small(a, factor) result(c)
      type(bigz), intent(in) :: a !! Arbitrary-precision multiplicand.
      integer(int64), intent(in) :: factor !! Signed machine-integer factor; magnitude must fit in int64.
      type(bigz) :: c
      integer(int64) :: fmag, carry, value
      integer :: i, n

      if (a%sign == 0 .or. factor == 0_int64) return
      if (factor == int(z'8000000000000000', int64)) then
         c = bigz_mul(a, bigz_from_int64(factor))
         return
      end if
      fmag = abs(factor)
      n = size(a%limb) + max(1, decimal_limbs_int64(fmag))
      allocate(c%limb(n))
      c%limb = 0_int32
      carry = 0_int64
      do i = 1, size(a%limb)
         value = int(a%limb(i), int64) * fmag + carry
         c%limb(i) = int(mod(value, limb_base), int32)
         carry = value / limb_base
      end do
      i = size(a%limb) + 1
      do while (carry > 0_int64)
         c%limb(i) = int(mod(carry, limb_base), int32)
         carry = carry / limb_base
         i = i + 1
      end do
      c%sign = a%sign * merge(-1, 1, factor < 0_int64)
      call normalize_bigz(c)
   end function bigz_mul_small

   pure subroutine bigz_divmod(a, b, quotient, remainder, info)
      type(bigz), intent(in) :: a !! Dividend.
      type(bigz), intent(in) :: b !! Nonzero divisor.
      type(bigz), intent(out) :: quotient !! Floor quotient, matching R integer-division semantics.
      type(bigz), intent(out) :: remainder !! Remainder satisfying a = quotient*b + remainder.
      integer, intent(out), optional :: info !! Zero on success and one when division by zero was requested.
      type(bigz) :: qa, ra, absa, absb
      integer :: stat

      stat = 0
      if (bigz_is_zero(b)) then
         quotient = bigz_from_int64(0_int64)
         remainder = bigz_from_int64(0_int64)
         stat = 1
            return
      end if
      if (bigz_is_zero(a)) then
         quotient = bigz_from_int64(0_int64)
         remainder = bigz_from_int64(0_int64)
            return
      end if

      absa = bigz_abs(a)
      absb = bigz_abs(b)
      call divmod_positive(absa, absb, qa, ra)
      if (a%sign == b%sign) then
         quotient = qa
         quotient%sign = merge(0, 1, bigz_is_zero(qa))
         remainder = ra
         if (.not. bigz_is_zero(remainder)) remainder%sign = b%sign
      else if (bigz_is_zero(ra)) then
         quotient = qa
         if (.not. bigz_is_zero(quotient)) quotient%sign = -1
         remainder = ra
      else
         quotient = bigz_add(qa, bigz_from_int64(1_int64))
         quotient%sign = -1
         remainder = bigz_sub(absb, ra)
         if (.not. bigz_is_zero(remainder)) remainder%sign = b%sign
      end if
      if (present(info)) info = stat
   end subroutine bigz_divmod

   pure elemental function bigz_quotient(a, b) result(q)
      type(bigz), intent(in) :: a !! Dividend.
      type(bigz), intent(in) :: b !! Nonzero divisor.
      type(bigz) :: q, r
      call bigz_divmod(a, b, q, r)
   end function bigz_quotient

   pure elemental function bigz_modulo(a, b) result(r)
      type(bigz), intent(in) :: a !! Dividend.
      type(bigz), intent(in) :: b !! Nonzero divisor.
      type(bigz) :: q, r
      call bigz_divmod(a, b, q, r)
   end function bigz_modulo

   pure subroutine bigz_div_small(a, divisor, quotient, remainder, info)
      type(bigz), intent(in) :: a !! Nonnegative arbitrary-precision dividend.
      integer(int64), intent(in) :: divisor !! Positive machine-integer divisor.
      type(bigz), intent(out) :: quotient !! Truncated nonnegative quotient.
      integer(int64), intent(out) :: remainder !! Nonnegative remainder less than divisor.
      integer, intent(out), optional :: info !! Zero on success and one when divisor is nonpositive.
      integer(int64) :: current
      integer :: i

      quotient = bigz_from_int64(0_int64)
      remainder = 0_int64
      if (divisor <= 0_int64) then
         if (present(info)) info = 1
         return
      end if
      if (a%sign == 0) then
         if (present(info)) info = 0
         return
      end if
      allocate(quotient%limb(size(a%limb)))
      quotient%limb = 0_int32
      do i = size(a%limb), 1, -1
         current = remainder * limb_base + int(a%limb(i), int64)
         quotient%limb(i) = int(current / divisor, int32)
         remainder = mod(current, divisor)
      end do
      quotient%sign = a%sign
      call normalize_bigz(quotient)
      if (present(info)) info = 0
   end subroutine bigz_div_small

   pure elemental function bigz_gcd(a, b) result(g)
      type(bigz), intent(in) :: a !! First integer in the greatest-common-divisor calculation.
      type(bigz), intent(in) :: b !! Second integer in the greatest-common-divisor calculation.
      type(bigz) :: g, x, y, r

      x = bigz_abs(a)
      y = bigz_abs(b)
      do while (.not. bigz_is_zero(y))
         r = bigz_modulo(x, y)
         x = y
         y = r
      end do
      g = x
   end function bigz_gcd

   pure elemental function bigz_lcm(a, b) result(l)
      type(bigz), intent(in) :: a !! First integer in the least-common-multiple calculation.
      type(bigz), intent(in) :: b !! Second integer in the least-common-multiple calculation.
      type(bigz) :: l, g, q

      if (bigz_is_zero(a) .or. bigz_is_zero(b)) return
      g = bigz_gcd(a, b)
      q = bigz_quotient(bigz_abs(a), g)
      l = bigz_mul(q, bigz_abs(b))
   end function bigz_lcm

   pure subroutine bigz_gcdex(a, b, g, s, t)
      type(bigz), intent(in) :: a !! First integer in the extended Euclidean algorithm.
      type(bigz), intent(in) :: b !! Second integer in the extended Euclidean algorithm.
      type(bigz), intent(out) :: g !! Nonnegative greatest common divisor.
      type(bigz), intent(out) :: s !! Bezout coefficient multiplying a.
      type(bigz), intent(out) :: t !! Bezout coefficient multiplying b.
      type(bigz) :: old_r, r, old_s, ss, old_t, tt, q, tmp

      old_r = a
      r = b
      old_s = bigz_from_int64(1_int64)
      ss = bigz_from_int64(0_int64)
      old_t = bigz_from_int64(0_int64)
      tt = bigz_from_int64(1_int64)
      do while (.not. bigz_is_zero(r))
         q = trunc_quotient(old_r, r)
         tmp = r
         r = bigz_sub(old_r, bigz_mul(q, r))
         old_r = tmp
         tmp = ss
         ss = bigz_sub(old_s, bigz_mul(q, ss))
         old_s = tmp
         tmp = tt
         tt = bigz_sub(old_t, bigz_mul(q, tt))
         old_t = tmp
      end do
      if (old_r%sign < 0) then
         old_r = bigz_neg(old_r)
         old_s = bigz_neg(old_s)
         old_t = bigz_neg(old_t)
      end if
      g = old_r
      s = old_s
      t = old_t
   end subroutine bigz_gcdex

   pure elemental function bigz_pow(base, exponent) result(value)
      type(bigz), intent(in) :: base !! Integer base.
      integer(int64), intent(in) :: exponent !! Nonnegative machine-integer exponent.
      type(bigz) :: value, factor
      integer(int64) :: e

      value = bigz_from_int64(1_int64)
      if (exponent < 0_int64) then
         value = bigz_from_int64(0_int64)
         return
      end if
      factor = base
      e = exponent
      do while (e > 0_int64)
         if (mod(e, 2_int64) == 1_int64) value = bigz_mul(value, factor)
         e = e / 2_int64
         if (e > 0_int64) factor = bigz_mul(factor, factor)
      end do
   end function bigz_pow

   pure elemental function bigz_powm(base, exponent, modulus) result(value)
      type(bigz), intent(in) :: base !! Integer base reduced modulo modulus during exponentiation.
      type(bigz), intent(in) :: exponent !! Nonnegative arbitrary-precision exponent.
      type(bigz), intent(in) :: modulus !! Nonzero modulus; its absolute value defines the residue ring.
      type(bigz) :: value, factor, e, q
      integer(int64) :: bit

      if (bigz_is_zero(modulus)) then
         value = bigz_from_int64(0_int64)
         return
      end if
      value = bigz_modulo(bigz_from_int64(1_int64), bigz_abs(modulus))
      if (exponent%sign < 0) then
         factor = bigz_invmod(base, bigz_abs(modulus))
         if (bigz_is_zero(factor) .and. .not. bigz_is_one(bigz_abs(modulus))) then
            value = bigz_from_int64(0_int64)
            return
         end if
         e = bigz_abs(exponent)
      else
         factor = bigz_modulo(base, bigz_abs(modulus))
         e = exponent
      end if
      do while (.not. bigz_is_zero(e))
         call bigz_div_small(e, 2_int64, q, bit)
         if (bit == 1_int64) value = bigz_modulo(bigz_mul(value, factor), bigz_abs(modulus))
         e = q
         if (.not. bigz_is_zero(e)) factor = bigz_modulo(bigz_mul(factor, factor), bigz_abs(modulus))
      end do
   end function bigz_powm

   pure elemental function bigz_invmod(a, modulus) result(inv)
      type(bigz), intent(in) :: a !! Integer whose multiplicative inverse is requested.
      type(bigz), intent(in) :: modulus !! Nonzero modulus.
      type(bigz) :: inv, g, s, t, m

      m = bigz_abs(modulus)
      if (bigz_is_zero(m)) then
         inv = bigz_from_int64(0_int64)
         return
      end if
      call bigz_gcdex(a, m, g, s, t)
      if (.not. bigz_is_one(g)) then
         inv = bigz_from_int64(0_int64)
         return
      end if
      inv = bigz_modulo(s, m)
   end function bigz_invmod

   pure elemental integer function bigz_sizeinbase(a, base) result(n_digits)
      type(bigz), intent(in) :: a !! Integer whose digit count is requested.
      integer, intent(in) :: base !! Radix from 2 through 36.
      type(bigz) :: work, q
      integer(int64) :: rem

      if (base < 2 .or. base > 36) then
         n_digits = 0
         return
      end if
      if (bigz_is_zero(a)) then
         n_digits = 1
         return
      end if
      work = bigz_abs(a)
      n_digits = 0
      do while (.not. bigz_is_zero(work))
         call bigz_div_small(work, int(base, int64), q, rem)
         work = q
         n_digits = n_digits + 1
      end do
   end function bigz_sizeinbase

   pure subroutine bigz_frexp(a, fraction, exponent)
      type(bigz), intent(in) :: a !! Integer to decompose as fraction times two raised to exponent.
      real(dp), intent(out) :: fraction !! Signed fraction in [0.5,1) for nonzero input, or zero for zero input.
      integer(int64), intent(out) :: exponent !! Binary exponent in the decomposition.
      type(bigz) :: work, q
      integer(int64) :: rem
      real(dp) :: scale
      integer :: k

      if (bigz_is_zero(a)) then
         fraction = 0.0_dp
         exponent = 0_int64
         return
      end if
      work = bigz_abs(a)
      exponent = 0_int64
      do while (bigz_sizeinbase(work, 2) > 52)
         call bigz_div_small(work, 2_int64, q, rem)
         work = q
         exponent = exponent + 1_int64
      end do
      scale = bigz_to_dp(work)
      k = 0
      do while (scale >= 1.0_dp)
         scale = scale * 0.5_dp
         k = k + 1
      end do
      fraction = real(a%sign, dp) * scale
      exponent = exponent + int(k, int64)
   end subroutine bigz_frexp

   pure elemental real(dp) function bigz_log(a) result(value)
      type(bigz), intent(in) :: a !! Positive integer whose natural logarithm is requested.
      real(dp) :: frac
      integer(int64) :: expo

      if (a%sign <= 0) then
         value = -huge(1.0_dp)
         return
      end if
      call bigz_frexp(a, frac, expo)
      value = log(frac) + real(expo, dp) * log(2.0_dp)
   end function bigz_log

   pure elemental real(dp) function bigz_log10(a) result(value)
      type(bigz), intent(in) :: a !! Positive integer whose base-10 logarithm is requested.
      value = bigz_log(a) / log(10.0_dp)
   end function bigz_log10

   pure elemental real(dp) function bigz_log2(a) result(value)
      type(bigz), intent(in) :: a !! Positive integer whose base-2 logarithm is requested.
      value = bigz_log(a) / log(2.0_dp)
   end function bigz_log2

   pure elemental logical function bigz_fits_int64(a) result(fits)
      type(bigz), intent(in) :: a !! Integer tested for exact representability in signed int64.
      character(len=:), allocatable :: s
      character(len=19), parameter :: pmax = '9223372036854775807'
      character(len=19), parameter :: nmax = '9223372036854775808'
      character(len=:), allocatable :: mag

      if (bigz_is_zero(a)) then
         fits = .true.
         return
      end if
      s = bigz_to_string(a)
      if (a%sign < 0) then
         mag = s(2:)
         fits = len(mag) < 19 .or. (len(mag) == 19 .and. mag <= nmax)
      else
         mag = s
         fits = len(mag) < 19 .or. (len(mag) == 19 .and. mag <= pmax)
      end if
   end function bigz_fits_int64

   pure elemental integer(int64) function bigz_to_int64(a) result(value)
      type(bigz), intent(in) :: a !! Integer to convert when exactly representable in signed int64.
      character(len=:), allocatable :: s
      integer :: stat

      if (.not. bigz_fits_int64(a)) then
         value = 0_int64
         return
      end if
      s = bigz_to_string(a)
      read(s, *, iostat=stat) value
   end function bigz_to_int64

   pure elemental integer function compare_magnitude(a, b) result(cmp)
      type(bigz), intent(in) :: a !! First nonzero magnitude to compare.
      type(bigz), intent(in) :: b !! Second nonzero magnitude to compare.
      integer :: i

      if (.not. allocated(a%limb)) then
         cmp = merge(0, -1, .not. allocated(b%limb))
         return
      end if
      if (.not. allocated(b%limb)) then
         cmp = 1
         return
      end if
      if (size(a%limb) < size(b%limb)) then
         cmp = -1
         return
      else if (size(a%limb) > size(b%limb)) then
         cmp = 1
         return
      end if
      do i = size(a%limb), 1, -1
         if (a%limb(i) < b%limb(i)) then
            cmp = -1
            return
         else if (a%limb(i) > b%limb(i)) then
            cmp = 1
            return
         end if
      end do
      cmp = 0
   end function compare_magnitude

   pure elemental function add_magnitude(a, b) result(c)
      type(bigz), intent(in) :: a !! First nonnegative magnitude.
      type(bigz), intent(in) :: b !! Second nonnegative magnitude.
      type(bigz) :: c
      integer(int64) :: carry, value
      integer :: i, n

      n = max(size(a%limb), size(b%limb))
      allocate(c%limb(n + 1))
      c%limb = 0_int32
      carry = 0_int64
      do i = 1, n
         value = carry
         if (i <= size(a%limb)) value = value + int(a%limb(i), int64)
         if (i <= size(b%limb)) value = value + int(b%limb(i), int64)
         c%limb(i) = int(mod(value, limb_base), int32)
         carry = value / limb_base
      end do
      if (carry > 0_int64) c%limb(n + 1) = int(carry, int32)
      c%sign = 1
      call normalize_bigz(c)
   end function add_magnitude

   pure elemental function sub_magnitude(a, b) result(c)
      type(bigz), intent(in) :: a !! Larger nonnegative magnitude.
      type(bigz), intent(in) :: b !! Smaller nonnegative magnitude.
      type(bigz) :: c
      integer(int64) :: borrow, value
      integer :: i

      allocate(c%limb(size(a%limb)))
      c%limb = 0_int32
      borrow = 0_int64
      do i = 1, size(a%limb)
         value = int(a%limb(i), int64) - borrow
         if (i <= size(b%limb)) value = value - int(b%limb(i), int64)
         if (value < 0_int64) then
            value = value + limb_base
            borrow = 1_int64
         else
            borrow = 0_int64
         end if
         c%limb(i) = int(value, int32)
      end do
      c%sign = 1
      call normalize_bigz(c)
   end function sub_magnitude

   pure elemental function bigz_add_small_nonnegative(a, value) result(c)
      type(bigz), intent(in) :: a !! Nonnegative arbitrary-precision addend.
      integer(int64), intent(in) :: value !! Nonnegative small addend less than the limb base.
      type(bigz) :: c
      integer(int64) :: carry, current
      integer :: i, n

      if (a%sign == 0) then
         c = bigz_from_int64(value)
         return
      end if
      n = size(a%limb)
      allocate(c%limb(n + 1))
      c%limb = 0_int32
      c%limb(1:n) = a%limb
      carry = value
      i = 1
      do while (carry > 0_int64 .and. i <= n)
         current = int(c%limb(i), int64) + carry
         c%limb(i) = int(mod(current, limb_base), int32)
         carry = current / limb_base
         i = i + 1
      end do
      if (carry > 0_int64) c%limb(n + 1) = int(carry, int32)
      c%sign = 1
      call normalize_bigz(c)
   end function bigz_add_small_nonnegative

   pure subroutine divmod_positive(a, b, q, r)
      type(bigz), intent(in) :: a !! Nonnegative dividend.
      type(bigz), intent(in) :: b !! Positive divisor.
      type(bigz), intent(out) :: q !! Nonnegative truncated quotient.
      type(bigz), intent(out) :: r !! Nonnegative remainder.
      integer :: i, lo, hi, mid, best
      type(bigz) :: prod

      if (compare_magnitude(a, b) < 0) then
         q = bigz_from_int64(0_int64)
         r = a
         return
      end if
      allocate(q%limb(size(a%limb)))
      q%limb = 0_int32
      q%sign = 1
      r = bigz_from_int64(0_int64)
      do i = size(a%limb), 1, -1
         r = shift_base_add(r, int(a%limb(i), int64))
         lo = 0
         hi = int(limb_base - 1_int64)
         best = 0
         do while (lo <= hi)
            mid = lo + (hi - lo) / 2
            prod = bigz_mul_small(b, int(mid, int64))
            if (bigz_compare(prod, r) <= 0) then
               best = mid
               lo = mid + 1
            else
               hi = mid - 1
            end if
         end do
         q%limb(i) = int(best, int32)
         if (best > 0) r = bigz_sub(r, bigz_mul_small(b, int(best, int64)))
      end do
      call normalize_bigz(q)
      call normalize_bigz(r)
   end subroutine divmod_positive

   pure elemental function shift_base_add(a, digit) result(c)
      type(bigz), intent(in) :: a !! Nonnegative value to multiply by the internal limb base.
      integer(int64), intent(in) :: digit !! New least-significant limb, in 0 through limb_base-1.
      type(bigz) :: c
      integer :: n

      if (a%sign == 0) then
         c = bigz_from_int64(digit)
         return
      end if
      n = size(a%limb)
      allocate(c%limb(n + 1))
      c%limb(1) = int(digit, int32)
      c%limb(2:n + 1) = a%limb
      c%sign = 1
      call normalize_bigz(c)
   end function shift_base_add

   pure elemental function trunc_quotient(a, b) result(q)
      type(bigz), intent(in) :: a !! Dividend for quotient truncated toward zero.
      type(bigz), intent(in) :: b !! Nonzero divisor for quotient truncated toward zero.
      type(bigz) :: q, r

      call divmod_positive(bigz_abs(a), bigz_abs(b), q, r)
      if (.not. bigz_is_zero(q)) q%sign = a%sign * b%sign
   end function trunc_quotient

   pure subroutine normalize_bigz(a)
      type(bigz), intent(inout) :: a !! Integer normalized by removing high zero limbs and canonicalizing zero.
      integer :: n
      integer(int32), allocatable :: tmp(:)

      if (.not. allocated(a%limb)) then
         a%sign = 0
         return
      end if
      n = size(a%limb)
      do while (n > 0)
         if (a%limb(n) /= 0_int32) exit
         n = n - 1
      end do
      if (n == 0) then
         deallocate(a%limb)
         a%sign = 0
      else if (n < size(a%limb)) then
         allocate(tmp(n))
         tmp = a%limb(1:n)
         call move_alloc(tmp, a%limb)
         if (a%sign == 0) a%sign = 1
      end if
   end subroutine normalize_bigz

   pure elemental integer function decimal_limbs_int64(value) result(n)
      integer(int64), intent(in) :: value !! Nonnegative int64 magnitude whose internal-limb count is requested.
      integer(int64) :: work

      if (value <= 0_int64) then
         n = 1
         return
      end if
      work = value
      n = 0
      do while (work > 0_int64)
         n = n + 1
         work = work / limb_base
      end do
   end function decimal_limbs_int64

end module gmp_bigz
