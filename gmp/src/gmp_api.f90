module gmp_api
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite
   use gmp_kinds, only : dp
   use gmp_bigz
   use gmp_bigq
   use gmp_number_theory
   use gmp_matrix
   use gmp_vector
   implicit none
   private

   public :: dp, bigz, bigq
   public :: bigz_from_int64, bigz_from_string, bigz_to_string, bigz_to_dp, bigz_to_int64, bigz_fits_int64
   public :: bigz_add, bigz_sub, bigz_neg, bigz_abs, bigz_mul, bigz_mul_small, bigz_quotient, bigz_modulo
   public :: bigz_compare, bigz_equal, bigz_is_zero, bigz_is_one, bigz_is_even
   public :: bigz_gcd, bigz_lcm, bigz_gcdex, bigz_pow, bigz_powm, bigz_invmod
   public :: bigz_sizeinbase, bigz_frexp, bigz_log, bigz_log10, bigz_log2
   public :: bigq_make, bigq_from_int64, bigq_from_string, bigq_to_string, bigq_to_dp
   public :: bigq_add, bigq_sub, bigq_neg, bigq_abs, bigq_mul, bigq_div, bigq_pow
   public :: bigq_compare, bigq_equal, bigq_is_zero, bigq_is_integer
   public :: bigq_floor, bigq_ceiling, bigq_trunc, bigq_round0, bigq_round_digits, bigq_numerator, bigq_denominator
   public :: factorial_z, choose_z, fibnum_z, fibnum2_z, lucnum_z, lucnum2_z
   public :: stirling1_z, stirling1_all_z, stirling2_z, stirling2_all_z, eulerian_z, eulerian_all_z
   public :: bernoulli_q, dbinom_q, isprime_z, nextprime_z, factorize_z, urand_bigz
   public :: bigz_matmul, bigq_matmul, bigz_crossprod, bigz_tcrossprod, bigq_crossprod, bigq_tcrossprod
   public :: bigz_transpose, bigq_transpose, bigq_solve, bigq_inverse, bigz_solve, bigz_inverse
   public :: bigz_sum, bigz_prod, bigz_cumsum, bigz_min, bigz_max, bigz_gamma
   public :: bigq_sum, bigq_prod, bigq_cumsum, bigq_min, bigq_max, bigq_mean
   public :: bigz_unique, bigz_duplicated, bigq_unique, bigq_duplicated, bigz_diff, bigq_diff
   public :: bigz_divide, bigq_from_bigz, bigz_from_bigq
   public :: bigz_lt, bigz_lte, bigz_gt, bigz_gte, bigz_neq
   public :: bigq_lt, bigq_lte, bigq_gt, bigq_gte, bigq_neq
   public :: bigz_signum, bigq_signum, bigz_round_digits
   public :: bigz_floor, bigz_trunc, bigz_is_whole, bigq_is_whole
   public :: bigz_is_finite, bigz_is_infinite, bigq_is_finite, bigq_is_infinite
   public :: real_is_whole
   public :: bigz_not, bigz_or, bigz_and, bigz_xor, bigq_not, bigq_or, bigq_and, bigq_xor
   public :: bigq_to_int64, bigq_set_numerator, bigq_set_denominator

contains

   pure elemental function bigz_divide(a, b) result(q)
      type(bigz), intent(in) :: a !! Integer numerator for exact rational division.
      type(bigz), intent(in) :: b !! Nonzero integer denominator for exact rational division.
      type(bigq) :: q
      q = bigq_make(a, b)
   end function bigz_divide

   pure elemental function bigq_from_bigz(a) result(q)
      type(bigz), intent(in) :: a !! Integer converted exactly to a rational with denominator one.
      type(bigq) :: q
      q = bigq_make(a, bigz_from_int64(1_int64))
   end function bigq_from_bigz

   pure elemental function bigz_from_bigq(a) result(z)
      type(bigq), intent(in) :: a !! Rational truncated toward zero for conversion to an integer.
      type(bigz) :: z
      z = bigq_trunc(a)
   end function bigz_from_bigq

   pure elemental logical function bigz_lt(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer comparison operand.
      type(bigz), intent(in) :: b !! Right integer comparison operand.
      value = bigz_compare(a, b) < 0
   end function bigz_lt

   pure elemental logical function bigz_lte(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer comparison operand.
      type(bigz), intent(in) :: b !! Right integer comparison operand.
      value = bigz_compare(a, b) <= 0
   end function bigz_lte

   pure elemental logical function bigz_gt(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer comparison operand.
      type(bigz), intent(in) :: b !! Right integer comparison operand.
      value = bigz_compare(a, b) > 0
   end function bigz_gt

   pure elemental logical function bigz_gte(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer comparison operand.
      type(bigz), intent(in) :: b !! Right integer comparison operand.
      value = bigz_compare(a, b) >= 0
   end function bigz_gte

   pure elemental logical function bigz_neq(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer comparison operand.
      type(bigz), intent(in) :: b !! Right integer comparison operand.
      value = .not. bigz_equal(a, b)
   end function bigz_neq

   pure elemental logical function bigq_lt(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational comparison operand.
      type(bigq), intent(in) :: b !! Right rational comparison operand.
      value = bigq_compare(a, b) < 0
   end function bigq_lt

   pure elemental logical function bigq_lte(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational comparison operand.
      type(bigq), intent(in) :: b !! Right rational comparison operand.
      value = bigq_compare(a, b) <= 0
   end function bigq_lte

   pure elemental logical function bigq_gt(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational comparison operand.
      type(bigq), intent(in) :: b !! Right rational comparison operand.
      value = bigq_compare(a, b) > 0
   end function bigq_gt

   pure elemental logical function bigq_gte(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational comparison operand.
      type(bigq), intent(in) :: b !! Right rational comparison operand.
      value = bigq_compare(a, b) >= 0
   end function bigq_gte

   pure elemental logical function bigq_neq(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational comparison operand.
      type(bigq), intent(in) :: b !! Right rational comparison operand.
      value = .not. bigq_equal(a, b)
   end function bigq_neq

   pure elemental integer function bigz_signum(a) result(value)
      type(bigz), intent(in) :: a !! Integer whose mathematical sign is requested.
      value = a%sign
   end function bigz_signum

   pure elemental integer function bigq_signum(a) result(value)
      type(bigq), intent(in) :: a !! Rational whose mathematical sign is requested.
      value = a%num%sign
   end function bigq_signum

   pure elemental function bigz_round_digits(a, digits) result(value)
      type(bigz), intent(in) :: a !! Integer value to round at a decimal position.
      integer, intent(in) :: digits !! Decimal digits; nonnegative values leave an integer unchanged.
      type(bigz) :: value, p10, rounded
      type(bigq) :: scaled
      integer :: i

      if (digits >= 0) then
         value = a
         return
      end if
      p10 = bigz_from_int64(1_int64)
      do i = 1, -digits
         p10 = bigz_mul(p10, bigz_from_int64(10_int64))
      end do
      scaled = bigq_make(a, p10)
      rounded = bigq_round0(scaled)
      value = bigz_mul(rounded, p10)
   end function bigz_round_digits

   pure elemental function bigz_floor(a) result(value)
      type(bigz), intent(in) :: a !! Integer whose mathematical floor is returned unchanged.
      type(bigz) :: value
      value = a
   end function bigz_floor

   pure elemental function bigz_trunc(a) result(value)
      type(bigz), intent(in) :: a !! Integer whose truncation toward zero is returned unchanged.
      type(bigz) :: value
      value = a
   end function bigz_trunc

   pure elemental logical function bigz_is_whole(a) result(value)
      type(bigz), intent(in) :: a !! Arbitrary-precision integer, always mathematically integral.
      value = bigz_equal(a, a)
   end function bigz_is_whole

   pure elemental logical function bigq_is_whole(a) result(value)
      type(bigq), intent(in) :: a !! Rational tested for an exact integer denominator after normalization.
      value = bigq_is_integer(a)
   end function bigq_is_whole

   pure elemental logical function bigz_is_finite(a) result(value)
      type(bigz), intent(in) :: a !! Representable arbitrary-precision integer, which is always finite.
      value = bigz_equal(a, a)
   end function bigz_is_finite

   pure elemental logical function bigz_is_infinite(a) result(value)
      type(bigz), intent(in) :: a !! Representable arbitrary-precision integer, which has no infinity encoding.
      value = .not. bigz_equal(a, a)
   end function bigz_is_infinite

   pure elemental logical function bigq_is_finite(a) result(value)
      type(bigq), intent(in) :: a !! Normalized rational value, which is always finite when representable.
      value = bigq_equal(a, a)
   end function bigq_is_finite

   pure elemental logical function bigq_is_infinite(a) result(value)
      type(bigq), intent(in) :: a !! Normalized rational value, which has no infinity encoding.
      value = .not. bigq_equal(a, a)
   end function bigq_is_infinite

   pure elemental logical function real_is_whole(x) result(value)
      real(dp), intent(in) :: x !! Real value tested for an exact integer value in the working real kind.
      if (ieee_is_nan(x)) then
         value = .false.
      else if (.not. ieee_is_finite(x)) then
         value = .false.
      else
         value = abs(x - anint(x)) <= 0.0_dp
      end if
   end function real_is_whole

   pure elemental logical function bigz_not(a) result(value)
      type(bigz), intent(in) :: a !! Integer interpreted as false only when exactly zero.
      value = bigz_is_zero(a)
   end function bigz_not

   pure elemental logical function bigz_or(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer truth-value operand.
      type(bigz), intent(in) :: b !! Right integer truth-value operand.
      value = .not. bigz_is_zero(a) .or. .not. bigz_is_zero(b)
   end function bigz_or

   pure elemental logical function bigz_and(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer truth-value operand.
      type(bigz), intent(in) :: b !! Right integer truth-value operand.
      value = .not. bigz_is_zero(a) .and. .not. bigz_is_zero(b)
   end function bigz_and

   pure elemental logical function bigz_xor(a, b) result(value)
      type(bigz), intent(in) :: a !! Left integer truth-value operand.
      type(bigz), intent(in) :: b !! Right integer truth-value operand.
      value = bigz_is_zero(a) .neqv. bigz_is_zero(b)
   end function bigz_xor

   pure elemental logical function bigq_not(a) result(value)
      type(bigq), intent(in) :: a !! Rational interpreted as false only when exactly zero.
      value = bigq_is_zero(a)
   end function bigq_not

   pure elemental logical function bigq_or(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational truth-value operand.
      type(bigq), intent(in) :: b !! Right rational truth-value operand.
      value = .not. bigq_is_zero(a) .or. .not. bigq_is_zero(b)
   end function bigq_or

   pure elemental logical function bigq_and(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational truth-value operand.
      type(bigq), intent(in) :: b !! Right rational truth-value operand.
      value = .not. bigq_is_zero(a) .and. .not. bigq_is_zero(b)
   end function bigq_and

   pure elemental logical function bigq_xor(a, b) result(value)
      type(bigq), intent(in) :: a !! Left rational truth-value operand.
      type(bigq), intent(in) :: b !! Right rational truth-value operand.
      value = bigq_is_zero(a) .neqv. bigq_is_zero(b)
   end function bigq_xor

   pure elemental integer(int64) function bigq_to_int64(a) result(value)
      type(bigq), intent(in) :: a !! Rational converted to a machine integer by truncation toward zero.
      type(bigz) :: z
      z = bigq_trunc(a)
      value = bigz_to_int64(z)
   end function bigq_to_int64

   pure elemental function bigq_set_numerator(a, numerator) result(value)
      type(bigq), intent(in) :: a !! Rational whose denominator is retained.
      type(bigz), intent(in) :: numerator !! Replacement exact numerator before normalization.
      type(bigq) :: value
      value = bigq_make(numerator, a%den)
   end function bigq_set_numerator

   pure elemental function bigq_set_denominator(a, denominator) result(value)
      type(bigq), intent(in) :: a !! Rational whose numerator is retained.
      type(bigz), intent(in) :: denominator !! Replacement nonzero denominator before normalization.
      type(bigq) :: value
      value = bigq_make(a%num, denominator)
   end function bigq_set_denominator

end module gmp_api
