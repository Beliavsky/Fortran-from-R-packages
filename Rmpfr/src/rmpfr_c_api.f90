module rmpfr_c_api
  use iso_c_binding, only: c_char, c_double, c_int, c_long, c_ptr, c_size_t
  implicit none
  private

  type, bind(C), public :: mpfr_c_struct
    integer(c_long) :: prec
    integer(c_int) :: sign
    integer(c_long) :: exponent
    type(c_ptr) :: limbs
  end type mpfr_c_struct

  integer(c_int), parameter, public :: mpfr_rndn = 0_c_int
  integer(c_int), parameter, public :: mpfr_rndz = 1_c_int
  integer(c_int), parameter, public :: mpfr_rndu = 2_c_int
  integer(c_int), parameter, public :: mpfr_rndd = 3_c_int
  integer(c_int), parameter, public :: mpfr_rnda = 4_c_int
  public :: c_mpfr_init2, c_mpfr_clear, c_mpfr_set, c_mpfr_set_d, c_mpfr_set_si, c_mpfr_set_str, c_mpfr_get_d
  public :: c_mpfr_get_si, c_mpfr_get_prec, c_mpfr_prec_round, c_mpfr_set_nan, c_mpfr_set_inf, c_mpfr_set_zero
  public :: c_mpfr_nan_p, c_mpfr_inf_p, c_mpfr_number_p, c_mpfr_integer_p, c_mpfr_zero_p, c_mpfr_sgn, c_mpfr_cmp
  public :: c_mpfr_cmp_d, c_mpfr_cmp_si, c_mpfr_add, c_mpfr_sub, c_mpfr_mul, c_mpfr_div, c_mpfr_pow, c_mpfr_atan2
  public :: c_mpfr_hypot, c_mpfr_gamma_inc, c_mpfr_min, c_mpfr_max, c_mpfr_pow_si, c_mpfr_neg, c_mpfr_abs, c_mpfr_sqrt
  public :: c_mpfr_exp, c_mpfr_expm1, c_mpfr_log, c_mpfr_log1p, c_mpfr_log2, c_mpfr_log10, c_mpfr_sin, c_mpfr_cos
  public :: c_mpfr_tan, c_mpfr_asin, c_mpfr_acos, c_mpfr_atan, c_mpfr_sinh, c_mpfr_cosh, c_mpfr_tanh, c_mpfr_asinh
  public :: c_mpfr_acosh, c_mpfr_atanh, c_mpfr_gamma, c_mpfr_lngamma, c_mpfr_digamma, c_mpfr_erf, c_mpfr_erfc
  public :: c_mpfr_zeta, c_mpfr_eint, c_mpfr_li2, c_mpfr_j0, c_mpfr_j1, c_mpfr_y0, c_mpfr_y1, c_mpfr_ai, c_mpfr_floor
  public :: c_mpfr_ceil, c_mpfr_trunc, c_mpfr_fmod, c_mpfr_jn, c_mpfr_yn, c_mpfr_const_pi, c_mpfr_const_euler
  public :: c_mpfr_const_catalan, c_mpfr_const_log2, c_mpfr_fac_ui, c_mpfr_get_version, c_mpfr_get_str
  public :: c_mpfr_free_str, c_mpfr_set_default_prec, c_mpfr_get_default_prec, c_mpfr_get_exp, c_mpfr_frexp
  public :: c_mpfr_get_emin, c_mpfr_get_emax, c_mpfr_set_emin, c_mpfr_set_emax
  public :: c_mpfr_get_emin_min, c_mpfr_get_emin_max, c_mpfr_get_emax_min, c_mpfr_get_emax_max
  public :: c_mpfr_mul_2si

  interface
    subroutine c_mpfr_init2(x, prec) bind(C, name="mpfr_init2")
      import :: c_long, c_ptr
      type(c_ptr), value :: x !! Address of the MPFR value to initialize.
      integer(c_long), value :: prec !! Binary precision in bits; must satisfy MPFR's precision limits.
    end subroutine c_mpfr_init2

    subroutine c_mpfr_clear(x) bind(C, name="mpfr_clear")
      import :: c_ptr
      type(c_ptr), value :: x !! Address of an initialized MPFR value whose limb storage is released.
    end subroutine c_mpfr_clear

    integer(c_int) function c_mpfr_set(rop, op, rnd) bind(C, name="mpfr_set")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_set

    integer(c_int) function c_mpfr_set_d(rop, op, rnd) bind(C, name="mpfr_set_d")
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      real(c_double), value :: op !! IEEE double-precision value to convert.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_set_d

    integer(c_int) function c_mpfr_set_si(rop, op, rnd) bind(C, name="mpfr_set_si")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_long), value :: op !! Signed C long integer to convert.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_set_si

    integer(c_int) function c_mpfr_set_str(rop, text, base, rnd) bind(C, name="mpfr_set_str")
      import :: c_char, c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      character(kind=c_char), intent(in) :: text(*) !! Null-terminated numeric text accepted by MPFR.
      integer(c_int), value :: base !! Radix from 2 through 62, or zero for prefix-based autodetection.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_set_str

    real(c_double) function c_mpfr_get_d(op, rnd) bind(C, name="mpfr_get_d")
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding mode used for conversion to double precision.
    end function c_mpfr_get_d

    integer(c_long) function c_mpfr_get_si(op, rnd) bind(C, name="mpfr_get_si")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding mode used for integer conversion.
    end function c_mpfr_get_si

    pure integer(c_long) function c_mpfr_get_prec(op) bind(C, name="mpfr_get_prec")
      import :: c_long, c_ptr
      type(c_ptr), value :: op !! Address of the initialized MPFR value whose precision is queried.
    end function c_mpfr_get_prec

    integer(c_int) function c_mpfr_prec_round(op, prec, rnd) bind(C, name="mpfr_prec_round")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: op !! Address of the initialized MPFR value to round in place.
      integer(c_long), value :: prec !! New binary precision in bits.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_prec_round

    subroutine c_mpfr_set_nan(op) bind(C, name="mpfr_set_nan")
      import :: c_ptr
      type(c_ptr), value :: op !! Address of the initialized MPFR value to set to NaN.
    end subroutine c_mpfr_set_nan

    subroutine c_mpfr_set_inf(op, sign) bind(C, name="mpfr_set_inf")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the initialized MPFR value to set to infinity.
      integer(c_int), value :: sign !! Positive for +Inf and negative for -Inf.
    end subroutine c_mpfr_set_inf

    subroutine c_mpfr_set_zero(op, sign) bind(C, name="mpfr_set_zero")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the initialized MPFR value to set to signed zero.
      integer(c_int), value :: sign !! Positive for +0 and negative for -0.
    end subroutine c_mpfr_set_zero

    pure integer(c_int) function c_mpfr_nan_p(op) bind(C, name="mpfr_nan_p")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the MPFR value queried by mpfr_nan_p.
    end function c_mpfr_nan_p

    pure integer(c_int) function c_mpfr_inf_p(op) bind(C, name="mpfr_inf_p")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the MPFR value queried by mpfr_inf_p.
    end function c_mpfr_inf_p

    pure integer(c_int) function c_mpfr_number_p(op) bind(C, name="mpfr_number_p")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the MPFR value queried by mpfr_number_p.
    end function c_mpfr_number_p

    pure integer(c_int) function c_mpfr_integer_p(op) bind(C, name="mpfr_integer_p")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the MPFR value queried by mpfr_integer_p.
    end function c_mpfr_integer_p

    pure integer(c_int) function c_mpfr_zero_p(op) bind(C, name="mpfr_zero_p")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the MPFR value queried by mpfr_zero_p.
    end function c_mpfr_zero_p

    pure integer(c_int) function c_mpfr_sgn(op) bind(C, name="mpfr_sgn")
      import :: c_int, c_ptr
      type(c_ptr), value :: op !! Address of the MPFR value queried by mpfr_sgn.
    end function c_mpfr_sgn

    pure integer(c_int) function c_mpfr_cmp(a, b) bind(C, name="mpfr_cmp")
      import :: c_int, c_ptr
      type(c_ptr), value :: a !! Address of the left MPFR comparison operand.
      type(c_ptr), value :: b !! Address of the right MPFR comparison operand.
    end function c_mpfr_cmp

    pure integer(c_int) function c_mpfr_cmp_d(a, b) bind(C, name="mpfr_cmp_d")
      import :: c_double, c_int, c_ptr
      type(c_ptr), value :: a !! Address of the MPFR comparison operand.
      real(c_double), value :: b !! Double-precision comparison operand.
    end function c_mpfr_cmp_d

    pure integer(c_int) function c_mpfr_cmp_si(a, b) bind(C, name="mpfr_cmp_si")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: a !! Address of the MPFR comparison operand.
      integer(c_long), value :: b !! Signed C long comparison operand.
    end function c_mpfr_cmp_si

    integer(c_int) function c_mpfr_add(rop, a, b, rnd) bind(C, name="mpfr_add")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_add

    integer(c_int) function c_mpfr_sub(rop, a, b, rnd) bind(C, name="mpfr_sub")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_sub

    integer(c_int) function c_mpfr_mul(rop, a, b, rnd) bind(C, name="mpfr_mul")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_mul

    integer(c_int) function c_mpfr_div(rop, a, b, rnd) bind(C, name="mpfr_div")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_div

    integer(c_int) function c_mpfr_pow(rop, a, b, rnd) bind(C, name="mpfr_pow")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_pow

    integer(c_int) function c_mpfr_atan2(rop, a, b, rnd) bind(C, name="mpfr_atan2")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_atan2

    integer(c_int) function c_mpfr_hypot(rop, a, b, rnd) bind(C, name="mpfr_hypot")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_hypot

    integer(c_int) function c_mpfr_gamma_inc(rop, a, b, rnd) bind(C, name="mpfr_gamma_inc")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_gamma_inc

    integer(c_int) function c_mpfr_min(rop, a, b, rnd) bind(C, name="mpfr_min")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_min

    integer(c_int) function c_mpfr_max(rop, a, b, rnd) bind(C, name="mpfr_max")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the first MPFR operand.
      type(c_ptr), value :: b !! Address of the second MPFR operand.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_max

    integer(c_int) function c_mpfr_pow_si(rop, a, n, rnd) bind(C, name="mpfr_pow_si")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the MPFR base.
      integer(c_long), value :: n !! Signed integer exponent.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_pow_si

    integer(c_int) function c_mpfr_neg(rop, op, rnd) bind(C, name="mpfr_neg")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_neg

    integer(c_int) function c_mpfr_abs(rop, op, rnd) bind(C, name="mpfr_abs")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_abs

    integer(c_int) function c_mpfr_sqrt(rop, op, rnd) bind(C, name="mpfr_sqrt")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_sqrt

    integer(c_int) function c_mpfr_exp(rop, op, rnd) bind(C, name="mpfr_exp")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_exp

    integer(c_int) function c_mpfr_expm1(rop, op, rnd) bind(C, name="mpfr_expm1")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_expm1

    integer(c_int) function c_mpfr_log(rop, op, rnd) bind(C, name="mpfr_log")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_log

    integer(c_int) function c_mpfr_log1p(rop, op, rnd) bind(C, name="mpfr_log1p")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_log1p

    integer(c_int) function c_mpfr_log2(rop, op, rnd) bind(C, name="mpfr_log2")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_log2

    integer(c_int) function c_mpfr_log10(rop, op, rnd) bind(C, name="mpfr_log10")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_log10

    integer(c_int) function c_mpfr_sin(rop, op, rnd) bind(C, name="mpfr_sin")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_sin

    integer(c_int) function c_mpfr_cos(rop, op, rnd) bind(C, name="mpfr_cos")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_cos

    integer(c_int) function c_mpfr_tan(rop, op, rnd) bind(C, name="mpfr_tan")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_tan

    integer(c_int) function c_mpfr_asin(rop, op, rnd) bind(C, name="mpfr_asin")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_asin

    integer(c_int) function c_mpfr_acos(rop, op, rnd) bind(C, name="mpfr_acos")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_acos

    integer(c_int) function c_mpfr_atan(rop, op, rnd) bind(C, name="mpfr_atan")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_atan

    integer(c_int) function c_mpfr_sinh(rop, op, rnd) bind(C, name="mpfr_sinh")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_sinh

    integer(c_int) function c_mpfr_cosh(rop, op, rnd) bind(C, name="mpfr_cosh")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_cosh

    integer(c_int) function c_mpfr_tanh(rop, op, rnd) bind(C, name="mpfr_tanh")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_tanh

    integer(c_int) function c_mpfr_asinh(rop, op, rnd) bind(C, name="mpfr_asinh")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_asinh

    integer(c_int) function c_mpfr_acosh(rop, op, rnd) bind(C, name="mpfr_acosh")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_acosh

    integer(c_int) function c_mpfr_atanh(rop, op, rnd) bind(C, name="mpfr_atanh")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_atanh

    integer(c_int) function c_mpfr_gamma(rop, op, rnd) bind(C, name="mpfr_gamma")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_gamma

    integer(c_int) function c_mpfr_lngamma(rop, op, rnd) bind(C, name="mpfr_lngamma")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_lngamma

    integer(c_int) function c_mpfr_digamma(rop, op, rnd) bind(C, name="mpfr_digamma")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_digamma

    integer(c_int) function c_mpfr_erf(rop, op, rnd) bind(C, name="mpfr_erf")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_erf

    integer(c_int) function c_mpfr_erfc(rop, op, rnd) bind(C, name="mpfr_erfc")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_erfc

    integer(c_int) function c_mpfr_zeta(rop, op, rnd) bind(C, name="mpfr_zeta")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_zeta

    integer(c_int) function c_mpfr_eint(rop, op, rnd) bind(C, name="mpfr_eint")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_eint

    integer(c_int) function c_mpfr_li2(rop, op, rnd) bind(C, name="mpfr_li2")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_li2

    integer(c_int) function c_mpfr_j0(rop, op, rnd) bind(C, name="mpfr_j0")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_j0

    integer(c_int) function c_mpfr_j1(rop, op, rnd) bind(C, name="mpfr_j1")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_j1

    integer(c_int) function c_mpfr_y0(rop, op, rnd) bind(C, name="mpfr_y0")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_y0

    integer(c_int) function c_mpfr_y1(rop, op, rnd) bind(C, name="mpfr_y1")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_y1

    integer(c_int) function c_mpfr_ai(rop, op, rnd) bind(C, name="mpfr_ai")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_ai

    integer(c_int) function c_mpfr_floor(rop, op) bind(C, name="mpfr_floor")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
    end function c_mpfr_floor

    integer(c_int) function c_mpfr_ceil(rop, op) bind(C, name="mpfr_ceil")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
    end function c_mpfr_ceil

    integer(c_int) function c_mpfr_trunc(rop, op) bind(C, name="mpfr_trunc")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
    end function c_mpfr_trunc

    integer(c_int) function c_mpfr_fmod(rop, a, b, rnd) bind(C, name="mpfr_fmod")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: a !! Address of the dividend MPFR value.
      type(c_ptr), value :: b !! Address of the divisor MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_fmod

    integer(c_int) function c_mpfr_jn(rop, n, op, rnd) bind(C, name="mpfr_jn")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_long), value :: n !! Signed Bessel-function order.
      type(c_ptr), value :: op !! Address of the MPFR argument.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_jn

    integer(c_int) function c_mpfr_yn(rop, n, op, rnd) bind(C, name="mpfr_yn")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_long), value :: n !! Signed Bessel-function order.
      type(c_ptr), value :: op !! Address of the MPFR argument.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_yn

    integer(c_int) function c_mpfr_const_pi(rop, rnd) bind(C, name="mpfr_const_pi")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_const_pi

    integer(c_int) function c_mpfr_const_euler(rop, rnd) bind(C, name="mpfr_const_euler")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_const_euler

    integer(c_int) function c_mpfr_const_catalan(rop, rnd) bind(C, name="mpfr_const_catalan")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_const_catalan

    integer(c_int) function c_mpfr_const_log2(rop, rnd) bind(C, name="mpfr_const_log2")
      import :: c_int, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_const_log2

    integer(c_int) function c_mpfr_fac_ui(rop, n, rnd) bind(C, name="mpfr_fac_ui")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      integer(c_long), value :: n !! Nonnegative factorial argument passed in an unsigned-long-sized slot.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_fac_ui

    type(c_ptr) function c_mpfr_get_version() bind(C, name="mpfr_get_version")
      import :: c_ptr
    end function c_mpfr_get_version

    type(c_ptr) function c_mpfr_get_str(text, exponent, base, digits, op, rnd) bind(C, name="mpfr_get_str")
      import :: c_int, c_ptr, c_size_t
      type(c_ptr), value :: text !! Optional destination character buffer; pass C_NULL_PTR for MPFR allocation.
      type(c_ptr), value :: exponent !! Address receiving the base-dependent exponent.
      integer(c_int), value :: base !! Output radix from 2 through 62.
      integer(c_size_t), value :: digits !! Number of significand digits; zero requests enough digits for exact recovery.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_get_str

    subroutine c_mpfr_free_str(text) bind(C, name="mpfr_free_str")
      import :: c_ptr
      type(c_ptr), value :: text !! Pointer returned by mpfr_get_str when MPFR allocated the character buffer.
    end subroutine c_mpfr_free_str

    subroutine c_mpfr_set_default_prec(prec) bind(C, name="mpfr_set_default_prec")
      import :: c_long
      integer(c_long), value :: prec !! New MPFR default binary precision in bits.
    end subroutine c_mpfr_set_default_prec

    integer(c_long) function c_mpfr_get_default_prec() bind(C, name="mpfr_get_default_prec")
      import :: c_long
    end function c_mpfr_get_default_prec

    integer(c_long) function c_mpfr_get_exp(op) bind(C, name="mpfr_get_exp")
      import :: c_long, c_ptr
      type(c_ptr), value :: op !! Address of a finite nonzero MPFR value whose binary exponent is queried.
    end function c_mpfr_get_exp

    integer(c_int) function c_mpfr_frexp(exponent, fraction, op, rnd) bind(C, name="mpfr_frexp")
      import :: c_int, c_ptr
      type(c_ptr), value :: exponent !! Address receiving the extracted binary exponent.
      type(c_ptr), value :: fraction !! Address of the initialized MPFR significand result.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_frexp

    integer(c_int) function c_mpfr_mul_2si(rop, op, exponent, rnd) bind(C, name="mpfr_mul_2si")
      import :: c_int, c_long, c_ptr
      type(c_ptr), value :: rop !! Address of the initialized destination MPFR value.
      type(c_ptr), value :: op !! Address of the source MPFR value.
      integer(c_long), value :: exponent !! Signed power-of-two exponent.
      integer(c_int), value :: rnd !! MPFR rounding-mode code.
    end function c_mpfr_mul_2si

    integer(c_long) function c_mpfr_get_emin() bind(C, name="mpfr_get_emin")
      import :: c_long
    end function c_mpfr_get_emin

    integer(c_long) function c_mpfr_get_emax() bind(C, name="mpfr_get_emax")
      import :: c_long
    end function c_mpfr_get_emax

    integer(c_int) function c_mpfr_set_emin(exponent) bind(C, name="mpfr_set_emin")
      import :: c_int, c_long
      integer(c_long), value :: exponent !! New minimum MPFR exponent, subject to the library's supported range.
    end function c_mpfr_set_emin

    integer(c_int) function c_mpfr_set_emax(exponent) bind(C, name="mpfr_set_emax")
      import :: c_int, c_long
      integer(c_long), value :: exponent !! New maximum MPFR exponent, subject to the library's supported range.
    end function c_mpfr_set_emax

    integer(c_long) function c_mpfr_get_emin_min() bind(C, name="mpfr_get_emin_min")
      import :: c_long
    end function c_mpfr_get_emin_min

    integer(c_long) function c_mpfr_get_emin_max() bind(C, name="mpfr_get_emin_max")
      import :: c_long
    end function c_mpfr_get_emin_max

    integer(c_long) function c_mpfr_get_emax_min() bind(C, name="mpfr_get_emax_min")
      import :: c_long
    end function c_mpfr_get_emax_min

    integer(c_long) function c_mpfr_get_emax_max() bind(C, name="mpfr_get_emax_max")
      import :: c_long
    end function c_mpfr_get_emax_max

  end interface
end module rmpfr_c_api
