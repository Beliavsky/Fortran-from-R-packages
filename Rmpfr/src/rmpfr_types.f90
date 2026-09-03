module rmpfr_types
  use iso_c_binding, only: c_char, c_f_pointer, c_int, c_loc, c_long, c_null_char, c_size_t, c_null_ptr, c_ptr
  use rmpfr_kinds, only: dp, i64
  use rmpfr_c_api
  implicit none
  private

  integer, parameter, public :: rnd_nearest = 0
  integer, parameter, public :: rnd_zero = 1
  integer, parameter, public :: rnd_up = 2
  integer, parameter, public :: rnd_down = 3
  integer, parameter, public :: rnd_away = 4

  type, public :: mpfr_real
    private
    type(mpfr_c_struct), pointer :: raw => null()
  contains
    final :: finalize_mpfr_real
  end type mpfr_real

  interface assignment(=)
    module procedure assign_mpfr_real
  end interface assignment(=)

  interface operator(+)
    module procedure add_mpfr_real
  end interface operator(+)
  interface operator(-)
    module procedure sub_mpfr_real
    module procedure neg_mpfr_real
  end interface operator(-)
  interface operator(*)
    module procedure mul_mpfr_real
  end interface operator(*)
  interface operator(/)
    module procedure div_mpfr_real
  end interface operator(/)
  interface operator(**)
    module procedure pow_mpfr_real
    module procedure pow_mpfr_integer
  end interface operator(**)
  interface operator(==)
    module procedure eq_mpfr_real
  end interface operator(==)
  interface operator(/=)
    module procedure ne_mpfr_real
  end interface operator(/=)
  interface operator(<)
    module procedure lt_mpfr_real
  end interface operator(<)
  interface operator(<=)
    module procedure le_mpfr_real
  end interface operator(<=)
  interface operator(>)
    module procedure gt_mpfr_real
  end interface operator(>)
  interface operator(>=)
    module procedure ge_mpfr_real
  end interface operator(>=)

  public :: assignment(=), operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)
  public :: mpfr_from_real, mpfr_from_integer, mpfr_from_string
  public :: mpfr_nan, mpfr_inf, mpfr_zero, mpfr_copy, mpfr_round
  public :: mpfr_to_real, mpfr_to_integer, mpfr_to_string, mpfr_precision
  public :: mpfr_set_default_precision, mpfr_get_default_precision, mpfr_version
  public :: mpfr_get_exponent_range, mpfr_set_exponent_range, mpfr_get_exponent_limits
  public :: mpfr_is_nan, mpfr_is_infinite, mpfr_is_finite, mpfr_is_integer, mpfr_is_zero, mpfr_sign
  public :: mpfr_abs, mpfr_sqrt, mpfr_floor, mpfr_ceiling, mpfr_trunc
  public :: mpfr_exp, mpfr_expm1, mpfr_log, mpfr_log1p, mpfr_log2, mpfr_log10
  public :: mpfr_sin, mpfr_cos, mpfr_tan, mpfr_asin, mpfr_acos, mpfr_atan
  public :: mpfr_sinh, mpfr_cosh, mpfr_tanh, mpfr_asinh, mpfr_acosh, mpfr_atanh
  public :: mpfr_gamma, mpfr_lgamma, mpfr_digamma, mpfr_erf, mpfr_erfc
  public :: mpfr_zeta, mpfr_eint, mpfr_li2, mpfr_j0, mpfr_j1, mpfr_y0, mpfr_y1, mpfr_ai
  public :: mpfr_jn, mpfr_yn, mpfr_atan2, mpfr_hypot, mpfr_igamma
  public :: mpfr_const_pi, mpfr_const_euler, mpfr_const_catalan, mpfr_const_log2
  public :: mpfr_sinpi, mpfr_cospi, mpfr_tanpi, mpfr_ldexp, mpfr_frexp, mpfr_factorial

contains

  subroutine finalize_mpfr_real(x)
    type(mpfr_real), intent(inout) :: x !! Arbitrary-precision value whose MPFR storage is finalized.

    if (associated(x%raw)) then
      call c_mpfr_clear(c_loc(x%raw))
      deallocate(x%raw)
      nullify(x%raw)
    end if
  end subroutine finalize_mpfr_real

  impure elemental subroutine assign_mpfr_real(lhs, rhs)
    type(mpfr_real), intent(out) :: lhs !! Destination receiving an independent deep copy of the source value.
    type(mpfr_real), intent(in) :: rhs !! Initialized arbitrary-precision source value to copy.
    integer(c_int) :: status

    if (.not. associated(rhs%raw)) return
    call initialize_raw(lhs, mpfr_precision(rhs))
    status = c_mpfr_set(c_loc(lhs%raw), c_loc(rhs%raw), mpfr_rndn)
  end subroutine assign_mpfr_real

  subroutine initialize_raw(x, prec_bits)
    type(mpfr_real), intent(inout) :: x !! Arbitrary-precision value to initialize or reinitialize.
    integer, intent(in) :: prec_bits !! Requested binary precision in bits; values below two are invalid.

    if (prec_bits < 2) error stop "Rmpfr: precision must be at least 2 bits"
    if (associated(x%raw)) then
      call c_mpfr_clear(c_loc(x%raw))
      deallocate(x%raw)
      nullify(x%raw)
    end if
    allocate(x%raw)
    call c_mpfr_init2(c_loc(x%raw), int(prec_bits, c_long))
  end subroutine initialize_raw

  pure subroutine require_initialized(x)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value required to have initialized MPFR storage.

    if (.not. associated(x%raw)) error stop "Rmpfr: uninitialized mpfr_real"
  end subroutine require_initialized

  integer(c_int) function rounding_code(rounding) result(code)
    integer, intent(in), optional :: rounding !! Public rounding code: nearest, zero, up, down, or away.
    integer :: r

    r = rnd_nearest
    if (present(rounding)) r = rounding
    select case (r)
    case (rnd_nearest)
      code = mpfr_rndn
    case (rnd_zero)
      code = mpfr_rndz
    case (rnd_up)
      code = mpfr_rndu
    case (rnd_down)
      code = mpfr_rndd
    case (rnd_away)
      code = mpfr_rnda
    case default
      error stop "Rmpfr: invalid rounding mode"
    end select
  end function rounding_code

  integer function common_precision(a, b) result(prec_bits)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision operand.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision operand.

    prec_bits = max(mpfr_precision(a), mpfr_precision(b))
  end function common_precision

  function mpfr_from_real(value, prec_bits, rounding) result(x)
    real(dp), intent(in) :: value !! Double-precision value to convert to arbitrary precision.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: rounding !! Rounding code used during conversion.
    type(mpfr_real) :: x
    integer :: p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(x, p)
    status = c_mpfr_set_d(c_loc(x%raw), value, rounding_code(rounding))
  end function mpfr_from_real

  function mpfr_from_integer(value, prec_bits, rounding) result(x)
    integer(i64), intent(in) :: value !! Signed integer value to convert to arbitrary precision.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: rounding !! Rounding code used during conversion.
    type(mpfr_real) :: x
    integer :: p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(x, p)
    status = c_mpfr_set_si(c_loc(x%raw), int(value, c_long), rounding_code(rounding))
  end function mpfr_from_integer

  function mpfr_from_string(text, prec_bits, base, rounding) result(x)
    character(len=*), intent(in) :: text !! Numeric character representation accepted by MPFR.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: base !! Radix from 2 through 62, or zero for prefix-based autodetection.
    integer, intent(in), optional :: rounding !! Rounding code used during conversion.
    type(mpfr_real) :: x
    character(kind=c_char, len=:), allocatable :: ctext
    integer :: b, p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    b = 10
    if (present(base)) b = base
    call initialize_raw(x, p)
    ctext = trim(text) // c_null_char
    status = c_mpfr_set_str(c_loc(x%raw), ctext, int(b, c_int), rounding_code(rounding))
    if (status /= 0_c_int) error stop "Rmpfr: invalid numeric string"
  end function mpfr_from_string

  function mpfr_nan(prec_bits) result(x)
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    type(mpfr_real) :: x
    integer :: p

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(x, p)
    call c_mpfr_set_nan(c_loc(x%raw))
  end function mpfr_nan

  function mpfr_inf(sign, prec_bits) result(x)
    integer, intent(in), optional :: sign !! Sign of infinity; negative selects -Inf and other values select +Inf.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    type(mpfr_real) :: x
    integer :: p, s

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    s = 1
    if (present(sign)) then
      if (sign < 0) s = -1
    end if
    call initialize_raw(x, p)
    call c_mpfr_set_inf(c_loc(x%raw), int(s, c_int))
  end function mpfr_inf

  function mpfr_zero(sign, prec_bits) result(x)
    integer, intent(in), optional :: sign !! Sign of zero; negative selects -0 and other values select +0.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    type(mpfr_real) :: x
    integer :: p, s

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    s = 1
    if (present(sign)) then
      if (sign < 0) s = -1
    end if
    call initialize_raw(x, p)
    call c_mpfr_set_zero(c_loc(x%raw), int(s, c_int))
  end function mpfr_zero

  function mpfr_copy(source, prec_bits, rounding) result(x)
    type(mpfr_real), intent(in) :: source !! Arbitrary-precision value to copy.
    integer, intent(in), optional :: prec_bits !! Destination precision in bits; defaults to the source precision.
    integer, intent(in), optional :: rounding !! Rounding code used if the destination precision differs.
    type(mpfr_real) :: x
    integer :: p
    integer(c_int) :: status

    call require_initialized(source)
    p = mpfr_precision(source)
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(x, p)
    status = c_mpfr_set(c_loc(x%raw), c_loc(source%raw), rounding_code(rounding))
  end function mpfr_copy

  function mpfr_round(source, prec_bits, rounding) result(x)
    type(mpfr_real), intent(in) :: source !! Arbitrary-precision value to round to a new binary precision.
    integer, intent(in) :: prec_bits !! New binary precision in bits.
    integer, intent(in), optional :: rounding !! MPFR rounding code applied while reducing precision.
    type(mpfr_real) :: x
    integer(c_int) :: status

    x = source
    status = c_mpfr_prec_round(c_loc(x%raw), int(prec_bits, c_long), rounding_code(rounding))
  end function mpfr_round

  real(dp) function mpfr_to_real(x, rounding) result(value)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision source value to convert to real(dp).
    integer, intent(in), optional :: rounding !! MPFR rounding code used for conversion.

    call require_initialized(x)
    value = c_mpfr_get_d(c_loc(x%raw), rounding_code(rounding))
  end function mpfr_to_real

  integer(i64) function mpfr_to_integer(x, rounding) result(value)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision source value to convert to a signed integer.
    integer, intent(in), optional :: rounding !! MPFR rounding code used for conversion.

    call require_initialized(x)
    value = int(c_mpfr_get_si(c_loc(x%raw), rounding_code(rounding)), i64)
  end function mpfr_to_integer

  pure elemental integer function mpfr_precision(x) result(prec_bits)
    type(mpfr_real), intent(in) :: x !! Initialized arbitrary-precision value whose precision is queried.

    call require_initialized(x)
    prec_bits = int(c_mpfr_get_prec(c_loc(x%raw)))
  end function mpfr_precision

  subroutine mpfr_set_default_precision(prec_bits)
    integer, intent(in) :: prec_bits !! New MPFR default precision in binary bits.

    if (prec_bits < 2) error stop "Rmpfr: precision must be at least 2 bits"
    call c_mpfr_set_default_prec(int(prec_bits, c_long))
  end subroutine mpfr_set_default_precision

  integer function mpfr_get_default_precision() result(prec_bits)
    prec_bits = int(c_mpfr_get_default_prec())
  end function mpfr_get_default_precision

  function mpfr_version() result(version)
    character(len=:), allocatable :: version
    type(c_ptr) :: ptr
    character(kind=c_char), pointer :: chars(:)
    integer :: i, n

    ptr = c_mpfr_get_version()
    call c_f_pointer(ptr, chars, [128])
    n = 0
    do i = 1, size(chars)
      if (chars(i) == c_null_char) exit
      n = n + 1
    end do
    allocate(character(len=n) :: version)
    do i = 1, n
      version(i:i) = achar(iachar(chars(i)))
    end do
  end function mpfr_version

  function mpfr_to_string(x, digits, base, rounding) result(text)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value to format without first converting to hardware real.
    integer, intent(in), optional :: digits !! Significant digits to request; defaults to enough decimal digits for the precision.
    integer, intent(in), optional :: base !! Output radix; this translation currently formats base 10 canonically.
    integer, intent(in), optional :: rounding !! MPFR rounding code used while creating the significand digits.
    character(len=:), allocatable :: text
    type(c_ptr) :: ptr
    character(kind=c_char), pointer :: chars(:)
    character(len=:), allocatable :: mantissa
    integer(c_long), target :: exponent
    integer :: b, i, n, nd
    character(len=64) :: exp_text
    logical :: negative

    call require_initialized(x)
    if (mpfr_is_nan(x)) then
      text = "NaN"
      return
    end if
    if (mpfr_is_infinite(x)) then
      if (mpfr_sign(x) < 0) then
        text = "-Inf"
      else
        text = "Inf"
      end if
      return
    end if
    if (mpfr_is_zero(x)) then
      if (mpfr_sign(x) < 0) then
        text = "-0"
      else
        text = "0"
      end if
      return
    end if

    b = 10
    if (present(base)) b = base
    if (b /= 10) error stop "Rmpfr: mpfr_to_string currently supports base 10 only"
    nd = ceiling(real(mpfr_precision(x), dp) * log10(2.0_dp)) + 2
    if (present(digits)) nd = max(1, digits)
    exponent = 0_c_long
    ptr = c_mpfr_get_str(c_null_ptr, c_loc(exponent), int(b, c_int), int(nd, c_size_t), &
                         c_loc(x%raw), rounding_code(rounding))
    if (.not. c_associated_local(ptr)) error stop "Rmpfr: mpfr_get_str failed"
    call c_f_pointer(ptr, chars, [nd + 4])
    n = 0
    do i = 1, size(chars)
      if (chars(i) == c_null_char) exit
      n = n + 1
    end do
    negative = n > 0 .and. chars(1) == '-'
    allocate(character(len=n) :: mantissa)
    do i = 1, n
      mantissa(i:i) = achar(iachar(chars(i)))
    end do
    call c_mpfr_free_str(ptr)
    write(exp_text, '(i0)') exponent
    if (negative) then
      if (n == 1) then
        text = "NaN"
      else
        text = "-0." // mantissa(2:) // "e" // trim(exp_text)
      end if
    else
      text = "0." // mantissa // "e" // trim(exp_text)
    end if
  end function mpfr_to_string

  logical function c_associated_local(ptr) result(ok)
    use iso_c_binding, only: c_associated
    type(c_ptr), intent(in) :: ptr !! C pointer whose association status is tested.

    ok = c_associated(ptr)
  end function c_associated_local

  pure elemental logical function mpfr_is_nan(x) result(is_nan)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value tested for NaN.

    call require_initialized(x)
    is_nan = c_mpfr_nan_p(c_loc(x%raw)) /= 0_c_int
  end function mpfr_is_nan

  pure elemental logical function mpfr_is_infinite(x) result(is_infinite)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value tested for positive or negative infinity.

    call require_initialized(x)
    is_infinite = c_mpfr_inf_p(c_loc(x%raw)) /= 0_c_int
  end function mpfr_is_infinite

  pure elemental logical function mpfr_is_finite(x) result(is_finite)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value tested for a finite number.

    call require_initialized(x)
    is_finite = c_mpfr_number_p(c_loc(x%raw)) /= 0_c_int
  end function mpfr_is_finite

  pure elemental logical function mpfr_is_integer(x) result(is_integer)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value tested for an exact integral value.

    call require_initialized(x)
    is_integer = c_mpfr_integer_p(c_loc(x%raw)) /= 0_c_int
  end function mpfr_is_integer

  pure elemental logical function mpfr_is_zero(x) result(is_zero)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value tested for signed or unsigned zero.

    call require_initialized(x)
    is_zero = c_mpfr_zero_p(c_loc(x%raw)) /= 0_c_int
  end function mpfr_is_zero

  pure elemental integer function mpfr_sign(x) result(sign_value)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value whose sign is queried.

    call require_initialized(x)
    sign_value = int(c_mpfr_sgn(c_loc(x%raw)))
  end function mpfr_sign

  function add_mpfr_real(a, b) result(r)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision addend.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision addend.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_add(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), mpfr_rndn)
  end function add_mpfr_real

  function sub_mpfr_real(a, b) result(r)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision subtraction operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision subtraction operand.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_sub(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), mpfr_rndn)
  end function sub_mpfr_real

  function neg_mpfr_real(a) result(r)
    type(mpfr_real), intent(in) :: a !! Arbitrary-precision value to negate.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(a))
    status = c_mpfr_neg(c_loc(r%raw), c_loc(a%raw), mpfr_rndn)
  end function neg_mpfr_real

  function mul_mpfr_real(a, b) result(r)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision multiplication operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision multiplication operand.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_mul(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), mpfr_rndn)
  end function mul_mpfr_real

  function div_mpfr_real(a, b) result(r)
    type(mpfr_real), intent(in) :: a !! Arbitrary-precision numerator.
    type(mpfr_real), intent(in) :: b !! Arbitrary-precision denominator.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_div(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), mpfr_rndn)
  end function div_mpfr_real

  function pow_mpfr_real(a, b) result(r)
    type(mpfr_real), intent(in) :: a !! Arbitrary-precision base.
    type(mpfr_real), intent(in) :: b !! Arbitrary-precision exponent.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_pow(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), mpfr_rndn)
  end function pow_mpfr_real

  function pow_mpfr_integer(a, n) result(r)
    type(mpfr_real), intent(in) :: a !! Arbitrary-precision base.
    integer, intent(in) :: n !! Signed integer exponent.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(a))
    status = c_mpfr_pow_si(c_loc(r%raw), c_loc(a%raw), int(n, c_long), mpfr_rndn)
  end function pow_mpfr_integer

  pure elemental logical function eq_mpfr_real(a, b) result(value)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision comparison operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision comparison operand.

    value = c_mpfr_cmp(c_loc(a%raw), c_loc(b%raw)) == 0_c_int
  end function eq_mpfr_real

  pure elemental logical function ne_mpfr_real(a, b) result(value)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision comparison operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision comparison operand.

    value = c_mpfr_cmp(c_loc(a%raw), c_loc(b%raw)) /= 0_c_int
  end function ne_mpfr_real

  pure elemental logical function lt_mpfr_real(a, b) result(value)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision comparison operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision comparison operand.

    value = c_mpfr_cmp(c_loc(a%raw), c_loc(b%raw)) < 0_c_int
  end function lt_mpfr_real

  pure elemental logical function le_mpfr_real(a, b) result(value)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision comparison operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision comparison operand.

    value = c_mpfr_cmp(c_loc(a%raw), c_loc(b%raw)) <= 0_c_int
  end function le_mpfr_real

  pure elemental logical function gt_mpfr_real(a, b) result(value)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision comparison operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision comparison operand.

    value = c_mpfr_cmp(c_loc(a%raw), c_loc(b%raw)) > 0_c_int
  end function gt_mpfr_real

  pure elemental logical function ge_mpfr_real(a, b) result(value)
    type(mpfr_real), intent(in) :: a !! Left arbitrary-precision comparison operand.
    type(mpfr_real), intent(in) :: b !! Right arbitrary-precision comparison operand.

    value = c_mpfr_cmp(c_loc(a%raw), c_loc(b%raw)) >= 0_c_int
  end function ge_mpfr_real

  function mpfr_abs(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument whose absolute value is required.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_abs(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_abs

  function mpfr_floor(x) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value to round toward negative infinity.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_floor(c_loc(r%raw), c_loc(x%raw))
  end function mpfr_floor

  function mpfr_ceiling(x) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value to round toward positive infinity.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_ceil(c_loc(r%raw), c_loc(x%raw))
  end function mpfr_ceiling

  function mpfr_trunc(x) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value to round toward zero to an integer.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_trunc(c_loc(r%raw), c_loc(x%raw))
  end function mpfr_trunc

  function mpfr_sqrt(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's sqrt function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_sqrt(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_sqrt

  function mpfr_exp(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's exp function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_exp(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_exp

  function mpfr_expm1(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's expm1 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_expm1(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_expm1

  function mpfr_log(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's log function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_log(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_log

  function mpfr_log1p(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's log1p function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_log1p(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_log1p

  function mpfr_log2(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's log2 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_log2(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_log2

  function mpfr_log10(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's log10 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_log10(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_log10

  function mpfr_sin(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's sin function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_sin(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_sin

  function mpfr_cos(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's cos function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_cos(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_cos

  function mpfr_tan(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's tan function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_tan(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_tan

  function mpfr_asin(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's asin function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_asin(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_asin

  function mpfr_acos(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's acos function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_acos(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_acos

  function mpfr_atan(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's atan function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_atan(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_atan

  function mpfr_sinh(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's sinh function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_sinh(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_sinh

  function mpfr_cosh(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's cosh function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_cosh(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_cosh

  function mpfr_tanh(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's tanh function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_tanh(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_tanh

  function mpfr_asinh(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's asinh function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_asinh(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_asinh

  function mpfr_acosh(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's acosh function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_acosh(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_acosh

  function mpfr_atanh(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's atanh function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_atanh(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_atanh

  function mpfr_gamma(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's gamma function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_gamma(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_gamma

  function mpfr_lgamma(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's lngamma function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_lngamma(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_lgamma

  function mpfr_digamma(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's digamma function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_digamma(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_digamma

  function mpfr_erf(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's erf function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_erf(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_erf

  function mpfr_erfc(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's erfc function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_erfc(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_erfc

  function mpfr_zeta(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's zeta function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_zeta(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_zeta

  function mpfr_eint(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's eint function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_eint(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_eint

  function mpfr_li2(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's li2 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_li2(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_li2

  function mpfr_j0(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's j0 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_j0(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_j0

  function mpfr_j1(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's j1 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_j1(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_j1

  function mpfr_y0(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's y0 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_y0(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_y0

  function mpfr_y1(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's y1 function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_y1(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_y1

  function mpfr_ai(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to MPFR's ai function.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_ai(c_loc(r%raw), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_ai

  function mpfr_jn(n, x, rounding) result(r)
    integer, intent(in) :: n !! Signed integer order of the Bessel J function.
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision Bessel-function argument.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_jn(c_loc(r%raw), int(n, c_long), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_jn

  function mpfr_yn(n, x, rounding) result(r)
    integer, intent(in) :: n !! Signed integer order of the Bessel Y function.
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision Bessel-function argument.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_yn(c_loc(r%raw), int(n, c_long), c_loc(x%raw), rounding_code(rounding))
  end function mpfr_yn

  function mpfr_atan2(a, b, rounding) result(r)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision argument to atan2.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision argument to atan2.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_atan2(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), rounding_code(rounding))
  end function mpfr_atan2

  function mpfr_hypot(a, b, rounding) result(r)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision argument to hypot.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision argument to hypot.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_hypot(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), rounding_code(rounding))
  end function mpfr_hypot

  function mpfr_igamma(a, b, rounding) result(r)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision argument to igamma.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision argument to igamma.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, common_precision(a, b))
    status = c_mpfr_gamma_inc(c_loc(r%raw), c_loc(a%raw), c_loc(b%raw), rounding_code(rounding))
  end function mpfr_igamma

  function mpfr_const_pi(prec_bits, rounding) result(r)
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the constant.
    type(mpfr_real) :: r
    integer :: p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(r, p)
    status = c_mpfr_const_pi(c_loc(r%raw), rounding_code(rounding))
  end function mpfr_const_pi

  function mpfr_const_euler(prec_bits, rounding) result(r)
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the constant.
    type(mpfr_real) :: r
    integer :: p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(r, p)
    status = c_mpfr_const_euler(c_loc(r%raw), rounding_code(rounding))
  end function mpfr_const_euler

  function mpfr_const_catalan(prec_bits, rounding) result(r)
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the constant.
    type(mpfr_real) :: r
    integer :: p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(r, p)
    status = c_mpfr_const_catalan(c_loc(r%raw), rounding_code(rounding))
  end function mpfr_const_catalan

  function mpfr_const_log2(prec_bits, rounding) result(r)
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to MPFR's current default.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the constant.
    type(mpfr_real) :: r
    integer :: p
    integer(c_int) :: status

    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(r, p)
    status = c_mpfr_const_log2(c_loc(r%raw), rounding_code(rounding))
  end function mpfr_const_log2

  function mpfr_sinpi(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to sin(pi*x).
    integer, intent(in), optional :: rounding !! MPFR rounding code for intermediate and final operations.
    type(mpfr_real) :: r, two, tmp
    integer(c_int) :: status

    two = mpfr_from_integer(2_i64, mpfr_precision(x), rounding)
    tmp = mpfr_copy(x)
    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_fmod(c_loc(r%raw), c_loc(tmp%raw), c_loc(two%raw), rounding_code(rounding))
    if (c_mpfr_cmp_si(c_loc(r%raw), -1_c_long) <= 0_c_int) then
      r = r + two
    else if (c_mpfr_cmp_si(c_loc(r%raw), 1_c_long) > 0_c_int) then
      r = r - two
    end if
    if (mpfr_is_integer(r)) then
      r = mpfr_zero(1, mpfr_precision(x))
    else if (c_mpfr_cmp_d(c_loc(r%raw), 0.5_dp) == 0_c_int) then
      r = mpfr_from_integer(1_i64, mpfr_precision(x), rounding)
    else if (c_mpfr_cmp_d(c_loc(r%raw), -0.5_dp) == 0_c_int) then
      r = mpfr_from_integer(-1_i64, mpfr_precision(x), rounding)
    else
      tmp = mpfr_const_pi(mpfr_precision(x), rounding)
      r = mpfr_sin(r * tmp, rounding)
    end if
  end function mpfr_sinpi

  function mpfr_cospi(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to cos(pi*x).
    integer, intent(in), optional :: rounding !! MPFR rounding code for intermediate and final operations.
    type(mpfr_real) :: r, two, tmp
    integer(c_int) :: status

    two = mpfr_from_integer(2_i64, mpfr_precision(x), rounding)
    tmp = mpfr_abs(x, rounding)
    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_fmod(c_loc(r%raw), c_loc(tmp%raw), c_loc(two%raw), rounding_code(rounding))
    if (c_mpfr_cmp_d(c_loc(r%raw), 0.5_dp) == 0_c_int .or. &
        c_mpfr_cmp_d(c_loc(r%raw), 1.5_dp) == 0_c_int) then
      r = mpfr_zero(1, mpfr_precision(x))
    else if (c_mpfr_cmp_si(c_loc(r%raw), 1_c_long) == 0_c_int) then
      r = mpfr_from_integer(-1_i64, mpfr_precision(x), rounding)
    else if (c_mpfr_cmp_si(c_loc(r%raw), 0_c_long) == 0_c_int) then
      r = mpfr_from_integer(1_i64, mpfr_precision(x), rounding)
    else
      tmp = mpfr_const_pi(mpfr_precision(x), rounding)
      r = mpfr_cos(r * tmp, rounding)
    end if
  end function mpfr_cospi

  function mpfr_tanpi(x, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision argument to tan(pi*x).
    integer, intent(in), optional :: rounding !! MPFR rounding code for intermediate and final operations.
    type(mpfr_real) :: r, one, tmp
    integer(c_int) :: status

    one = mpfr_from_integer(1_i64, mpfr_precision(x), rounding)
    call initialize_raw(r, mpfr_precision(x))
    status = c_mpfr_fmod(c_loc(r%raw), c_loc(x%raw), c_loc(one%raw), rounding_code(rounding))
    if (c_mpfr_cmp_d(c_loc(r%raw), -0.5_dp) <= 0_c_int) then
      r = r + one
    else if (c_mpfr_cmp_d(c_loc(r%raw), 0.5_dp) > 0_c_int) then
      r = r - one
    end if
    if (mpfr_is_zero(r)) then
      r = mpfr_zero(1, mpfr_precision(x))
    else if (c_mpfr_cmp_d(c_loc(r%raw), 0.5_dp) == 0_c_int .or. &
             c_mpfr_cmp_d(c_loc(r%raw), -0.5_dp) == 0_c_int) then
      r = mpfr_inf(mpfr_sign(r), mpfr_precision(x))
    else
      tmp = mpfr_const_pi(mpfr_precision(x), rounding)
      r = mpfr_tan(r * tmp, rounding)
    end if
  end function mpfr_tanpi

  function mpfr_ldexp(fraction, exponent, rounding) result(r)
    type(mpfr_real), intent(in) :: fraction !! Arbitrary-precision significand to multiply by 2**exponent.
    integer, intent(in) :: exponent !! Signed binary exponent.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer(c_int) :: status

    call initialize_raw(r, mpfr_precision(fraction))
    status = c_mpfr_mul_2si(c_loc(r%raw), c_loc(fraction%raw), int(exponent, c_long), rounding_code(rounding))
  end function mpfr_ldexp

  subroutine mpfr_frexp(x, fraction, exponent, rounding)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision value to decompose as fraction*2**exponent.
    type(mpfr_real), intent(out) :: fraction !! Returned arbitrary-precision significand.
    integer, intent(out) :: exponent !! Returned signed binary exponent.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the significand.
    integer(c_long), target :: c_exponent
    integer(c_int) :: status

    call initialize_raw(fraction, mpfr_precision(x))
    c_exponent = 0_c_long
    status = c_mpfr_frexp(c_loc(c_exponent), c_loc(fraction%raw), c_loc(x%raw), rounding_code(rounding))
    exponent = int(c_exponent)
  end subroutine mpfr_frexp

  function mpfr_factorial(n, prec_bits, rounding) result(r)
    integer, intent(in) :: n !! Nonnegative integer whose factorial is required.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to enough bits for exact n! or MPFR default.
    integer, intent(in), optional :: rounding !! MPFR rounding code for the result.
    type(mpfr_real) :: r
    integer :: p
    integer(c_int) :: status

    if (n < 0) error stop "Rmpfr: factorial requires n >= 0"
    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    call initialize_raw(r, p)
    status = c_mpfr_fac_ui(c_loc(r%raw), int(n, c_long), rounding_code(rounding))
  end function mpfr_factorial

  subroutine mpfr_get_exponent_range(emin, emax)
    integer(i64), intent(out) :: emin !! Current minimum MPFR exponent used for underflow detection.
    integer(i64), intent(out) :: emax !! Current maximum MPFR exponent used for overflow detection.

    emin = int(c_mpfr_get_emin(), i64)
    emax = int(c_mpfr_get_emax(), i64)
  end subroutine mpfr_get_exponent_range

  subroutine mpfr_set_exponent_range(emin, emax)
    integer(i64), intent(in) :: emin !! Requested minimum MPFR exponent within the library-supported range.
    integer(i64), intent(in) :: emax !! Requested maximum MPFR exponent within the library-supported range.
    integer(c_int) :: status

    if (emin >= emax) error stop "Rmpfr: exponent range requires emin < emax"
    status = c_mpfr_set_emin(int(emin, c_long))
    if (status /= 0_c_int) error stop "Rmpfr: invalid minimum MPFR exponent"
    status = c_mpfr_set_emax(int(emax, c_long))
    if (status /= 0_c_int) error stop "Rmpfr: invalid maximum MPFR exponent"
  end subroutine mpfr_set_exponent_range

  subroutine mpfr_get_exponent_limits(emin_min, emin_max, emax_min, emax_max)
    integer(i64), intent(out) :: emin_min !! Smallest minimum-exponent setting supported by this MPFR build.
    integer(i64), intent(out) :: emin_max !! Largest minimum-exponent setting supported by this MPFR build.
    integer(i64), intent(out) :: emax_min !! Smallest maximum-exponent setting supported by this MPFR build.
    integer(i64), intent(out) :: emax_max !! Largest maximum-exponent setting supported by this MPFR build.

    emin_min = int(c_mpfr_get_emin_min(), i64)
    emin_max = int(c_mpfr_get_emin_max(), i64)
    emax_min = int(c_mpfr_get_emax_min(), i64)
    emax_max = int(c_mpfr_get_emax_max(), i64)
  end subroutine mpfr_get_exponent_limits

end module rmpfr_types
