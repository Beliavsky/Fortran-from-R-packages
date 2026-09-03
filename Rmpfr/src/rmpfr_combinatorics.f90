module rmpfr_combinatorics
  use rmpfr_kinds, only: i64
  use rmpfr_types
  implicit none
  private

  public :: mpfr_choose, mpfr_pochhammer, mpfr_beta, mpfr_lbeta, mpfr_bernoulli

contains

  function mpfr_choose(x, n, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision upper argument of the generalized binomial coefficient.
    integer, intent(in) :: n !! Integer lower argument; negative values return zero as in upstream Rmpfr.
    integer, intent(in), optional :: rounding !! MPFR rounding code used for intermediate operations.
    type(mpfr_real) :: r, work, factor, denom, x_minus_n, two_n
    integer :: i, k, p

    p = mpfr_precision(x)
    if (n < 0) then
      r = mpfr_zero(1, p)
      return
    end if
    if (n == 0) then
      r = mpfr_from_integer(1_i64, p, rounding)
      return
    end if

    k = n
    if (mpfr_is_integer(x)) then
      x_minus_n = x - mpfr_from_integer(int(n, i64), p, rounding)
      two_n = mpfr_from_integer(int(2 * n, i64), p, rounding)
      if (x >= mpfr_from_integer(int(n, i64), p, rounding) .and. x < two_n) then
        if (mpfr_is_integer(x_minus_n)) k = int(mpfr_to_integer(x_minus_n))
      end if
    end if

    work = x
    r = x
    do i = 1, k - 1
      factor = x - mpfr_from_integer(int(i, i64), p, rounding)
      denom = mpfr_from_integer(int(i + 1, i64), p, rounding)
      r = (r * factor) / denom
      work = factor
    end do
  end function mpfr_choose

  function mpfr_pochhammer(x, n, rounding) result(r)
    type(mpfr_real), intent(in) :: x !! Arbitrary-precision starting value of the rising factorial.
    integer, intent(in) :: n !! Number of factors; zero returns one and negative values are unsupported.
    integer, intent(in), optional :: rounding !! MPFR rounding code used for intermediate operations.
    type(mpfr_real) :: r, factor
    integer :: i, p

    if (n < 0) error stop "Rmpfr: negative Pochhammer order is not supported by upstream kernel"
    p = mpfr_precision(x)
    if (n == 0) then
      r = mpfr_from_integer(1_i64, p, rounding)
      return
    end if
    r = x
    do i = 1, n - 1
      factor = x + mpfr_from_integer(int(i, i64), p, rounding)
      r = r * factor
    end do
  end function mpfr_pochhammer

  function mpfr_lbeta(a, b, rounding) result(r)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision beta-function argument.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision beta-function argument.
    integer, intent(in), optional :: rounding !! MPFR rounding code used for the calculation.
    type(mpfr_real) :: r, aa, bb, s, comb, denom
    integer(i64) :: bi
    integer :: p

    p = max(mpfr_precision(a), mpfr_precision(b))
    aa = mpfr_copy(a, p, rounding)
    bb = mpfr_copy(b, p, rounding)
    s = aa + bb
    if (mpfr_is_integer(s) .and. mpfr_sign(s) <= 0) then
      if (.not. mpfr_is_integer(aa) .and. .not. mpfr_is_integer(bb)) then
        r = mpfr_inf(-1, p)
        return
      end if
      if (mpfr_sign(aa) * mpfr_sign(bb) < 0) then
        if (mpfr_sign(bb) < 0) then
          s = aa
          aa = bb
          bb = s
          s = aa + bb
        end if
        if (.not. mpfr_is_integer(bb)) then
          r = mpfr_nan(p)
          return
        end if
        bi = mpfr_to_integer(bb)
        if (bi < 0_i64 .or. bi > int(huge(1), i64)) then
          r = mpfr_nan(p)
          return
        end if
        comb = mpfr_choose(s - mpfr_from_integer(1_i64, p, rounding), int(bi), rounding)
        denom = bb * comb
        r = -mpfr_log(mpfr_abs(denom, rounding), rounding)
        return
      end if
    end if
    r = mpfr_lgamma(aa, rounding) + mpfr_lgamma(bb, rounding) - mpfr_lgamma(s, rounding)
  end function mpfr_lbeta

  function mpfr_beta(a, b, rounding) result(r)
    type(mpfr_real), intent(in) :: a !! First arbitrary-precision beta-function argument.
    type(mpfr_real), intent(in) :: b !! Second arbitrary-precision beta-function argument.
    integer, intent(in), optional :: rounding !! MPFR rounding code used for the calculation.
    type(mpfr_real) :: r, aa, bb, s, comb, denom
    integer(i64) :: bi
    integer :: p, parity

    p = max(mpfr_precision(a), mpfr_precision(b))
    aa = mpfr_copy(a, p, rounding)
    bb = mpfr_copy(b, p, rounding)
    s = aa + bb
    if (mpfr_is_integer(s) .and. mpfr_sign(s) <= 0) then
      if (.not. mpfr_is_integer(aa) .and. .not. mpfr_is_integer(bb)) then
        r = mpfr_zero(1, p)
        return
      end if
      if (mpfr_sign(aa) * mpfr_sign(bb) < 0) then
        if (mpfr_sign(bb) < 0) then
          s = aa
          aa = bb
          bb = s
          s = aa + bb
        end if
        bi = mpfr_to_integer(bb)
        if (bi < 0_i64 .or. bi > int(huge(1), i64)) then
          r = mpfr_nan(p)
          return
        end if
        comb = mpfr_choose(s - mpfr_from_integer(1_i64, p, rounding), int(bi), rounding)
        denom = bb * comb
        r = mpfr_from_integer(1_i64, p, rounding) / denom
        parity = modulo(int(bi), 2)
        if (parity /= 0) r = -r
        return
      end if
    end if
    r = (mpfr_gamma(aa, rounding) * mpfr_gamma(bb, rounding)) / mpfr_gamma(s, rounding)
  end function mpfr_beta

  function mpfr_bernoulli(k, prec_bits, rounding) result(r)
    integer, intent(in) :: k !! Nonnegative Bernoulli-number index using the upstream Rmpfr B1=+1/2 convention.
    integer, intent(in), optional :: prec_bits !! Binary precision in bits; defaults to the current MPFR default.
    integer, intent(in), optional :: rounding !! MPFR rounding code for zeta and multiplication operations.
    type(mpfr_real) :: r, kk, argument
    integer :: p

    if (k < 0) error stop "Rmpfr: Bernoulli index must be nonnegative"
    p = mpfr_get_default_precision()
    if (present(prec_bits)) p = prec_bits
    if (k == 0) then
      r = mpfr_from_integer(1_i64, p, rounding)
      return
    end if
    kk = mpfr_from_integer(int(k, i64), p, rounding)
    argument = mpfr_from_integer(int(1 - k, i64), p, rounding)
    r = -kk * mpfr_zeta(argument, rounding)
  end function mpfr_bernoulli

end module rmpfr_combinatorics
