module test_functions
  use rmpfr
  implicit none
  private

  public :: sine_fun, sqrt2_fun, quadratic_fun, sphere_fun

contains

  function sine_fun(arg) result(value)
    type(mpfr_real), intent(in) :: arg !! Integration argument for sin(x).
    type(mpfr_real) :: value

    value = mpfr_sin(arg)
  end function sine_fun

  function sqrt2_fun(arg) result(value)
    type(mpfr_real), intent(in) :: arg !! Root-finding argument for x**2 - 2.
    type(mpfr_real) :: value

    value = arg * arg - mpfr_from_integer(2_i64, mpfr_precision(arg))
  end function sqrt2_fun

  function quadratic_fun(arg) result(value)
    type(mpfr_real), intent(in) :: arg !! Optimization argument for a shifted quadratic.
    type(mpfr_real) :: value, center, two
    integer :: prec_bits

    prec_bits = mpfr_precision(arg)
    center = mpfr_from_string('1.2345', prec_bits)
    two = mpfr_from_integer(2_i64, prec_bits)
    value = (arg - center) * (arg - center) + two
  end function quadratic_fun

  function sphere_fun(arg) result(value)
    type(mpfr_real), intent(in) :: arg(:) !! Two-dimensional Hooke-Jeeves test point.
    type(mpfr_real) :: value, one, two
    integer :: prec_bits

    prec_bits = max(mpfr_precision(arg(1)), mpfr_precision(arg(2)))
    one = mpfr_from_integer(1_i64, prec_bits)
    two = mpfr_from_integer(2_i64, prec_bits)
    value = (arg(1) - one) * (arg(1) - one) + (arg(2) + two) * (arg(2) + two)
  end function sphere_fun

end module test_functions
