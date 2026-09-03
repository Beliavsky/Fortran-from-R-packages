program test_rmpfr
  use rmpfr
  use test_functions
  implicit none

  integer, parameter :: p = 256
  integer :: exponent
  integer(i64) :: emin, emax, emin_min, emin_max, emax_min, emax_max
  type(mpfr_real) :: x, fpart
  type(mpfr_real), allocatable :: values(:), a(:, :), b(:, :), c(:, :)
  type(mpfr_integral_result) :: integ
  type(mpfr_optimize_result) :: opt
  type(mpfr_root_result) :: root, quant
  type(mpfr_hjk_result) :: hj

  call test_construction_and_math()
  call test_combinatorics()
  call test_probability()
  call test_summary_matrix()
  call test_algorithms()
  print '(a)', 'All Rmpfr deterministic tests passed.'

contains

  subroutine test_construction_and_math()
    type(mpfr_real) :: down, up, pi_value, half

    x = mpfr_from_string('3.1415926535897932384626433832795028841971693993751058209749445923', p)
    call assert_close(x, &
      '3.1415926535897932384626433832795028841971693993751058209749445923', '1e-70', 'string conversion')
    if (mpfr_precision(x) /= p) error stop 'precision preservation failed'
    if (len_trim(mpfr_version()) == 0) error stop 'MPFR version string is empty'

    down = mpfr_from_string('0.1', 24, rounding=rnd_down)
    up = mpfr_from_string('0.1', 24, rounding=rnd_up)
    if (.not. down < up) error stop 'directed string rounding failed'

    pi_value = mpfr_const_pi(p)
    call assert_close(pi_value, &
      '3.141592653589793238462643383279502884197169399375105820974944592307816406286', &
      '1e-75', 'pi constant')
    half = mpfr_from_string('0.5', p)
    call assert_close(mpfr_sinpi(half), '1', '1e-75', 'sinpi half')
    call assert_close(mpfr_cospi(mpfr_from_integer(1_i64, p)), '-1', '1e-75', 'cospi one')
    if (.not. mpfr_is_infinite(mpfr_tanpi(half))) error stop 'tanpi half should be infinite'

    call assert_close(mpfr_gamma(mpfr_from_string('0.5', p)), &
      '1.772453850905516027298167483341145182797549456122387128213807789852911284591', &
      '1e-74', 'gamma half')
    call assert_close(mpfr_zeta(mpfr_from_integer(2_i64, p)), &
      '1.644934066848226436472415166646025189218949901206798437735558229370007470404', &
      '1e-74', 'zeta two')

    call mpfr_frexp(pi_value, fpart, exponent)
    call assert_close(mpfr_ldexp(fpart, exponent), &
      '3.141592653589793238462643383279502884197169399375105820974944592307816406286', &
      '1e-74', 'frexp ldexp')

    call mpfr_get_exponent_range(emin, emax)
    call mpfr_get_exponent_limits(emin_min, emin_max, emax_min, emax_max)
    if (emin < emin_min .or. emin > emin_max) error stop 'emin outside MPFR limits'
    if (emax < emax_min .or. emax > emax_max) error stop 'emax outside MPFR limits'
  end subroutine test_construction_and_math

  subroutine test_combinatorics()
    type(mpfr_real) :: n100

    call assert_close(mpfr_factorial(50, p), &
      '30414093201713378043612608166064768844377641568960512000000000000', '0', 'factorial 50')
    n100 = mpfr_from_integer(100_i64, p)
    call assert_close(mpfr_choose(n100, 50), '100891344545564193334812497256', '0', 'choose 100 50')
    call assert_close(mpfr_pochhammer(mpfr_from_string('1.5', p), 3), '13.125', '0', 'Pochhammer')
    call assert_close(mpfr_beta(mpfr_from_integer(2_i64, p), mpfr_from_integer(3_i64, p)), &
      '0.08333333333333333333333333333333333333333333333333333333333333333333333333333', &
      '1e-75', 'beta 2 3')
    call assert_close(mpfr_bernoulli(1, p), '0.5', '1e-75', 'Bernoulli 1')
    call assert_close(mpfr_bernoulli(2, p), &
      '0.1666666666666666666666666666666666666666666666666666666666666666666666666667', &
      '1e-74', 'Bernoulli 2')
  end subroutine test_combinatorics

  subroutine test_probability()
    type(mpfr_real) :: zero, one

    zero = mpfr_zero(1, p)
    one = mpfr_from_integer(1_i64, p)
    call assert_close(mpfr_pnorm(zero), '0.5', '1e-74', 'pnorm zero')
    call assert_close(mpfr_dnorm(zero), &
      '0.3989422804014326779399460599343818684758586311649346576659258296706579258993', &
      '1e-74', 'dnorm zero')
    call assert_close(mpfr_dt(one, mpfr_from_integer(7_i64, p)), &
      '0.2256749202754574839303391114936004922118830510090708570646796310838953745879', &
      '1e-72', 'Student t density')
    call assert_close(mpfr_dpois(mpfr_from_integer(3_i64, p), mpfr_from_string('2.5', p)), &
      '0.2137630172497364457539809230915619995776148984776996063691625273833978043761', &
      '1e-73', 'Poisson density')
    call assert_close(mpfr_dbinom(mpfr_from_integer(3_i64, p), mpfr_from_integer(10_i64, p), &
      mpfr_from_string('0.2', p)), '0.201326592', '1e-74', 'binomial density')
    call assert_close(mpfr_dnbinom(mpfr_from_integer(4_i64, p), mpfr_from_integer(2_i64, p), &
      mpfr_from_string('0.3', p)), '0.108045', '1e-74', 'negative binomial density')
    call assert_close(mpfr_dgamma(mpfr_from_integer(2_i64, p), mpfr_from_integer(3_i64, p), &
      scale=mpfr_from_string('1.5', p)), &
      '0.1562057114759862341209090048199525330581751599620003317595493300493813695771', &
      '1e-72', 'gamma density')
    call assert_close(mpfr_pgamma(mpfr_from_integer(2_i64, p), mpfr_from_integer(3_i64, p), &
      scale=mpfr_from_string('1.5', p)), &
      '0.1506314438493248519675572862915081014961725677066231960574505178564888029243', &
      '1e-72', 'gamma CDF')
    call assert_close(mpfr_pbeta_integer(mpfr_from_string('0.4', p), 2, 3), &
      '0.5248', '1e-74', 'integer beta CDF')
    call assert_close_value(mpfr_log1mexp(mpfr_const_log2(p)), mpfr_log(mpfr_from_string('0.5', p)), &
      '1e-74', 'log1mexp log2')
    call assert_close_value(mpfr_log1pexp(zero), mpfr_const_log2(p), '1e-74', 'log1pexp zero')
  end subroutine test_probability

  subroutine test_summary_matrix()
    allocate(values(3))
    values(1) = mpfr_from_string('1.25', p)
    values(2) = mpfr_from_string('-0.5', p)
    values(3) = mpfr_from_string('2.0', p)
    call assert_close(mpfr_sum(values), '2.75', '0', 'sum')
    call assert_close(mpfr_product(values), '-1.25', '0', 'product')
    call assert_close(mpfr_minimum(values), '-0.5', '0', 'min')
    call assert_close(mpfr_maximum(values), '2', '0', 'max')
    deallocate(values)

    allocate(a(2, 2), b(2, 2))
    a(1, 1) = mpfr_from_integer(1_i64, p)
    a(2, 1) = mpfr_from_integer(3_i64, p)
    a(1, 2) = mpfr_from_integer(2_i64, p)
    a(2, 2) = mpfr_from_integer(4_i64, p)
    b(1, 1) = mpfr_from_integer(5_i64, p)
    b(2, 1) = mpfr_from_integer(7_i64, p)
    b(1, 2) = mpfr_from_integer(6_i64, p)
    b(2, 2) = mpfr_from_integer(8_i64, p)
    call mpfr_matmul(a, b, c)
    call assert_close(c(1, 1), '19', '0', 'matmul 11')
    call assert_close(c(2, 2), '50', '0', 'matmul 22')
    deallocate(c)
    call mpfr_crossprod(a, c)
    call assert_close(c(1, 1), '10', '0', 'crossprod 11')
    call assert_close(c(2, 2), '20', '0', 'crossprod 22')
    deallocate(c)
    call mpfr_tcrossprod(a, c)
    call assert_close(c(1, 1), '5', '0', 'tcrossprod 11')
    call assert_close(c(2, 2), '25', '0', 'tcrossprod 22')
    deallocate(c, a, b)
  end subroutine test_summary_matrix

  subroutine test_algorithms()
    type(mpfr_real), allocatable :: start(:)

    integ = mpfr_integrate_romberg(sine_fun, mpfr_zero(1, p), mpfr_const_pi(p), &
      rel_tol=mpfr_from_string('1e-50', p), abs_tol=mpfr_from_string('1e-50', p), max_order=18)
    if (.not. integ%converged) error stop 'Romberg integration did not converge'
    call assert_close(integ%value, '2', '1e-45', 'Romberg sine integral')

    root = mpfr_uniroot(sqrt2_fun, mpfr_from_integer(1_i64, p), mpfr_from_integer(2_i64, p), &
      tol=mpfr_from_string('1e-60', p), maxiter=300)
    if (.not. root%converged) error stop 'Brent root finder did not converge'
    call assert_close(root%root, &
      '1.414213562373095048801688724209698078569671875376948073176679737990732478462', &
      '1e-55', 'Brent sqrt2 root')

    opt = mpfr_optimize(quadratic_fun, mpfr_from_integer(-5_i64, p), mpfr_from_integer(5_i64, p), &
      tol=mpfr_from_string('1e-50', p), method=optimize_brent, maxiter=500)
    if (.not. opt%converged) error stop 'Brent optimization did not converge'
    call assert_close(opt%location, '1.2345', '1e-38', 'Brent minimum location')
    call assert_close(opt%objective, '2', '1e-70', 'Brent minimum objective')

    opt = mpfr_optimize(quadratic_fun, mpfr_from_integer(-5_i64, p), mpfr_from_integer(5_i64, p), &
      tol=mpfr_from_string('1e-30', p), method=optimize_golden, maxiter=500)
    if (.not. opt%converged) error stop 'golden optimization did not converge'
    call assert_close(opt%location, '1.2345', '1e-28', 'golden minimum location')

    quant = mpfr_qnorm(mpfr_from_string('0.975', p), tol=mpfr_from_string('1e-60', p), maxiter=500)
    if (.not. quant%converged) error stop 'qnorm inversion did not converge'
    call assert_close(quant%root, &
      '1.959963984540054235524594430520551527955550077869548398476952646361635274145', &
      '1e-52', 'qnorm 0.975')

    allocate(start(2))
    start(1) = mpfr_zero(1, p)
    start(2) = mpfr_zero(1, p)
    call mpfr_hjk(start, sphere_fun, hj, tol=mpfr_from_string('1e-8', p), maxfeval=20000)
    if (.not. hj%converged) error stop 'Hooke-Jeeves did not converge'
    call assert_close(hj%par(1), '1', '2e-8', 'Hooke-Jeeves x1')
    call assert_close(hj%par(2), '-2', '2e-8', 'Hooke-Jeeves x2')
    call assert_close(hj%value, '0', '1e-14', 'Hooke-Jeeves objective')
    deallocate(start)
  end subroutine test_algorithms

  subroutine assert_close(actual, expected_text, tolerance_text, label)
    type(mpfr_real), intent(in) :: actual !! Computed arbitrary-precision value under test.
    character(len=*), intent(in) :: expected_text !! Decimal reference value parsed at the actual value's precision.
    character(len=*), intent(in) :: tolerance_text !! Nonnegative absolute tolerance in decimal notation.
    character(len=*), intent(in) :: label !! Short test label printed on failure.
    type(mpfr_real) :: expected, tolerance, difference
    integer :: prec_bits

    prec_bits = mpfr_precision(actual)
    expected = mpfr_from_string(expected_text, prec_bits)
    tolerance = mpfr_from_string(tolerance_text, prec_bits)
    difference = mpfr_abs(actual - expected)
    if (difference > tolerance) then
      write (*, '(a)') 'FAIL: ' // trim(label)
      write (*, '(a)') ' actual   = ' // trim(mpfr_to_string(actual, 60))
      write (*, '(a)') ' expected = ' // trim(mpfr_to_string(expected, 60))
      write (*, '(a)') ' |error|  = ' // trim(mpfr_to_string(difference, 20))
      error stop 'Rmpfr deterministic test failure'
    end if
  end subroutine assert_close

  subroutine assert_close_value(actual, expected, tolerance_text, label)
    type(mpfr_real), intent(in) :: actual !! Computed arbitrary-precision value under test.
    type(mpfr_real), intent(in) :: expected !! Arbitrary-precision reference value.
    character(len=*), intent(in) :: tolerance_text !! Nonnegative absolute tolerance in decimal notation.
    character(len=*), intent(in) :: label !! Short test label printed on failure.
    type(mpfr_real) :: tolerance, difference
    integer :: prec_bits

    prec_bits = max(mpfr_precision(actual), mpfr_precision(expected))
    tolerance = mpfr_from_string(tolerance_text, prec_bits)
    difference = mpfr_abs(mpfr_copy(actual, prec_bits) - mpfr_copy(expected, prec_bits))
    if (difference > tolerance) then
      write (*, '(a)') 'FAIL: ' // trim(label)
      error stop 'Rmpfr deterministic test failure'
    end if
  end subroutine assert_close_value

end program test_rmpfr
