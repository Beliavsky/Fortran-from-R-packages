program test_gcd_summary
  use polynom
  implicit none
  type(polynomial_t) :: a, b, c, g, l
  type(polynomial_summary_t) :: s
  type(polylist_t) :: list

  a = poly_from_roots([-1.0_dp])
  b = poly_from_roots([-1.0_dp, -1.0_dp])
  c = poly_from_roots([1.0_dp])
  g = polynomial_gcd(a, b)
  call assert_close(g%coef, a%coef, 1.0e-10_dp)
  l = polynomial_lcm(a, b)
  call assert_close(l%coef, b%coef, 1.0e-10_dp)

  allocate(list%item(3))
  list%item = [a, b, c]
  g = polylist_gcd(list)
  call assert_close(g%coef, [1.0_dp], 1.0e-10_dp)
  l = polylist_lcm(list)
  call assert_close(l%coef, [-1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp], 1.0e-9_dp)

  s = summarize_polynomial(polynomial([-1.0_dp, 0.0_dp, 1.0_dp]))
  if (size(s%zeros) /= 2 .or. size(s%stationary_points) /= 1) error stop 'summary sizes failed'
  if (maxval(abs(abs(s%zeros) - 1.0_dp)) > 1.0e-8_dp) error stop 'summary roots failed'
  if (abs(s%stationary_points(1)) > 1.0e-8_dp) error stop 'stationary point failed'

  print *, 'test_gcd_summary: PASS'
contains
  subroutine assert_close(x, y, tol)
    real(dp), intent(in) :: x(:), y(:), tol
    if (size(x) /= size(y)) error stop 'shape mismatch'
    if (maxval(abs(x-y)) > tol) error stop 'values differ'
  end subroutine assert_close
end program test_gcd_summary
