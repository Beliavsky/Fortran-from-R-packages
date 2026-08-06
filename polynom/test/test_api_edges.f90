program test_api_edges
  use polynom
  implicit none
  type(polynomial_t) :: p, q
  type(polylist_t) :: list, fits, transformed
  type(poly_status_t) :: status
  real(dp) :: x(4), y(4,2)

  p = polynomial([1.23456_dp, -2.71828_dp, 0.0_dp])
  q = round_coefficients(p, 2)
  call assert_close(q%coef, [1.23_dp, -2.72_dp], 1.0e-12_dp)
  q = monic(polynomial(0.0_dp), status)
  if (status%succeeded()) error stop 'zero monic should fail'

  x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
  y(:,1) = x
  y(:,2) = x*x + 1.0_dp
  fits = poly_from_values_matrix(x, y)
  if (fits%size() /= 2) error stop 'matrix interpolation failed'
  call assert_close(fits%item(1)%coef, [0.0_dp, 1.0_dp], 1.0e-11_dp)
  call assert_close(fits%item(2)%coef, [1.0_dp, 0.0_dp, 1.0_dp], 1.0e-11_dp)

  allocate(list%item(3))
  list%item(1) = polynomial(1.0_dp)
  list%item(2) = polynomial([0.0_dp, 1.0_dp])
  list%item(3) = polynomial([1.0_dp, 1.0_dp])
  q = sum_polynomials(list)
  call assert_close(q%coef, [2.0_dp, 2.0_dp], 1.0e-12_dp)
  q = product_polynomials(list)
  call assert_close(q%coef, [0.0_dp, 1.0_dp, 1.0_dp], 1.0e-12_dp)
  transformed = derivative_polylist(list)
  call assert_close(transformed%item(3)%coef, [1.0_dp], 1.0e-12_dp)
  transformed = integral_polylist(transformed)
  call assert_close(transformed%item(3)%coef, [0.0_dp, 1.0_dp], 1.0e-12_dp)
  if (len(p%to_string()) == 0) error stop 'string formatting failed'

  print *, 'test_api_edges: PASS'
contains
  subroutine assert_close(a, b, tol)
    real(dp), intent(in) :: a(:), b(:), tol
    if (size(a) /= size(b)) error stop 'shape mismatch'
    if (maxval(abs(a-b)) > tol) error stop 'values differ'
  end subroutine assert_close
end program test_api_edges
