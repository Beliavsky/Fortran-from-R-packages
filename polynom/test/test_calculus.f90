program test_calculus
  use polynom
  implicit none
  type(polynomial_t) :: p, d, a, shifted
  real(dp) :: values(8)

  p = poly_from_roots([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])
  call assert_close(p%coef, [-120.0_dp, 274.0_dp, -225.0_dp, 85.0_dp, -15.0_dp, 1.0_dp], 1.0e-12_dp)
  d = derivative(p)
  call assert_close(d%coef, [274.0_dp, -450.0_dp, 255.0_dp, -60.0_dp, 5.0_dp], 1.0e-12_dp)
  a = integral_polynomial(d, -120.0_dp)
  call assert_close(a%coef, p%coef, 1.0e-12_dp)
  if (abs(definite_integral(d, 0.0_dp, 1.0_dp) - (p%evaluate(1.0_dp) - p%evaluate(0.0_dp))) > 1.0e-12_dp) then
    error stop 'definite integral failed'
  end if

  shifted = change_origin(p, 3.0_dp)
  call assert_close(shifted%coef, [0.0_dp, 4.0_dp, 0.0_dp, -5.0_dp, 0.0_dp, 1.0_dp], 1.0e-11_dp)
  values = p%evaluate([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp])
  if (maxval(abs(values - [-120.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 120.0_dp, 720.0_dp])) > 1.0e-8_dp) then
    error stop 'evaluation failed'
  end if

  print *, 'test_calculus: PASS'
contains
  subroutine assert_close(x, y, tol)
    real(dp), intent(in) :: x(:), y(:), tol
    if (size(x) /= size(y)) error stop 'shape mismatch'
    if (maxval(abs(x-y)) > tol) error stop 'values differ'
  end subroutine assert_close
end program test_calculus
