program test_arithmetic
  use polynom
  implicit none
  type(polynomial_t) :: p, q, r, quotient, remainder
  complex(dp), allocatable :: roots(:)
  real(dp) :: x(5), y(5)
  integer :: i

  p = polynomial([1.0_dp, 2.0_dp, 1.0_dp])
  q = polynomial([0.0_dp, -1.0_dp, 0.0_dp, 1.0_dp])
  r = (q - 2.0_dp * p)**2
  call assert_close(r%coef, [4.0_dp, 20.0_dp, 33.0_dp, 16.0_dp, -6.0_dp, -4.0_dp, 1.0_dp], 1.0e-12_dp)

  call poly_divmod(q, p, quotient, remainder)
  call assert_close(quotient%coef, [-2.0_dp, 1.0_dp], 1.0e-12_dp)
  call assert_close(remainder%coef, [2.0_dp, 2.0_dp], 1.0e-12_dp)
  r = p * quotient + remainder
  call assert_close(r%coef, q%coef, 1.0e-12_dp)

  x = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  y = x*x + 1.0_dp
  r = poly_from_values(x, y)
  call assert_close(r%coef, [1.0_dp, 0.0_dp, 1.0_dp], 1.0e-11_dp)

  p = poly_from_roots([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])
  roots = polynomial_roots(p)
  if (size(roots) /= 5) error stop 'wrong root count'
  do i = 1, size(roots)
    if (abs(p%evaluate(real(roots(i), dp))) > 1.0e-7_dp .or. abs(aimag(roots(i))) > 1.0e-7_dp) then
      error stop 'root verification failed'
    end if
  end do
  roots = polynomial_roots(polynomial([1.0_dp, 0.0_dp, 1.0_dp]))
  if (size(roots) /= 2 .or. maxval(abs(abs(roots) - 1.0_dp)) > 1.0e-9_dp) then
    error stop 'complex roots failed'
  end if
  if (maxval(abs(real(roots, dp))) > 1.0e-9_dp) error stop 'complex root real parts failed'

  print *, 'test_arithmetic: PASS'
contains
  subroutine assert_close(a, b, tol)
    real(dp), intent(in) :: a(:), b(:), tol
    if (size(a) /= size(b)) error stop 'shape mismatch'
    if (maxval(abs(a - b)) > tol) error stop 'values differ'
  end subroutine assert_close
end program test_arithmetic
