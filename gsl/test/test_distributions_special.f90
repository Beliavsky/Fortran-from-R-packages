program test_distributions_special
  use iso_c_binding, only: c_double, c_int
  use gsl_special, only: beta_inc, hyperg_2f1, ellint_kcomp, legendre_p2
  implicit none
  real(c_double), target :: a(2), b(2), c(2), x(2), v(2), e(2)
  integer(c_int), target :: st(2), mode

  a = [2.0_c_double, 1.0_c_double]
  b = [3.0_c_double, 1.0_c_double]
  x = [0.5_c_double, 0.25_c_double]
  call beta_inc(a, b, x, v, e, st)
  call assert_close(v(1), 0.6875_c_double, 3.0e-14_c_double)
  call assert_close(v(2), 0.25_c_double, 3.0e-14_c_double)

  a = [1.0_c_double, 2.0_c_double]
  b = [1.0_c_double, 3.0_c_double]
  c = [2.0_c_double, 4.0_c_double]
  x = 0.0_c_double
  call hyperg_2f1(a, b, c, x, v, e, st)
  if (any(abs(v - 1.0_c_double) > 3.0e-14_c_double)) error stop 2

  x = [0.0_c_double, 0.5_c_double]
  mode = 0_c_int
  call ellint_kcomp(x, mode, v, e, st)
  call assert_close(v(1), acos(-1.0_c_double) / 2.0_c_double, 3.0e-14_c_double)

  x = [-0.4_c_double, 0.25_c_double]
  call legendre_p2(x, v, e, st)
  call assert_close(v(1), 0.5_c_double * (3.0_c_double * x(1)**2 - 1.0_c_double), 3.0e-14_c_double)

  print '(a)', 'test_distributions_special: PASS'
contains
  subroutine assert_close(u, w, tol)
    real(c_double), intent(in) :: u, w, tol
    if (abs(u - w) > tol) error stop 1
  end subroutine assert_close
end program test_distributions_special
