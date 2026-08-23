program test_special
  use iso_c_binding, only: c_double, c_int
  use gsl_special, only: airy_ai, gsl_sf_gamma, erf, bessel_cyl_j0, lambert_w0, zeta
  implicit none
  real(c_double), target :: x(3), v(3), e(3)
  integer(c_int), target :: st(3), mode

  x = [0.0_c_double, 1.0_c_double, 2.0_c_double]
  mode = 0_c_int
  call airy_ai(x, mode, v, e, st)
  call assert_close(v(1), 0.3550280538878172_c_double, 2.0e-14_c_double)

  call gsl_sf_gamma(x + 1.0_c_double, v, e, st)
  call assert_close(v(1), 1.0_c_double, 2.0e-14_c_double)
  call assert_close(v(3), 2.0_c_double, 2.0e-14_c_double)

  call erf(x, mode, v, e, st)
  call assert_close(v(2), 0.8427007929497149_c_double, 2.0e-14_c_double)

  call bessel_cyl_j0(x, v, e, st)
  call assert_close(v(1), 1.0_c_double, 2.0e-14_c_double)

  x = [0.0_c_double, 1.0_c_double, 2.0_c_double]
  call lambert_w0(x, v, e, st)
  call assert_close(v(2), 0.5671432904097839_c_double, 2.0e-14_c_double)

  x = [2.0_c_double, 3.0_c_double, 4.0_c_double]
  call zeta(x, v, e, st)
  call assert_close(v(1), 1.6449340668482264_c_double, 3.0e-14_c_double)

  print '(a)', 'test_special: PASS'
contains
  subroutine assert_close(a, b, tol)
    real(c_double), intent(in) :: a, b, tol
    if (abs(a - b) > tol) error stop 1
  end subroutine assert_close
end program test_special
