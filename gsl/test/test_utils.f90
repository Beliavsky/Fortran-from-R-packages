program test_utils
  use iso_c_binding, only: c_double, c_int
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use gsl_utils, only: gsl_poly, strictify, gsl_sn, gsl_cn, gsl_dn
  implicit none
  real(c_double) :: c(3), x(3), y(3)
  integer(c_int) :: st(3)
  complex(c_double) :: sn, cn, dn

  c = [1.0_c_double, 2.0_c_double, 3.0_c_double]
  x = [-1.0_c_double, 0.0_c_double, 2.0_c_double]
  y = gsl_poly(c, x)
  if (maxval(abs(y - [2.0_c_double, 1.0_c_double, 17.0_c_double])) > 1.0e-14_c_double) error stop 1

  st = [0_c_int, 1_c_int, 0_c_int]
  call strictify(y, st)
  if (.not. ieee_is_nan(y(2))) error stop 2

  sn = gsl_sn(cmplx(0.4_c_double, 0.2_c_double, c_double), 0.3_c_double)
  cn = gsl_cn(cmplx(0.4_c_double, 0.2_c_double, c_double), 0.3_c_double)
  dn = gsl_dn(cmplx(0.4_c_double, 0.2_c_double, c_double), 0.3_c_double)
  if (abs(sn**2 + cn**2 - 1.0_c_double) > 2.0e-12_c_double) error stop 3
  if (abs(dn**2 + 0.3_c_double * sn**2 - 1.0_c_double) > 2.0e-12_c_double) error stop 4

  print '(a)', 'test_utils: PASS'
end program test_utils
