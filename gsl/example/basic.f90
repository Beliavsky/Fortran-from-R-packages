program basic
  use iso_c_binding, only: c_double, c_int
  use gsl, only: airy_ai, bessel_cyl_j0, gsl_sf_gamma
  implicit none
  real(c_double), target :: x(1), value(1), error(1)
  integer(c_int), target :: status(1), mode

  x = 1.0_c_double
  mode = 0_c_int
  call airy_ai(x, mode, value, error, status)
  print '(a,f14.10)', 'Ai(1)      = ', value(1)
  call bessel_cyl_j0(x, value, error, status)
  print '(a,f14.10)', 'J0(1)      = ', value(1)
  call gsl_sf_gamma(x + 4.0_c_double, value, error, status)
  print '(a,f14.6)', 'Gamma(5)   = ', value(1)
end program basic
