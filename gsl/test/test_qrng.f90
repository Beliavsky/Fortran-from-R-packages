program test_qrng
  use iso_c_binding, only: c_double
  use gsl_qrng
  implicit none
  type(gsl_qrng_type) :: q
  real(c_double) :: x(3,8)
  integer :: status

  q = qrng_alloc(qrng_sobol, 3)
  if (.not. qrng_is_allocated(q)) error stop 1
  status = qrng_get_n(q, x)
  if (status /= 0) error stop 2
  if (any(x < 0.0_c_double) .or. any(x >= 1.0_c_double)) error stop 3
  if (maxval(abs(x(:,1) - 0.5_c_double)) > 1.0e-15_c_double) error stop 4
  if (len_trim(qrng_name(q)) == 0) error stop 5
  call qrng_free(q)

  print '(a)', 'test_qrng: PASS'
end program test_qrng
