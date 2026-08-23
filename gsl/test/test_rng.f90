program test_rng
  use iso_c_binding, only: c_long, c_double
  use gsl_rng
  implicit none
  type(gsl_rng_type) :: r
  real(c_double) :: u(1000)
  integer(c_long) :: k(100)

  r = rng_alloc(rng_mt19937)
  if (.not. rng_is_allocated(r)) error stop 1
  call rng_set(r, 12345_c_long)
  call rng_uniform_array(r, u)
  if (any(u < 0.0_c_double) .or. any(u >= 1.0_c_double)) error stop 2
  if (abs(sum(u) / size(u) - 0.5_c_double) > 0.05_c_double) error stop 3
  call rng_uniform_int_array(r, 17_c_long, k)
  if (any(k < 0_c_long) .or. any(k >= 17_c_long)) error stop 4
  if (len_trim(rng_name(r)) == 0) error stop 5
  call rng_free(r)

  print '(a)', 'test_rng: PASS'
end program test_rng
