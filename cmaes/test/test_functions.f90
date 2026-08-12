program test_functions
  use iso_fortran_env, only : int64
  use cmaes_kinds, only : dp
  use cmaes_functions, only : f_sphere, f_rand, seed_f_rand, f_rosenbrock, f_rastrigin, &
    shifted_value, rotated_value, biased_value
  implicit none
  real(dp) :: x(3), off(3), m(3, 3), r1, r2

  x = [1.0_dp, 2.0_dp, 3.0_dp]
  if (abs(f_sphere(x) - 14.0_dp) > 1.0e-14_dp) error stop "sphere helper"
  if (abs(f_rastrigin([0.0_dp, 0.0_dp])) > 1.0e-14_dp) error stop "rastrigin helper"
  if (abs(f_rosenbrock([0.0_dp, 0.0_dp])) > 1.0e-14_dp) error stop "rosenbrock helper"
  off = 1.0_dp
  if (abs(shifted_value(f_sphere, x, off) - 5.0_dp) > 1.0e-14_dp) error stop "shift helper"
  m = 0.0_dp
  m(1, 2) = 1.0_dp
  m(2, 1) = 1.0_dp
  m(3, 3) = 1.0_dp
  if (abs(rotated_value(f_sphere, x, m) - 14.0_dp) > 1.0e-14_dp) error stop "rotate helper"
  if (abs(biased_value(f_sphere, x, 2.5_dp) - 16.5_dp) > 1.0e-14_dp) error stop "bias helper"
  call seed_f_rand(99_int64)
  r1 = f_rand(x)
  call seed_f_rand(99_int64)
  r2 = f_rand(x)
  if (r1 < r2 .or. r1 > r2) error stop "f_rand seed reproducibility"
  print '(a)', 'function helpers: PASS'
end program test_functions
