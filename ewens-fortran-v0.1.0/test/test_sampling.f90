! SPDX-License-Identifier: MIT
program test_sampling
  use ewens, only : dp, i8, ewens_seed, rewens, gcrp, number_of_classes, ewens_k_exact
  implicit none
  integer, allocatable :: a(:), b(:), x(:)
  integer :: i, reps
  real(dp) :: average_k, target

  call ewens_seed(1234567_i8)
  a = rewens(80, 1.7_dp)
  call ewens_seed(1234567_i8)
  b = gcrp(80, 0.0_dp, 1.7_dp)
  if (any(a /= b)) error stop 1

  a = rewens(10, 0.0_dp)
  if (any(a /= 1)) error stop 2

  x = gcrp(200, 0.35_dp, 1.2_dp)
  if (minval(x) < 1 .or. maxval(x) > 200) error stop 3
  if (number_of_classes(x) /= maxval(x)) error stop 4

  call ewens_seed(998877_i8)
  reps = 8000
  average_k = 0.0_dp
  do i = 1, reps
    x = rewens(30, 1.5_dp)
    average_k = average_k + real(number_of_classes(x), dp)
  end do
  average_k = average_k / real(reps, dp)
  target = ewens_k_exact(30, 1.5_dp)
  if (abs(average_k - target) > 0.07_dp) then
    print '(a,2f14.7)', 'sampling mean mismatch: ', average_k, target
    error stop 5
  end if

  x = rewens(0, 1.0_dp)
  if (size(x) /= 0) error stop 6

  print '(a)', 'test_sampling: PASS'
end program test_sampling
