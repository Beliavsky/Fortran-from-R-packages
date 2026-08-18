! SPDX-License-Identifier: MIT
program test_gem
  use ewens, only : dp, i8, ewens_seed, rgem
  implicit none
  real(dp), allocatable :: s(:)
  real(dp) :: average_first, expected_first
  integer :: i, reps

  call ewens_seed(424242_i8)
  reps = 10000
  average_first = 0.0_dp
  do i = 1, reps
    s = rgem(alpha=0.25_dp, theta=1.5_dp, trunc_at=20)
    if (any(s < 0.0_dp) .or. sum(s) > 1.0_dp + 2.0e-14_dp) error stop 1
    average_first = average_first + s(1)
  end do
  average_first = average_first / real(reps, dp)
  expected_first = (1.0_dp - 0.25_dp) / (1.0_dp + 1.5_dp)
  if (abs(average_first - expected_first) > 0.012_dp) then
    print '(a,2f14.7)', 'GEM first-share mismatch: ', average_first, expected_first
    error stop 2
  end if

  s = rgem(alpha=0.0_dp, theta=0.0_dp, trunc_at=5)
  if (abs(s(1) - 1.0_dp) > 0.0_dp) error stop 3
  if (any(abs(s(2:)) > 0.0_dp)) error stop 4

  print '(a)', 'test_gem: PASS'
end program test_gem
