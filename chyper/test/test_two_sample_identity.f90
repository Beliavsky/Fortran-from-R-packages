! SPDX-License-Identifier: MIT
program test_two_sample_identity
  use chyper, only : dp, chyper_probabilities
  implicit none
  integer, parameter :: n(2) = [9, 11]
  integer, parameter :: m(2) = [6, 7]
  real(dp), allocatable :: p(:), pr(:)
  integer :: status

  call chyper_probabilities(8, n, m, p, status)
  if (status /= 0) error stop 1
  call chyper_probabilities(8, n(2:1:-1), m(2:1:-1), pr, status)
  if (size(p) /= size(pr)) error stop 2
  if (maxval(abs(p - pr)) > 5.0e-14_dp) error stop 3
  if (abs(sum(p) - 1.0_dp) > 2.0e-15_dp) error stop 4
  print *, 'test_two_sample_identity: PASS'
end program test_two_sample_identity
