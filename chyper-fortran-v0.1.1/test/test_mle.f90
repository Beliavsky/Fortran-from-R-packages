! SPDX-License-Identifier: MIT
program test_mle
  use chyper, only : dp, mle_s, mle_n, mle_m
  implicit none
  integer, parameter :: k(7) = [0, 0, 1, 1, 0, 2, 0]
  integer, parameter :: n(3) = [12, 13, 14]
  integer, parameter :: m(3) = [7, 8, 9]
  integer :: mt(3), nt(3)
  real(dp) :: nhat

  if (mle_s(k, n, m) /= 7) error stop 1
  nt = [0, 13, 14]
  nhat = mle_n(1, k, 8, nt, m)
  if (abs(nhat - 7.0_dp) > 0.5_dp) error stop 2
  mt = [0, 8, 9]
  if (mle_m(1, k, 8, n, mt) /= 9) error stop 3
  print *, 'test_mle: PASS'
end program test_mle
