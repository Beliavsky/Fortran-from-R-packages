! SPDX-License-Identifier: MIT
program test_rng
  use chyper, only : dp, chyper_probabilities, chyper_seed, rchyper
  implicit none
  integer, parameter :: n(3) = [12, 13, 14]
  integer, parameter :: m(3) = [7, 8, 9]
  integer, allocatable :: x(:)
  real(dp), allocatable :: p(:)
  real(dp) :: exact_mean, sample_mean
  integer :: i, status

  allocate(x(40000))
  call chyper_probabilities(10, n, m, p, status)
  exact_mean = 0.0_dp
  do i = 0, ubound(p, 1)
    exact_mean = exact_mean + real(i, dp) * p(i)
  end do
  call chyper_seed(12345)
  call rchyper(size(x), 10, n, m, x)
  if (minval(x) < 0 .or. maxval(x) > 7) error stop 1
  sample_mean = sum(real(x, dp)) / real(size(x), dp)
  if (abs(sample_mean - exact_mean) > 0.02_dp) error stop 2
  print *, 'test_rng: PASS'
end program test_rng
