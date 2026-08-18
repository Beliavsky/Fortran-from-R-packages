! SPDX-License-Identifier: MIT
program test_distribution
  use chyper, only : dp, chyper_probabilities, dchyper, pchyper
  implicit none
  integer, parameter :: n(3) = [12, 13, 14]
  integer, parameter :: m(3) = [7, 8, 9]
  real(dp), allocatable :: p(:)
  real(dp), parameter :: ref(0:7) = [ &
      0.6405485878096133_dp, 0.3072766936904082_dp, &
      0.0488723178299738_dp, 0.00321214699093537_dp, &
      8.92792711896144e-5_dp, 9.71011736475252e-7_dp, &
      3.39375471325742e-9_dp, 2.38828621622620e-12_dp ]
  integer :: i, status

  call chyper_probabilities(10, n, m, p, status)
  if (status /= 0) error stop 1
  if (lbound(p, 1) /= 0 .or. ubound(p, 1) /= 7) error stop 2
  if (maxval(abs(p - ref)) > 2.0e-14_dp) error stop 3
  if (abs(sum(p) - 1.0_dp) > 2.0e-15_dp) error stop 4
  do i = 0, 7
    if (abs(dchyper(i, 10, n, m) - ref(i)) > 2.0e-14_dp) error stop 5
  end do
  if (abs(pchyper(3, 10, n, m) - 0.9999097463209307_dp) > 3.0e-14_dp) error stop 6
  if (abs(dchyper(-1, 10, n, m)) > tiny(1.0_dp)) error stop 7
  if (abs(pchyper(-1, 10, n, m)) > tiny(1.0_dp)) error stop 8
  if (abs(pchyper(99, 10, n, m) - 1.0_dp) > tiny(1.0_dp)) error stop 9
  print *, 'test_distribution: PASS'
end program test_distribution
