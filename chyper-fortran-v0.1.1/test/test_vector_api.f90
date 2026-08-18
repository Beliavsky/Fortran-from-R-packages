! SPDX-License-Identifier: MIT
program test_vector_api
  use chyper, only : dp, dchyper_vec, pchyper_vec, qchyper_vec, mle_n
  implicit none
  integer, parameter :: n(3) = [12, 13, 14]
  integer, parameter :: m(3) = [7, 8, 9]
  integer, parameter :: k(4) = [-1, 0, 3, 8]
  real(dp), parameter :: pp(4) = [-0.1_dp, 0.5_dp, 0.9_dp, 1.1_dp]
  real(dp) :: d(4), c(4), inf_est
  integer :: q(4)
  integer, parameter :: zeros(4) = [0, 0, 0, 0]
  integer :: nt(3)

  call dchyper_vec(k, 10, n, m, d)
  call pchyper_vec(k, 10, n, m, c)
  call qchyper_vec(pp, 10, n, m, q)
  if (abs(d(1)) > tiny(1.0_dp)) error stop 1
  if (abs(d(2) - 0.6405485878096133_dp) > 2.0e-14_dp) error stop 2
  if (abs(d(3) - 0.00321214699093537_dp) > 2.0e-14_dp) error stop 3
  if (abs(d(4)) > tiny(1.0_dp)) error stop 4
  if (abs(c(1)) > tiny(1.0_dp)) error stop 5
  if (abs(c(4) - 1.0_dp) > tiny(1.0_dp)) error stop 6
  if (any(q /= [0, 0, 1, 7])) error stop 7
  nt = [0, 13, 14]
  inf_est = mle_n(1, zeros, 8, nt, m)
  if (inf_est < 0.5_dp * huge(1.0_dp)) error stop 8
  print *, 'test_vector_api: PASS'
end program test_vector_api
