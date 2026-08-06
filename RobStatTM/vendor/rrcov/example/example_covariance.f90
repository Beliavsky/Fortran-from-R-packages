! SPDX-License-Identifier: GPL-3.0-or-later
program example_covariance
  use rrcov, only : dp, covariance_result, robust_covariance
  implicit none
  real(dp) :: x(30, 2)
  type(covariance_result) :: estimate
  integer :: i
  do i = 1, 30
    x(i, 1) = sin(0.2_dp * real(i, dp))
    x(i, 2) = 0.8_dp * x(i, 1) + 0.2_dp * cos(0.3_dp * real(i, dp))
  end do
  x(1, :) = [20.0_dp, -20.0_dp]
  call robust_covariance(x, "mcd", estimate, nsamp=100, seed=123)
  print '(a,2f12.6)', "center: ", estimate%center
  print '(a)', "covariance:"
  print '(2f12.6)', estimate%covariance(1, :)
  print '(2f12.6)', estimate%covariance(2, :)
end program example_covariance
