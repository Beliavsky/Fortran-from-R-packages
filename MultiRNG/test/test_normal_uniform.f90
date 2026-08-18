! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_normal_uniform
  use multirng, only : dp, seed_rng, draw_d_variate_normal, draw_d_variate_uniform
  use test_support, only : assert_true, column_mean, sample_covariance
  implicit none
  real(dp) :: mu(3), sigma(3,3)
  real(dp), allocatable :: x(:,:), u(:,:), m(:), c(:,:)
  integer, parameter :: n = 50000

  call seed_rng(12345)
  mu = [1.0_dp, -2.0_dp, 0.5_dp]
  sigma = reshape([1.0_dp, 0.3_dp, -0.2_dp, &
                   0.3_dp, 2.0_dp, 0.4_dp, &
                  -0.2_dp, 0.4_dp, 1.5_dp], [3,3])
  x = draw_d_variate_normal(n, 3, mu, sigma)
  m = column_mean(x)
  c = sample_covariance(x)
  call assert_true(maxval(abs(m - mu)) < 0.025_dp, 'multivariate normal means')
  call assert_true(maxval(abs(c - sigma)) < 0.045_dp, 'multivariate normal covariance')

  u = draw_d_variate_uniform(n, 3, sigma)
  m = column_mean(u)
  call assert_true(maxval(abs(m - 0.5_dp)) < 0.01_dp, 'Gaussian-copula uniform margins')
  call assert_true(minval(u) > 0.0_dp .and. maxval(u) < 1.0_dp, 'uniform support')
  write(*,'(a)') 'test_normal_uniform: PASS'
end program test_normal_uniform
