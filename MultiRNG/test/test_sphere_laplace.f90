! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_sphere_laplace
  use multirng, only : dp, seed_rng, generate_point_in_sphere, draw_multivariate_laplace
  use test_support, only : assert_true, column_mean, sample_covariance
  implicit none
  real(dp) :: mu(3), sigma(3,3)
  real(dp), allocatable :: s(:,:), x(:,:), m(:), c(:,:)
  integer :: i
  integer, parameter :: nrep=50000

  call seed_rng(271828)
  s = generate_point_in_sphere(10000, 5)
  do i=1,size(s,1)
    call assert_true(abs(sum(s(i,:)*s(i,:)) - 1.0_dp) < 2.0e-12_dp, 'unit sphere norm')
  end do

  mu = [0.0_dp, 3.0_dp, 7.0_dp]
  sigma = reshape([1.0_dp,0.2_dp,0.3_dp, 0.2_dp,1.0_dp,0.2_dp, 0.3_dp,0.2_dp,1.0_dp], [3,3])
  x = draw_multivariate_laplace(nrep, 3, 2.0_dp, mu, sigma)
  m = column_mean(x)
  c = sample_covariance(x)
  call assert_true(maxval(abs(m-mu)) < 0.025_dp, 'multivariate Laplace means')
  call assert_true(maxval(abs(c-sigma)) < 0.045_dp, 'gamma=2 multivariate Laplace covariance')
  write(*,'(a)') 'test_sphere_laplace: PASS'
end program test_sphere_laplace
