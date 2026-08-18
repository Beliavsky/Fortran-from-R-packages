! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_t_locmin
  use multirng, only : dp, seed_rng, draw_d_variate_t, loc_min
  use test_support, only : assert_true, column_mean
  implicit none
  real(dp) :: mu(2), sigma(2,2), a(3,3)
  real(dp), allocatable :: x(:,:), m(:)
  integer :: where(2)

  call seed_rng(9876)
  mu=[2.0_dp,-1.0_dp]
  sigma=reshape([1.0_dp,0.4_dp,0.4_dp,1.5_dp],[2,2])
  x=draw_d_variate_t(7.0_dp,40000,2,mu,sigma)
  m=column_mean(x)
  call assert_true(maxval(abs(m-mu)) < 0.04_dp, 'multivariate t means')

  a=reshape([0.0_dp,0.4_dp,0.3_dp, 0.4_dp,0.2_dp,0.6_dp, 0.3_dp,0.6_dp,0.5_dp],[3,3])
  where=loc_min(a,3)
  call assert_true(all(where == [2,2]), 'loc_min ignores zeros and finds minimum')
  write(*,'(a)') 'test_t_locmin: PASS'
end program test_t_locmin
