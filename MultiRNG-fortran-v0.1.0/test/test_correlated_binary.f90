! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_correlated_binary
  use multirng, only : dp, seed_rng, draw_correlated_binary
  use test_support, only : assert_true, column_mean, sample_covariance
  implicit none
  real(dp) :: p(3), rho(3,3), denom
  integer, allocatable :: y(:,:)
  real(dp), allocatable :: yr(:,:), m(:), c(:,:)
  integer :: i,j
  integer, parameter :: nrep=120000

  call seed_rng(777)
  p = [0.3_dp, 0.5_dp, 0.7_dp]
  rho = reshape([1.0_dp,0.2_dp,0.3_dp, 0.2_dp,1.0_dp,0.2_dp, 0.3_dp,0.2_dp,1.0_dp], [3,3])
  y = draw_correlated_binary(nrep, 3, p, rho)
  allocate(yr(nrep,3))
  yr = real(y,dp)
  m = column_mean(yr)
  c = sample_covariance(yr)
  call assert_true(maxval(abs(m-p)) < 0.007_dp, 'correlated binary means')
  do i=1,3
    do j=1,3
      denom = sqrt(c(i,i)*c(j,j))
      call assert_true(abs(c(i,j)/denom-rho(i,j)) < 0.018_dp, 'correlated binary correlation')
    end do
  end do
  call assert_true(minval(y) == 0 .and. maxval(y) == 1, 'binary support')
  write(*,'(a)') 'test_correlated_binary: PASS'
end program test_correlated_binary
