! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_wishart
  use multirng, only : dp, seed_rng, draw_wishart, draw_inv_wishart, draw_inv_wishart_legacy
  use test_support, only : assert_true
  implicit none
  real(dp) :: sigma(2,2), inv_sigma(2,2), meanw(2,2), meaniw(2,2), targetiw(2,2)
  real(dp), allocatable :: w(:,:,:), iw(:,:,:), legacy(:,:,:)
  integer :: i
  integer, parameter :: nrep=25000, nu=8

  call seed_rng(314159)
  sigma = reshape([1.0_dp, 0.25_dp, 0.25_dp, 1.5_dp], [2,2])
  inv_sigma = reshape([1.0434782608695652_dp, -0.1739130434782609_dp, &
                      -0.1739130434782609_dp, 0.6956521739130435_dp], [2,2])
  w = draw_wishart(nrep, 2, nu, sigma)
  meanw = 0.0_dp
  do i=1,nrep
    meanw = meanw + w(i,:,:)
  end do
  meanw = meanw / real(nrep,dp)
  call assert_true(maxval(abs(meanw - real(nu,dp)*sigma)) < 0.11_dp, 'Wishart mean')

  iw = draw_inv_wishart(nrep, 2, nu, inv_sigma)
  meaniw = 0.0_dp
  do i=1,nrep
    meaniw = meaniw + iw(i,:,:)
  end do
  meaniw = meaniw / real(nrep,dp)
  targetiw = sigma / real(nu - 2 - 1, dp)
  call assert_true(maxval(abs(meaniw - targetiw)) < 0.015_dp, 'inverse-Wishart mean')

  legacy = draw_inv_wishart_legacy(20, 2, nu, inv_sigma)
  call assert_true(all(legacy(:,1,1) > 0.0_dp), 'legacy inverse-Wishart compatibility path')
  write(*,'(a)') 'test_wishart: PASS'
end program test_wishart
