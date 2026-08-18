! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program example_multirng
  use multirng, only : dp, seed_rng
  use multirng, only : draw_d_variate_normal, draw_dirichlet
  use multirng, only : draw_correlated_binary, draw_wishart
  implicit none
  real(dp) :: mu(3), sigma(3,3), alpha(3), p(3), rho(3,3)
  real(dp), allocatable :: x(:,:), dir(:,:), w(:,:,:)
  integer, allocatable :: y(:,:)
  integer :: j

  call seed_rng(20260817)

  mu = [0.0_dp, 3.0_dp, 7.0_dp]
  sigma = reshape([1.0_dp,0.2_dp,0.3_dp, 0.2_dp,1.0_dp,0.2_dp, 0.3_dp,0.2_dp,1.0_dp], [3,3])
  x = draw_d_variate_normal(10000, 3, mu, sigma)
  write(*,'(a)') 'Multivariate normal sample means:'
  write(*,'(3f12.6)') [(sum(x(:,j))/real(size(x,1),dp), j=1,3)]

  alpha = [1.0_dp, 3.0_dp, 4.0_dp]
  dir = draw_dirichlet(5, 3, alpha, 2.0_dp)
  write(*,'(/,a)') 'Five Dirichlet draws:'
  do j=1,size(dir,1)
    write(*,'(3f12.6)') dir(j,:)
  end do

  p = [0.3_dp, 0.5_dp, 0.7_dp]
  rho = sigma
  y = draw_correlated_binary(5, 3, p, rho)
  write(*,'(/,a)') 'Five correlated Bernoulli draws:'
  do j=1,size(y,1)
    write(*,'(3i4)') y(j,:)
  end do

  w = draw_wishart(1, 3, 5, sigma)
  write(*,'(/,a)') 'One Wishart draw:'
  do j=1,3
    write(*,'(3f12.6)') w(1,j,:)
  end do
end program example_multirng
