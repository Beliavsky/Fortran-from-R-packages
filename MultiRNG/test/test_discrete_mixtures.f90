! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_discrete_mixtures
  use multirng, only : dp, seed_rng, draw_dirichlet, draw_multinomial
  use multirng, only : draw_dirichlet_multinomial, draw_multivariate_hypergeometric
  use test_support, only : assert_true, column_mean
  implicit none
  real(dp) :: alpha(4), theta(4)
  real(dp), allocatable :: dir(:,:), m(:), xr(:,:)
  integer, allocatable :: multi(:,:), dm(:,:), hg(:,:)
  integer :: counts(4), i
  integer, parameter :: nrep = 40000

  call seed_rng(24680)
  alpha = [1.0_dp, 3.0_dp, 4.0_dp, 4.0_dp]
  theta = alpha / sum(alpha)
  dir = draw_dirichlet(nrep, 4, alpha, 2.0_dp)
  m = column_mean(dir)
  call assert_true(maxval(abs(m - theta)) < 0.006_dp, 'Dirichlet means')
  call assert_true(maxval(abs(sum(dir, dim=2) - 1.0_dp)) < 1.0e-12_dp, 'Dirichlet row sums')

  multi = draw_multinomial(nrep, 4, theta, 12)
  allocate(xr(nrep,4))
  xr = real(multi, dp)
  m = column_mean(xr)
  call assert_true(maxval(abs(m - 12.0_dp * theta)) < 0.04_dp, 'multinomial means')
  call assert_true(all(sum(multi, dim=2) == 12), 'multinomial row sums')

  dm = draw_dirichlet_multinomial(nrep, 4, alpha, 2.0_dp, 8)
  xr = real(dm, dp)
  m = column_mean(xr)
  call assert_true(maxval(abs(m - 8.0_dp * theta)) < 0.05_dp, 'upstream Dirichlet-multinomial means')
  call assert_true(all(sum(dm, dim=2) == 8), 'Dirichlet-multinomial row sums')

  counts = [10, 20, 30, 40]
  hg = draw_multivariate_hypergeometric(nrep, 4, counts, 15)
  xr = real(hg, dp)
  m = column_mean(xr)
  do i = 1, 4
    call assert_true(abs(m(i) - 15.0_dp * real(counts(i),dp) / 100.0_dp) < 0.035_dp, 'hypergeometric means')
  end do
  call assert_true(all(sum(hg, dim=2) == 15), 'hypergeometric row sums')
  write(*,'(a)') 'test_discrete_mixtures: PASS'
end program test_discrete_mixtures
