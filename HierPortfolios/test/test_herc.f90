! SPDX-License-Identifier: GPL-2.0-only
program test_herc
  use hierportfolios, only: dp, portfolio_result, HERC_Portfolio
  implicit none

  real(dp) :: covar(6, 6), ratio12, ratio23
  type(portfolio_result) :: result

  call make_block_covariance(covar)
  call HERC_Portfolio(covar, result, clusters=2)
  call check(result%ok(), 'HERC status')
  call check(result%n_clusters == 2, 'HERC cluster count')
  call check(abs(sum(result%weights) - 1.0_dp) < 1.0e-12_dp, 'HERC sum')
  call check(all(result%weights > 0.0_dp), 'HERC positivity')
  call check(sum(result%weights(1:3)) > sum(result%weights(4:6)), &
    'HERC inverse cluster-risk allocation')

  ratio12 = result%weights(1) / result%weights(2)
  ratio23 = result%weights(2) / result%weights(3)
  call check(abs(ratio12 - covar(2, 2) / covar(1, 1)) < 1.0e-10_dp, &
    'HERC inverse variance ratio 1-2')
  call check(abs(ratio23 - covar(3, 3) / covar(2, 2)) < 1.0e-10_dp, &
    'HERC inverse variance ratio 2-3')

  print '(a)', 'test_herc: PASS'

contains

  subroutine make_block_covariance(c)
    real(dp), intent(out) :: c(:, :)
    integer :: i, j

    c = 0.0_dp
    do i = 1, 6
      c(i, i) = 0.02_dp + 0.005_dp * real(i, dp)
    end do
    do j = 1, 3
      do i = 1, 3
        if (i /= j) c(i, j) = 0.60_dp * sqrt(c(i, i) * c(j, j))
      end do
    end do
    do j = 4, 6
      do i = 4, 6
        if (i /= j) c(i, j) = 0.50_dp * sqrt(c(i, i) * c(j, j))
      end do
    end do
    do j = 4, 6
      do i = 1, 3
        c(i, j) = 0.05_dp * sqrt(c(i, i) * c(j, j))
        c(j, i) = c(i, j)
      end do
    end do
  end subroutine make_block_covariance

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write (*, '(a)') 'FAIL: '//label
      error stop 1
    end if
  end subroutine check

end program test_herc
