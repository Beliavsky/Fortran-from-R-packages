! SPDX-License-Identifier: GPL-2.0-only
program test_hrp
  use hierportfolios, only: dp, portfolio_result, HRP_Portfolio, hp_invalid_argument
  implicit none

  real(dp) :: covar(6, 6), bad(2, 3)
  type(portfolio_result) :: result

  call make_block_covariance(covar)
  call HRP_Portfolio(covar, result, 'single')
  call check(result%ok(), 'HRP status')
  call check(abs(sum(result%weights) - 1.0_dp) < 1.0e-12_dp, 'HRP sum')
  call check(all(result%weights > 0.0_dp), 'HRP positivity')
  call check(size(result%order) == 6, 'HRP order')
  call check(sum(result%weights(1:3)) > sum(result%weights(4:6)), &
    'HRP lower-risk block allocation')

  call HRP_Portfolio(covar, result, 'invalid')
  call check(result%status == hp_invalid_argument, 'invalid linkage')
  bad = 0.0_dp
  call HRP_Portfolio(bad, result)
  call check(result%status == hp_invalid_argument, 'nonsquare covariance')

  print '(a)', 'test_hrp: PASS'

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

end program test_hrp
