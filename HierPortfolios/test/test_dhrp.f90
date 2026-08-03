! SPDX-License-Identifier: GPL-2.0-only
program test_dhrp
  use hierportfolios, only: dp, portfolio_result, DHRP_Portfolio, hp_invalid_argument
  implicit none

  real(dp) :: covar(6, 6), lower(6), upper(6), small(2, 2)
  type(portfolio_result) :: middle, flexible, invalid, two_asset

  call make_block_covariance(covar)
  lower = 0.05_dp
  upper = 0.22_dp
  call DHRP_Portfolio(covar, middle, tau=0.0_dp, lb=lower, ub=upper)
  call check(middle%ok(), 'DHRP status')
  call check(abs(sum(middle%weights) - 1.0_dp) < 1.0e-11_dp, 'DHRP sum')
  call check(all(middle%weights >= lower - 1.0e-12_dp), 'DHRP lower bounds')
  call check(all(middle%weights <= upper + 1.0e-12_dp), 'DHRP upper bounds')

  call DHRP_Portfolio(covar, flexible, tau=1.0_dp, lb=lower, ub=upper)
  call check(flexible%ok(), 'DHRP flexible status')
  call check(maxval(abs(flexible%weights - middle%weights)) > 1.0e-5_dp, &
    'DHRP tau changes splits')

  small = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], [2, 2])
  call DHRP_Portfolio(small, two_asset, tau=0.5_dp)
  call check(two_asset%ok(), 'DHRP two-asset projection')
  call check(abs(sum(two_asset%weights) - 1.0_dp) < 1.0e-12_dp, &
    'DHRP two-asset sum')

  lower = 0.20_dp
  upper = 0.30_dp
  call DHRP_Portfolio(covar, invalid, lb=lower, ub=upper)
  call check(invalid%status == hp_invalid_argument, 'DHRP infeasible bounds')

  print '(a)', 'test_dhrp: PASS'

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

end program test_dhrp
