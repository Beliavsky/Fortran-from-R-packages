! SPDX-License-Identifier: AGPL-3.0-or-later
program basic_portfolios
  use ren, only : dp, po_jm, po_avg
  implicit none
  real(dp) :: returns(40, 4), zero(40)
  real(dp), allocatable :: weights(:)
  integer :: i, j, status
  do j = 1, 4
    do i = 1, 40
      returns(i, j) = sin(0.1_dp * real(i * j, dp)) + 0.03_dp * real(j, dp)
    end do
  end do
  zero = 0.0_dp
  call po_jm(returns, weights, status)
  print '(a,*(1x,f9.5))', 'long-only minimum variance:', weights
  call po_avg(zero, returns, 'LASSO', weights, status, seed=123)
  print '(a,*(1x,f9.5))', 'averaged lasso:', weights
end program basic_portfolios
