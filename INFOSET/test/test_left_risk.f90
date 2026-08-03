! SPDX-License-Identifier: GPL-2.0-or-later
program test_left_risk
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use infoset, only : dp, lr_cp, left_risk_result
  use infoset, only : infoset_success, infoset_no_split
  implicit none
  integer, parameter :: n = 180, p = 3
  real(dp) :: prices(n,p), return_value
  type(left_risk_result) :: risk
  integer :: i, j

  prices(1,:) = [100.0_dp, 80.0_dp, 120.0_dp]
  do i = 2, n
    do j = 1, p
      return_value = 0.0004_dp + 0.003_dp*sin(0.13_dp*real(i*j,dp))
      if (mod(i + 7*j, 41) == 0) return_value = return_value - 0.045_dp
      prices(i,j) = prices(i-1,j)*exp(return_value)
    end do
  end do
  call lr_cp(prices, 80, 20, risk)
  call check(risk%status == infoset_success .or. &
    risk%status == infoset_no_split, 'left risk status')
  call check(all(shape(risk%values) == [p,6]), 'left risk shape')
  call check(all(ieee_is_finite(risk%values)), 'left risk finite')
  call check(all(risk%first_change_point > 0.0_dp), 'change points positive')
  print '(a)', 'test_left_risk: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_left_risk
