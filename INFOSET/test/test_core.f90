! SPDX-License-Identifier: GPL-2.0-or-later
program test_core
  use infoset, only : dp, g_ret, create_overlapping_windows
  use infoset, only : window_collection, infoset_success
  implicit none
  real(dp) :: prices(3), matrix_prices(5,2)
  real(dp), allocatable :: r(:), rm(:,:)
  type(window_collection) :: windows
  integer :: i

  prices = [100.0_dp, 110.0_dp, 99.0_dp]
  r = g_ret(prices)
  call check(size(r) == 2, 'gross return size')
  call check(abs(r(1) - 0.9_dp) < 1.0e-12_dp, 'gross return low')
  call check(abs(r(2) - 1.1_dp) < 1.0e-12_dp, 'gross return high')

  do i = 1, 5
    matrix_prices(i,1) = real(i,dp)
    matrix_prices(i,2) = 2.0_dp*real(i,dp)
  end do
  rm = g_ret(matrix_prices)
  call check(all(shape(rm) == [4,2]), 'matrix gross return shape')
  call check(maxval(abs(rm(:,1)-rm(:,2))) < 1.0e-12_dp, 'matrix gross returns')

  windows = create_overlapping_windows(matrix_prices, 3, 1)
  call check(windows%status == infoset_success, 'window status')
  call check(all(shape(windows%values) == [3,2,3]), 'window shape')
  call check(abs(windows%values(1,1,3)-3.0_dp) < 1.0e-12_dp, 'window start')

  windows = create_overlapping_windows(matrix_prices, 8, 1)
  call check(windows%status /= infoset_success, 'invalid window')
  print '(a)', 'test_core: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//message
      error stop 1
    end if
  end subroutine check
end program test_core
