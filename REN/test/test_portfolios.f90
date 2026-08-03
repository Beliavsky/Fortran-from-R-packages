! SPDX-License-Identifier: AGPL-3.0-or-later
program test_portfolios
  use ren, only : dp, po_cols, po_jm, po_avg, po_cov_shrink, po_tzt, po_sw
  implicit none
  real(dp) :: x(30, 4), y(30)
  real(dp), allocatable :: w(:)
  integer :: i, j, status
  do j = 1, 4
    do i = 1, 30
      x(i, j) = 0.03_dp * real(i, dp) + sin(0.17_dp * real(i * j, dp)) + 0.2_dp * real(j, dp)
    end do
  end do
  y = 0.0_dp
  call po_cols(y, x, w, status)
  call assert_close(sum(w), 1.0_dp, 1.0e-8_dp, 'po_cols sum')
  call po_jm(x, w, status)
  call assert_close(sum(w), 1.0_dp, 1.0e-8_dp, 'po_jm sum')
  if (minval(w) < -1.0e-10_dp) error stop 'po_jm negative weight'
  call po_avg(y, x, 'RIDGE', w, status, seed=11)
  call assert_close(sum(w), 1.0_dp, 1.0e-7_dp, 'po_avg sum')
  call po_cov_shrink(y, x, w, status)
  call assert_close(sum(w), 1.0_dp, 1.0e-8_dp, 'po_cov_shrink sum')
  call po_tzt(x, 3.0_dp, w, status)
  call assert_close(sum(w), 1.0_dp, 1.0e-7_dp, 'po_tzt sum')
  call po_sw(x, 3, 20, w, status, seed=22)
  call assert_close(sum(w), 1.0_dp, 1.0e-8_dp, 'po_sw sum')
  print '(a)', 'test_portfolios: PASS'
contains
  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      print '(a,2(1x,es16.8))', trim(label), actual, expected
      error stop 1
    end if
  end subroutine assert_close
end program test_portfolios
