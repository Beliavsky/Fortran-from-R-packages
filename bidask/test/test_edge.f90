! SPDX-License-Identifier: MIT
program test_edge
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
  use bidask, only: dp, edge_estimate
  implicit none
  integer, parameter :: n = 80
  real(dp) :: o(n), h(n), l(n), c(n), value

  call make_prices(o, h, l, c)
  value = edge_estimate(o, h, l, c)
  call assert_close(value, 0.0016767995272337732_dp, 2.0e-13_dp)
  value = edge_estimate(o, h, l, c, signed=.true.)
  call assert_close(value, 0.0016767995272337732_dp, 2.0e-13_dp)

  value = edge_estimate([18.21_dp, 17.61_dp, 17.61_dp], &
    [18.21_dp, 17.61_dp, 17.61_dp], [17.61_dp, 17.61_dp, 17.61_dp], &
    [17.61_dp, 17.61_dp, 17.61_dp])
  call assert_true(ieee_is_nan(value))

  c(20) = ieee_value(0.0_dp, ieee_quiet_nan)
  call assert_true(ieee_is_nan(edge_estimate(o, h, l, c, na_rm=.false.)))
  call assert_true(.not. ieee_is_nan(edge_estimate(o, h, l, c, na_rm=.true.)))
  print '(a)', 'test_edge: PASS'

contains

  subroutine make_prices(open, high, low, close)
    real(dp), intent(out) :: open(:), high(:), low(:), close(:)
    real(dp) :: previous
    integer :: i
    previous = 100.0_dp
    do i = 1, size(open)
      close(i) = 100.0_dp * exp(0.0007_dp * real(i, dp) + &
        0.012_dp * sin(0.37_dp * real(i, dp)) + &
        0.003_dp * cos(0.11_dp * real(i, dp)))
      if (mod(i - 1, 2) == 0) then
        open(i) = previous * exp(0.0025_dp * cos(0.29_dp * real(i, dp))) * 1.0025_dp
      else
        open(i) = previous * exp(0.0025_dp * cos(0.29_dp * real(i, dp))) * 0.9975_dp
      end if
      high(i) = max(open(i), close(i)) * exp(0.004_dp + &
        0.0015_dp * abs(sin(0.51_dp * real(i, dp))))
      low(i) = min(open(i), close(i)) * exp(-0.0045_dp - &
        0.001_dp * abs(cos(0.43_dp * real(i, dp))))
      previous = close(i)
    end do
  end subroutine make_prices

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print '(a,3es25.16)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_edge
