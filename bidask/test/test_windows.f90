! SPDX-License-Identifier: MIT
program test_windows
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  use bidask, only: dp, ohlc_data, spread_result, spread_series_result, &
    edge_estimate, edge_rolling, edge_expanding, spread, spread_expanding, spread_endpoints
  implicit none
  integer, parameter :: n = 80
  type(ohlc_data) :: data
  type(spread_result) :: one
  type(spread_series_result) :: series
  real(dp), allocatable :: values(:)
  integer :: widths(n)
  character(len=16) :: methods(4)
  integer :: i

  allocate(data%open(n), data%high(n), data%low(n), data%close(n))
  call make_prices(data%open, data%high, data%low, data%close)
  values = edge_rolling(data%open, data%high, data%low, data%close, 20, na_rm=.true.)
  call assert_true(ieee_is_nan(values(19)))
  call assert_close(values(n), 0.0007430099922199551_dp, 3.0e-13_dp)

  values = edge_expanding(data%open, data%high, data%low, data%close)
  call assert_close(values(40), 0.0009897266336830544_dp, 3.0e-13_dp)
  call assert_close(values(n), edge_estimate(data%open, data%high, data%low, data%close), 1.0e-15_dp)

  widths = [(i, i = 1, n)]
  values = edge_rolling(data%open, data%high, data%low, data%close, widths, na_rm=.true.)
  call assert_close(values(n), edge_estimate(data%open, data%high, data%low, data%close), 1.0e-15_dp)

  methods = ['EDGE            ', 'AR              ', 'CS              ', 'OHLC.CHLO       ']
  one = spread(data, methods, signed=.true., na_rm=.true.)
  call assert_true(one%ok)
  call assert_close(one%value(1), 0.0016767995272337732_dp, 3.0e-13_dp)
  call assert_close(one%value(2), -0.0035365413824787898_dp, 3.0e-13_dp)

  series = spread(data, 20, methods, signed=.true., na_rm=.true.)
  call assert_close(series%value(n, 1), values_from_last20(data), 1.0e-15_dp)

  series = spread_expanding(data, methods, na_rm=.true.)
  call assert_close(series%value(n, 1), one%value(1), 3.0e-13_dp)

  series = spread_endpoints(data, [1, 40, 80], methods, na_rm=.true.)
  call assert_close(series%value(40, 1), edge_estimate(data%open(1:40), data%high(1:40), &
    data%low(1:40), data%close(1:40)), 1.0e-15_dp)
  call assert_close(series%value(80, 1), edge_estimate(data%open(40:80), data%high(40:80), &
    data%low(40:80), data%close(40:80)), 1.0e-15_dp)
  print '(a)', 'test_windows: PASS'

contains

  real(dp) function values_from_last20(x) result(value)
    type(ohlc_data), intent(in) :: x
    value = edge_estimate(x%open(61:80), x%high(61:80), x%low(61:80), x%close(61:80), &
      signed=.true., na_rm=.true.)
  end function values_from_last20

  subroutine make_prices(open, high, low, close)
    real(dp), intent(out) :: open(:), high(:), low(:), close(:)
    real(dp) :: previous
    integer :: j
    previous = 100.0_dp
    do j = 1, size(open)
      close(j) = 100.0_dp * exp(0.0007_dp * real(j, dp) + &
        0.012_dp * sin(0.37_dp * real(j, dp)) + 0.003_dp * cos(0.11_dp * real(j, dp)))
      if (mod(j - 1, 2) == 0) then
        open(j) = previous * exp(0.0025_dp * cos(0.29_dp * real(j, dp))) * 1.0025_dp
      else
        open(j) = previous * exp(0.0025_dp * cos(0.29_dp * real(j, dp))) * 0.9975_dp
      end if
      high(j) = max(open(j), close(j)) * exp(0.004_dp + 0.0015_dp * abs(sin(0.51_dp * real(j, dp))))
      low(j) = min(open(j), close(j)) * exp(-0.0045_dp - 0.001_dp * abs(cos(0.43_dp * real(j, dp))))
      previous = close(j)
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

end program test_windows
