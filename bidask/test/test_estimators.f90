! SPDX-License-Identifier: MIT
program test_estimators
  use bidask, only: dp, ar_estimate, cs_estimate, roll_estimate, ohlc_estimate, estimate_method
  implicit none
  integer, parameter :: n = 80
  real(dp) :: o(n), h(n), l(n), c(n)

  call make_prices(o, h, l, c)
  call assert_close(ar_estimate(h, l, c, signed=.true., na_rm=.true.), &
    -0.0035365413824787898_dp, 3.0e-13_dp)
  call assert_close(ar_estimate(h, l, c, two_period=.true., na_rm=.true.), &
    0.0005610368063406199_dp, 3.0e-13_dp)
  call assert_close(cs_estimate(h, l, c, signed=.true., na_rm=.true.), &
    0.005276368508308839_dp, 3.0e-13_dp)
  call assert_close(cs_estimate(h, l, c, two_period=.true., na_rm=.true.), &
    0.00554050413846179_dp, 3.0e-13_dp)
  call assert_close(roll_estimate(c, signed=.true., na_rm=.true.), &
    -0.006081489001682237_dp, 3.0e-13_dp)
  call assert_close(ohlc_estimate(o, h, l, c, 'OHL', signed=.true., na_rm=.true.), &
    0.0035475807562462495_dp, 3.0e-13_dp)
  call assert_close(ohlc_estimate(o, h, l, c, 'OHLC', signed=.true., na_rm=.true.), &
    0.004304587537858832_dp, 3.0e-13_dp)
  call assert_close(ohlc_estimate(o, h, l, c, 'CHL', signed=.true., na_rm=.true.), &
    -0.003528934222275473_dp, 3.0e-13_dp)
  call assert_close(ohlc_estimate(o, h, l, c, 'CHLO', signed=.true., na_rm=.true.), &
    -0.002551319677372114_dp, 3.0e-13_dp)
  call assert_close(ohlc_estimate(o, h, l, c, 'OHLC.CHLO', signed=.true., na_rm=.true.), &
    0.0024515547898163945_dp, 3.0e-13_dp)
  call assert_close(ohlc_estimate(o, h, l, c, 'OHL.CHL', signed=.true., na_rm=.true.), &
    0.0002568584015967925_dp, 3.0e-13_dp)
  call assert_close(estimate_method(o, h, l, c, 'roll', signed=.true., na_rm=.true.), &
    -0.006081489001682237_dp, 3.0e-13_dp)
  print '(a)', 'test_estimators: PASS'

contains

  subroutine make_prices(open, high, low, close)
    real(dp), intent(out) :: open(:), high(:), low(:), close(:)
    real(dp) :: previous
    integer :: i
    previous = 100.0_dp
    do i = 1, size(open)
      close(i) = 100.0_dp * exp(0.0007_dp * real(i, dp) + &
        0.012_dp * sin(0.37_dp * real(i, dp)) + 0.003_dp * cos(0.11_dp * real(i, dp)))
      if (mod(i - 1, 2) == 0) then
        open(i) = previous * exp(0.0025_dp * cos(0.29_dp * real(i, dp))) * 1.0025_dp
      else
        open(i) = previous * exp(0.0025_dp * cos(0.29_dp * real(i, dp))) * 0.9975_dp
      end if
      high(i) = max(open(i), close(i)) * exp(0.004_dp + 0.0015_dp * abs(sin(0.51_dp * real(i, dp))))
      low(i) = min(open(i), close(i)) * exp(-0.0045_dp - 0.001_dp * abs(cos(0.43_dp * real(i, dp))))
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

end program test_estimators
