! SPDX-License-Identifier: MIT
program basic_estimators
  use bidask, only: dp, edge, ar_estimate, cs_estimate, roll_estimate, ohlc_estimate
  implicit none
  integer, parameter :: n = 60
  real(dp) :: o(n), h(n), l(n), c(n), previous
  integer :: i

  previous = 100.0_dp
  do i = 1, n
    c(i) = 100.0_dp * exp(0.0005_dp * real(i, dp) + 0.01_dp * sin(0.3_dp * real(i, dp)))
    o(i) = previous * exp(0.002_dp * cos(0.2_dp * real(i, dp)))
    h(i) = max(o(i), c(i)) * exp(0.005_dp)
    l(i) = min(o(i), c(i)) * exp(-0.005_dp)
    previous = c(i)
  end do

  print '(a,f12.8)', 'EDGE:       ', edge(o, h, l, c)
  print '(a,f12.8)', 'AR:         ', ar_estimate(h, l, c, na_rm=.true.)
  print '(a,f12.8)', 'CS:         ', cs_estimate(h, l, c, na_rm=.true.)
  print '(a,f12.8)', 'ROLL:       ', roll_estimate(c, na_rm=.true.)
  print '(a,f12.8)', 'OHLC.CHLO:  ', ohlc_estimate(o, h, l, c, 'OHLC.CHLO', na_rm=.true.)
end program basic_estimators
