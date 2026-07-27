! SPDX-License-Identifier: MIT
program bidask_demo
  use iso_fortran_env, only: int64
  use bidask, only: dp, ohlc_data, spread_result, sim, spread
  implicit none
  type(ohlc_data) :: prices
  type(spread_result) :: estimates
  character(len=16) :: methods(7)
  integer :: i

  methods = ['EDGE            ', 'AR              ', 'AR2             ', &
    'CS              ', 'CS2             ', 'ROLL            ', 'OHLC.CHLO       ']
  prices = sim(1000, 100, probability=1.0_dp, spread=0.01_dp, &
    volatility=0.03_dp, seed=20260726_int64)
  estimates = spread(prices, methods, na_rm=.true.)

  print '(a)', 'Bid-ask spread estimates'
  do i = 1, size(estimates%method)
    print '(2x,a16,2x,f12.8)', trim(estimates%method(i)), estimates%value(i)
  end do
end program bidask_demo
