! SPDX-License-Identifier: MIT
program rolling_and_simulation
  use iso_fortran_env, only: int64
  use bidask, only: dp, ohlc_data, sim, edge_rolling, edge_expanding
  implicit none
  type(ohlc_data) :: prices
  real(dp), allocatable :: rolling(:), expanding(:)
  integer :: n

  prices = sim(250, 100, probability=0.9_dp, spread=0.01_dp, &
    volatility=0.025_dp, overnight=0.005_dp, seed=123456_int64)
  rolling = edge_rolling(prices%open, prices%high, prices%low, prices%close, &
    21, na_rm=.true.)
  expanding = edge_expanding(prices%open, prices%high, prices%low, prices%close)
  n = prices%size()
  print '(a,f12.8)', 'Last 21-period EDGE estimate: ', rolling(n)
  print '(a,f12.8)', 'Full expanding EDGE estimate: ', expanding(n)
end program rolling_and_simulation
