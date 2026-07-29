! SPDX-License-Identifier: GPL-3.0-or-later
program series_tools
  use frapo
  implicit none

  real(dp) :: prices(6)
  real(dp), allocatable :: returns(:), smooth(:), hp(:)

  prices = [100.0_dp, 101.0_dp, 100.5_dp, 103.0_dp, 104.0_dp, 105.5_dp]
  returns = returnseries(prices, method=returns_discrete, &
                         percentage=.false., trim=.true.)
  smooth = trdes(returns, lambda=0.25_dp)
  hp = trdhp(prices, lambda=1600.0_dp)

  write(*, '(a,*(f11.6,1x))') 'Discrete returns: ', returns
  write(*, '(a,*(f11.6,1x))') 'Exponential trend: ', smooth
  write(*, '(a,*(f11.6,1x))') 'HP price trend: ', hp
end program series_tools
