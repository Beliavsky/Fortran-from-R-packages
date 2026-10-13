! SPDX-License-Identifier: GPL-2.0-or-later
program technical_indicators
   use ttr, only: dp, ema, rsi, bbands, volatility
   implicit none
   real(dp) :: price(60), average(60), strength(60), bands(60, 4), ohlc(60, 4), risk(60)
   integer :: i
   do i = 1, size(price)
      price(i) = 100.0_dp + 0.2_dp * real(i, dp) + 2.0_dp * sin(real(i, dp) / 3.0_dp)
   end do
   average = ema(price, 10)
   strength = rsi(price, 14)
   bands = bbands(reshape(price, [60, 1]), 20)
   ohlc(:, 1) = price - 0.2_dp
   ohlc(:, 2) = price + 1.0_dp
   ohlc(:, 3) = price - 1.0_dp
   ohlc(:, 4) = price
   risk = volatility(ohlc, 10, calc = 'garman.klass')
   print '(a,f12.6)', 'Final 10-period EMA: ', average(60)
   print '(a,f12.6)', 'Final 14-period RSI: ', strength(60)
   print '(a,3f12.6)', '20-period lower/middle/upper bands: ', bands(60, :3)
   print '(a,f12.6)', 'Annualized Garman-Klass volatility: ', risk(60)
end program technical_indicators
