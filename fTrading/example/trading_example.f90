! SPDX-License-Identifier: LGPL-2.0-or-later
program trading_example
   use ftrading, only: dp, ema_ta, max_draw_down, sharpe_ratio, bollinger_ta, roll_fun
   implicit none
   real(dp) :: levels(40), average(40), bands(40, 3), maximum
   real(dp), allocatable :: rolling(:)
   integer, allocatable :: from(:), to(:)
   integer :: i
   do i = 1, size(levels)
      levels(i) = 100.0_dp + 0.1_dp * real(i, dp) + 2.0_dp * sin(real(i, dp) / 3.0_dp)
   end do
   average = ema_ta(levels, 5.0_dp)
   bands = bollinger_ta(levels, 5)
   rolling = roll_fun(levels, 5, 'mean')
   call max_draw_down(levels, maximum, from, to)
   print '(a,f10.5)', 'Last startup-seeded EMA: ', average(40)
   print '(a,f10.5)', 'Last complete rolling mean: ', rolling(size(rolling))
   print '(a,3f10.5)', 'Last upper / original / lower band: ', bands(40, :)
   print '(a,f10.5)', 'Maximum absolute drawdown: ', maximum
   print '(a,2i5)', 'Last peak and trough: ', from(size(from)), to(size(to))
   print '(a,f10.5)', 'Annualized Sharpe of level differences: ', sharpe_ratio(levels)
end program trading_example
