! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program opthedging_demo
   use opthedging, only : dp, hedging_iid, hedging_result, random_normal, seed_random
   implicit none

   integer, parameter :: n_periods = 5
   integer, parameter :: n_grid = 1001
   integer, parameter :: n_returns = 20000
   real(dp), parameter :: maturity = 1.0_dp
   real(dp), parameter :: strike = 100.0_dp
   real(dp), parameter :: rate = 0.05_dp
   real(dp), parameter :: sigma = 0.20_dp
   real(dp), parameter :: mu = 0.09_dp
   real(dp), parameter :: spot = 100.0_dp
   real(dp), allocatable :: returns(:)
   real(dp) :: period
   real(dp) :: period_rate
   real(dp) :: period_sigma
   type(hedging_result) :: fit

   allocate(returns(n_returns))
   call seed_random(20260725)
   call random_normal(returns)
   period = maturity / real(n_periods, dp)
   period_rate = rate * period
   period_sigma = sigma * sqrt(period)
   returns = mu * period - 0.5_dp * period_sigma**2 - period_rate + &
      period_sigma * returns

   fit = hedging_iid(returns, maturity, strike, rate, .true., n_periods, &
      n_grid, 40.0_dp, 160.0_dp)
   if (.not. fit%ok) error stop fit%message

   print '(a,f12.6)', 'discounted option value: ', fit%option_value_at(1, spot)
   print '(a,f12.6)', 'initial shares:          ', fit%initial_hedge_at(spot)
   print '(a,f12.6)', 'rho:                     ', fit%rho
end program opthedging_demo
