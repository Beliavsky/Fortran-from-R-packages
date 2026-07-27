! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program interpolation_and_rehedging
   use opthedging, only : dp, hedging_iid, hedging_result
   implicit none

   real(dp) :: log_returns(5)
   real(dp) :: portfolio
   real(dp) :: spot
   type(hedging_result) :: result

   log_returns = log([0.90_dp, 0.96_dp, 1.00_dp, 1.05_dp, 1.12_dp])
   result = hedging_iid(log_returns, 1.0_dp, 100.0_dp, 0.02_dp, .false., &
      3, 301, 50.0_dp, 170.0_dp)
   if (.not. result%ok) error stop result%message

   spot = 103.25_dp
   portfolio = result%option_value_at(1, spot)
   print '(a,f10.5)', 'interpolated call value: ', portfolio
   print '(a,f10.5)', 'shares at time zero:      ', &
      result%shares_at(1, spot, portfolio)
end program interpolation_and_rehedging
