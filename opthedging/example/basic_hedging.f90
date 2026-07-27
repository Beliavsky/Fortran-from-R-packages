! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
program basic_hedging
   use opthedging, only : dp, hedging_iid, hedging_result
   implicit none

   real(dp) :: log_returns(6)
   type(hedging_result) :: result

   log_returns = log([0.85_dp, 0.93_dp, 0.99_dp, 1.02_dp, 1.08_dp, 1.18_dp])
   result = hedging_iid(log_returns, 0.5_dp, 100.0_dp, 0.03_dp, .true., &
      4, 401, 50.0_dp, 150.0_dp)
   if (.not. result%ok) error stop result%message

   print '(a,f10.5)', 'put value at S=100: ', result%option_value_at(1, 100.0_dp)
   print '(a,f10.5)', 'initial hedge:       ', result%initial_hedge_at(100.0_dp)
end program basic_hedging
