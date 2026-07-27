! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program risk_measures
   use parma
   implicit none
   real(dp) :: data(5,2),weights(2)

   data = reshape([0.01_dp,-0.03_dp,0.02_dp,0.00_dp,0.015_dp, &
                   0.005_dp,0.01_dp,-0.01_dp,0.02_dp,-0.005_dp],[5,2])
   weights = [0.6_dp,0.4_dp]
   print '(a,f12.7)', 'MAD:     ',mad_risk(weights,data)
   print '(a,f12.7)', 'Variance:',variance_risk(weights,data)
   print '(a,f12.7)', 'CVaR:    ',cvar_risk(weights,data,0.20_dp)
   print '(a,f12.7)', 'CDaR:    ',cdar_risk(weights,data,0.20_dp)
   print '(a,f12.7)', 'LPM(1):  ',lpm_risk(weights,data,0.0_dp,1.0_dp)
end program risk_measures
