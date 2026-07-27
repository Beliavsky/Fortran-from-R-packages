! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program parma_demo
   use parma
   implicit none
   type(parma_spec) :: spec
   type(parma_port) :: solution
   type(parma_options) :: options
   real(dp) :: data(8,3),lb(3),ub(3)
   integer :: info

   data = reshape([ &
      0.010_dp,-0.012_dp,0.018_dp,0.004_dp,-0.006_dp,0.013_dp,0.002_dp,0.008_dp, &
      0.006_dp,0.004_dp,-0.003_dp,0.007_dp,0.005_dp,-0.002_dp,0.009_dp,0.003_dp, &
      0.015_dp,-0.020_dp,0.025_dp,-0.008_dp,0.011_dp,0.018_dp,-0.010_dp,0.014_dp],[8,3])
   lb = 0.0_dp
   ub = 0.8_dp
   call parmaspec(spec,data=data,risk=risk_cvar,objective=solve_min_risk,lb=lb,ub=ub, &
      alpha=0.25_dp,target=0.004_dp,info=info)
   options%max_iter = 1200
   options%seed = 24680
   call parmasolve(spec,solution,options)

   print '(a,i0)', 'status: ',solution%status
   print '(a,*(f10.6,1x))', 'weights: ',solution%weights
   print '(a,f12.7)', 'reward: ',solution%reward
   print '(a,f12.7)', 'CVaR:  ',solution%risk
end program parma_demo
