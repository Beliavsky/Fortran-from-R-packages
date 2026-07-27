! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
program efficient_frontier
   use parma
   implicit none
   type(parma_spec) :: spec
   type(parma_port), allocatable :: frontier(:)
   type(parma_options) :: options
   real(dp) :: data(6,2),targets(3),lb(2),ub(2)
   integer :: i,info

   data = reshape([0.01_dp,-0.01_dp,0.02_dp,0.00_dp,0.015_dp,-0.005_dp, &
                   0.02_dp,0.01_dp,-0.02_dp,0.015_dp,0.005_dp,0.025_dp],[6,2])
   lb = 0.0_dp
   ub = 1.0_dp
   targets = [0.006_dp,0.008_dp,0.009_dp]
   call parmaspec(spec,data=data,risk=risk_ev,objective=solve_min_risk,lb=lb,ub=ub,info=info)
   options%max_iter = 1500
   call parmafrontier(spec,targets,frontier,options)
   print '(a)', ' target       reward         risk       weights'
   do i = 1,size(frontier)
      print '(3f13.7,2f11.6)',targets(i),frontier(i)%reward,frontier(i)%risk,frontier(i)%weights
   end do
end program efficient_frontier
