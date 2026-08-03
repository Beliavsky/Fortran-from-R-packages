! SPDX-License-Identifier: GPL-3.0-or-later
program test_cla_bounds
   use kind_mod, only: dp
   use cla, only: cla_result_t, critical_line, cla_infeasible_bounds
   implicit none
   real(dp) :: mu(4), covar(4,4), lower(4), upper(4)
   type(cla_result_t) :: result
   integer :: i

   mu = [0.09_dp,0.07_dp,0.05_dp,0.02_dp]
   covar = reshape([0.040_dp,0.006_dp,0.004_dp,0.002_dp, &
                    0.006_dp,0.025_dp,0.005_dp,0.003_dp, &
                    0.004_dp,0.005_dp,0.018_dp,0.002_dp, &
                    0.002_dp,0.003_dp,0.002_dp,0.010_dp],[4,4])
   lower = [0.05_dp,0.0_dp,0.0_dp,0.10_dp]
   upper = [0.50_dp,0.60_dp,0.70_dp,0.80_dp]
   result = critical_line(mu,covar,lower,upper)
   if(result%info/=0)error stop 'bounded CLA failed'
   do i=1,result%n_turning
      if(abs(sum(result%weights(:,i))-1.0_dp)>1.0e-11_dp)error stop 'budget failed'
      if(any(result%weights(:,i)<lower-1.0e-10_dp))error stop 'lower bound failed'
      if(any(result%weights(:,i)>upper+1.0e-10_dp))error stop 'upper bound failed'
   end do
   lower = 0.4_dp
   upper = 0.5_dp
   result = critical_line(mu,covar,lower,upper)
   if(result%info/=cla_infeasible_bounds)error stop 'infeasible bounds not detected'
   print '(a)', 'test_cla_bounds: PASS'
end program test_cla_bounds
