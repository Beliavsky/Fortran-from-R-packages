! SPDX-License-Identifier: GPL-3.0-or-later
program cla_demo
   use kind_mod, only: dp
   use cla, only: cla_result_t, critical_line
   implicit none
   real(dp) :: mu(3), covar(3,3), lower(3), upper(3)
   type(cla_result_t) :: result
   integer :: j

   mu = [0.0408_dp,0.102_dp,-0.023_dp]
   covar = reshape([0.00648_dp,0.00792_dp,0.00473_dp, &
                    0.00792_dp,0.0334_dp,0.0121_dp, &
                    0.00473_dp,0.0121_dp,0.0793_dp],[3,3])
   lower = 0.0_dp
   upper = 1.0_dp
   result = critical_line(mu,covar,lower,upper)
   if(result%info/=0)error stop 'CLA failed'
   write(*,'(a)') ' turning       lambda         sigma            mu       weights'
   do j=1,result%n_turning
      write(*,'(i8,3f14.8,3f11.6)')j,result%lambdas(j),result%sigma(j), &
         result%mu(j),result%weights(:,j)
   end do
end program cla_demo
