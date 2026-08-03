! SPDX-License-Identifier: GPL-3.0-or-later
program frontier_query
   use kind_mod, only: dp
   use cla, only: cla_result_t, cla_path_query_t, critical_line, find_sigma, find_mu
   implicit none
   real(dp) :: mu(3), covar(3,3)
   type(cla_result_t) :: frontier
   type(cla_path_query_t) :: by_mu, by_sigma
   mu=[0.0408_dp,0.102_dp,-0.023_dp]
   covar=reshape([0.00648_dp,0.00792_dp,0.00473_dp,0.00792_dp,0.0334_dp, &
      0.0121_dp,0.00473_dp,0.0121_dp,0.0793_dp],[3,3])
   frontier=critical_line(mu,covar,[0.0_dp,0.0_dp,0.0_dp],[1.0_dp,1.0_dp,1.0_dp])
   by_mu=find_sigma([0.04_dp],frontier,covar)
   by_sigma=find_mu(by_mu%value,frontier,covar,tolerance=1.0e-11_dp)
   write(*,'(a,f12.8,a,f12.8)')'At mean ',by_sigma%value(1),', sigma = ',by_mu%value(1)
   write(*,'(a,3f12.8)')'Weights: ',by_mu%weights(:,1)
end program frontier_query
