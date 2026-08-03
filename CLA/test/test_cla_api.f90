! SPDX-License-Identifier: GPL-3.0-or-later
program test_cla_api
   use kind_mod, only: dp
   use cla_api, only: CLA, MS, findSig, findMu, cla_result_t, cla_path_query_t
   implicit none
   real(dp) :: mu(2), covar(2,2), sig(2), means(2)
   type(cla_result_t) :: result
   type(cla_path_query_t) :: a, b
   mu=[0.08_dp,0.03_dp]
   covar=reshape([0.04_dp,0.005_dp,0.005_dp,0.01_dp],[2,2])
   result=CLA(mu,covar,0.0_dp,1.0_dp)
   if(result%info/=0 .or. result%n_turning<2)error stop 'CLA compatibility failed'
   call MS(result%weights,mu,covar,sig,means)
   if(maxval(abs(sig-result%sigma))>1.0e-12_dp)error stop 'MS compatibility failed'
   a=findSig([sum(result%mu)/real(result%n_turning,dp)],result,covar)
   if(a%info/=0)error stop 'findSig compatibility failed'
   b=findMu(a%value,result,covar,tolerance=1.0e-10_dp)
   if(b%info/=0)error stop 'findMu compatibility failed'
   print '(a)', 'test_cla_api: PASS'
end program test_cla_api
