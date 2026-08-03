! SPDX-License-Identifier: Artistic-2.0
program volatility_analysis
   use mts
   implicit none
   integer,parameter::n=300,k=2
   real(dp)::returns(n,k),draw(k),lambda,loglik
   real(dp),allocatable::covariance(:,:,:),correlation(:,:,:)
   type(dcc_model)::dcc
   integer::t,status

   call set_random_seed(3303)
   do t=1,n
      call random_multivariate_normal([0.0_dp,0.0_dp],reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[k,k]),draw)
      returns(t,:)=draw
   end do
   call fit_ewma_lambda(returns,lambda,covariance,loglik,status)
   call fit_dcc(returns,dcc,max_iterations=40)
   correlation=dcc%correlations

   print '(a,f8.5)','EWMA lambda: ',lambda
   print '(a,2(1x,f9.5))','last EWMA variances:',covariance(1,1,n),covariance(2,2,n)
   print '(a,2(1x,f9.5))','DCC alpha beta:',dcc%alpha,dcc%beta
   print '(a,f9.5)','last conditional correlation:',correlation(1,2,n)
end program volatility_analysis
