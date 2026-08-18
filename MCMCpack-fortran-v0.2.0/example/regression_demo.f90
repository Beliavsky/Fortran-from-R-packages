program regression_demo
   use mcmcpack, only : dp, set_seed, mcmc_result, mcmc_regress
   implicit none
   real(dp) :: x(8,2),y(8),beta_start(2),b0(2),b0prec(2,2),postmean(3)
   type(mcmc_result) :: fit
   integer :: i

   x(:,1)=1.0_dp
   x(:,2)=[(-2.0_dp+real(i-1,dp),i=1,8)]
   y=1.25_dp+0.6_dp*x(:,2)+[0.10_dp,-0.10_dp,0.05_dp,0.00_dp,0.03_dp,-0.04_dp,0.08_dp,-0.02_dp]
   beta_start=0.0_dp;b0=0.0_dp;b0prec=0.0_dp;b0prec(1,1)=0.01_dp;b0prec(2,2)=0.01_dp

   call set_seed(12345)
   fit=mcmc_regress(y,x,beta_start,b0,b0prec,2.0_dp,1.0_dp,200,2000,2)
   if(fit%status/=0) error stop 'mcmc_regress failed'
   postmean=sum(fit%draws,dim=1)/real(size(fit%draws,1),dp)
   write(*,'(a,3f12.6)') 'posterior means (beta0,beta1,sigma2): ',postmean
end program regression_demo
