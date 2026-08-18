program test_samplers
   use mcmcpack
   implicit none
   real(dp) :: x(8,2),y(8),b0(2),p0(2,2),bstart(2),v(2,2),tune(2,2),yb(8)
   real(dp) :: gamma(4)
   integer :: i
   integer :: yo(8)
   type(mcmc_result) :: fit
   call set_seed(24680)
   do i=1,8
      x(i,1)=1.0_dp;x(i,2)=real(i-4,dp)/3.0_dp;y(i)=1.0_dp+1.5_dp*x(i,2)+0.2_dp*sin(real(i,dp))
   end do
   b0=0.0_dp;p0=0.0_dp;p0(1,1)=0.01_dp;p0(2,2)=0.01_dp;bstart=0.0_dp
   fit=mcmc_regress(y,x,bstart,b0,p0,0.01_dp,0.01_dp,100,400,2)
   if(fit%status/=0.or.any(shape(fit%draws)/=[200,3])) error stop 1
   yb=merge(1.0_dp,0.0_dp,y>1.0_dp)
   fit=mcmc_probit(yb,x,bstart,b0,p0,100,200,2)
   if(fit%status/=0.or.size(fit%draws,2)/=2) error stop 2
   v=0.0_dp;v(1,1)=1.0_dp;v(2,2)=1.0_dp;tune=0.0_dp;tune(1,1)=0.4_dp;tune(2,2)=0.4_dp
   fit=mcmc_logit(yb,x,bstart,b0,p0,v,tune,50,100,2)
   if(fit%status/=0.or.fit%accept_rate<0.0_dp.or.fit%accept_rate>1.0_dp) error stop 3
   fit=mcmc_tobit(max(y,0.0_dp),x,bstart,b0,p0,0.01_dp,0.01_dp,0.0_dp,huge(1.0_dp)/10.0_dp,50,100,2)
   if(fit%status/=0) error stop 4
   fit=mcmc_quantreg(0.5_dp,y,x,bstart,b0,p0,50,100,2)
   if(fit%status/=0) error stop 5
   yo=[1,1,1,2,2,2,3,3];gamma=[-300.0_dp,0.0_dp,1.0_dp,300.0_dp]
   fit=mcmc_oprobit(yo,x,bstart,gamma,b0,p0,0.1_dp,50,100,2)
   if(fit%status/=0.or.size(fit%draws,2)/=3) error stop 6
   print *, 'test_samplers: PASS'
end program test_samplers
