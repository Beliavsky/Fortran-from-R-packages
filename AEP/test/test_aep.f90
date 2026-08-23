program test_aep
   use aep
   implicit none
   integer,parameter::n=12000
   real(dp)::x(n),u(7),q(7),p(7),m,v,mu,sigma,alpha,epsi
   type(aep_fit_result)::fit
   type(aep_reg_result)::rf
   real(dp)::yr(400),xr(400,1),er(400)
   integer::i
   alpha=2.0_dp;sigma=1.3_dp;mu=-0.4_dp;epsi=0.0_dp
   if(abs(daep(mu,2.0_dp,sigma,mu,0.0_dp)-1.0_dp/(sigma*sqrt(acos(-1.0_dp))))>2e-12_dp) error stop 1
   if(abs(paep(mu,alpha,sigma,mu,epsi)-0.5_dp)>1e-13_dp) error stop 2
   u=[0.01_dp,0.1_dp,0.25_dp,0.5_dp,0.75_dp,0.9_dp,0.99_dp]
   q=qaep(u,1.25_dp,0.8_dp,0.3_dp,0.35_dp);p=paep(q,1.25_dp,0.8_dp,0.3_dp,0.35_dp)
   if(maxval(abs(p-u))>2e-9_dp) error stop 3
   call raep(x,1.5_dp,1.1_dp,0.2_dp,0.25_dp)
   m=sum(x)/n;v=sum((x-m)**2)/(n-1)
   if(.not.(m>-0.1_dp.and.m<0.8_dp.and.v>0.2_dp.and.v<5.0_dp)) error stop 4
   call raep(x,1.25_dp,0.9_dp,-0.2_dp,0.2_dp)
   call fitaep(x(1:2500),fit,max_iter=1000)
   if(abs(fit%mu+0.2_dp)>0.18_dp.or.abs(fit%epsilon-0.2_dp)>0.16_dp) error stop 5
   call raep(er,1.4_dp,0.5_dp,0.0_dp,0.15_dp)
   do i=1,400
      xr(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/399.0_dp
      yr(i)=1.2_dp+0.8_dp*xr(i,1)+er(i)
   end do
   call regaep(yr,xr,rf,max_iter=1000)
   if(abs(rf%beta(1)-1.2_dp)>0.25_dp.or.abs(rf%beta(2)-0.8_dp)>0.18_dp) error stop 6
   if(rf%r2<0.5_dp) error stop 7
   print '(a)', 'test_aep: PASS'
end program
