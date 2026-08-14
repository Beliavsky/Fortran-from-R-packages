program test_rlr
   use rebayes_kinds, only : dp
   use rebayes_rlr
   implicit none
   real(dp)::x(8,2),y(8),d(2,2)
   type(rlr_result)::r
   x(:,1)=1.0_dp
   x(:,2)=[-2.0_dp,-1.5_dp,-1.0_dp,-0.5_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp]
   y=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp]
   d=0.0_dp;d(1,1)=1.0_dp;d(2,2)=1.0_dp
   call rlr_fit(x,y,d,0.2_dp,r)
   if(r%coef(2)<=0.0_dp)error stop "rlr slope"
   if(.not.(r%loglik<0.0_dp))error stop "rlr loglik"
   print *,"test_rlr: PASS"
end program test_rlr
