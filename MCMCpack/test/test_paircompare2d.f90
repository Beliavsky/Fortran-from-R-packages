program test_paircompare2d
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_paircompare2d, only : paircompare2d_result,mcmc_paircompare2d
   implicit none
   integer :: md(8,4)
   real(dp) :: theta0(3,2),gamma0(2),eq(3,2),ineq(3,2)
   type(paircompare2d_result) :: fit
   md=reshape([ &
      1,1,2,1, 1,1,3,1, 1,2,3,2, 1,1,2,1, &
      2,1,2,2, 2,1,3,1, 2,2,3,3, 2,1,3,1],[8,4],order=[2,1])
   theta0=0.0_dp;theta0(1,1)=0.5_dp;theta0(2,1)=-0.2_dp;theta0(3,2)=0.3_dp
   gamma0=[0.4_dp,1.0_dp]
   eq=-999.0_dp;ineq=0.0_dp
   eq(1,2)=0.0_dp
   ineq(1,1)=1.0_dp
   call set_seed(4811)
   fit=mcmc_paircompare2d(md,theta0,gamma0,eq,ineq,0.15_dp,20,40,2,.true.,.true.)
   if(fit%status/=0)error stop 1
   if(any(shape(fit%draws)/=[20,8]))error stop 2
   if(any(fit%draws(:,1)<0.0_dp))error stop 3
   if(maxval(abs(fit%draws(:,4)))>1.0e-12_dp)error stop 4
   if(any(fit%draws(:,7:8)<0.0_dp).or.any(fit%draws(:,7:8)>0.5_dp*acos(-1.0_dp)))error stop 5
   if(any(fit%gamma_accept_rate<0.0_dp).or.any(fit%gamma_accept_rate>1.0_dp))error stop 6
   print '(a)', 'test_paircompare2d: PASS'
end program test_paircompare2d
