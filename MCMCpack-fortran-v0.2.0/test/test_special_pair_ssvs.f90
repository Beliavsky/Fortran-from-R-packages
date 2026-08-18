program test_special_pair_ssvs
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_paircompare2d_dp, only : paircompare2d_dp_result,mcmc_paircompare2d_dp
   use mcmcpack_ssvs_quantreg, only : ssvs_quantreg_result,ssvs_quantreg
   implicit none
   type(paircompare2d_dp_result)::pc
   type(ssvs_quantreg_result)::qr
   integer::md(8,4),membership(2),i
   real(dp)::theta(3,2),cg(3),eq(3,2),ineq(3,2)
   real(dp)::x(24,3),y(24),z
   call set_seed(24680)
   md=reshape([1,1,2,1, 1,1,3,3, 1,2,3,2, 1,1,2,2, &
               2,1,2,1, 2,1,3,1, 2,2,3,3, 2,1,2,1],[8,4],order=[2,1])
   theta=reshape([-1.0_dp,0.0_dp, 1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[3,2],order=[2,1])
   cg=[0.2_dp,0.9_dp,1.2_dp];membership=[1,2];eq=-999.0_dp;ineq=0.0_dp
   pc=mcmc_paircompare2d_dp(md,theta,cg,membership,eq,ineq,0.08_dp,2,1.0_dp,.false.,2.0_dp,2.0_dp,5,10,1,.true.,.true.)
   if(pc%status/=0.or.size(pc%draws,1)/=10.or.any(pc%gamma_accept_rate<0.0_dp).or.any(pc%gamma_accept_rate>1.0_dp))error stop 1
   do i=1,24
      z=(real(i,dp)-12.5_dp)/6.0_dp;x(i,:)=[1.0_dp,z,sin(real(i,dp))];y(i)=1.0_dp+1.8_dp*z+0.2_dp*cos(real(i,dp))
   end do
   qr=ssvs_quantreg(y,x,0.5_dp,1,5,12,1,1.0_dp,1.0_dp)
   if(qr%status/=0.or.size(qr%gamma,1)/=12.or.any(qr%gamma(:,1)/=1))error stop 2
   if(any(.not.(qr%beta==qr%beta)).or.any(qr%pi0<=0.0_dp).or.any(qr%pi0>=1.0_dp))error stop 3
   print *,'test_special_pair_ssvs: PASS'
end program
