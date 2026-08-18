program test_special_ei_irt
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_ei, only : ei_result,mcmc_hier_ei,mcmc_dynamic_ei
   use mcmcpack_dynamic_irt, only : dynamic_irt_result,mcmc_dynamic_irt1d
   use mcmcpack_irt_robust, only : irt_robust_result,mcmc_irtkd_rob
   implicit none
   type(ei_result)::eh,ed
   type(dynamic_irt_result)::di
   type(irt_robust_result)::ir
   real(dp)::r0(3),r1(3),c0e(3),c1e(3),w(3,3)
   integer::xx(3,4),it(4)
   real(dp)::ths(3,2),as(4),bs(4),ts(3),c0(3),d0(3),im(3),iv(3),eq(3),ineq(3)
   integer::xr(4,3)
   real(dp)::lam(3,3),thr(4,2),leq(3,3),lineq(3,3),teq(4,2),tineq(4,2),lpm(3,3),lpp(3,3)
   call set_seed(12345)
   r0=[60.0_dp,50.0_dp,70.0_dp];r1=[40.0_dp,50.0_dp,30.0_dp]
   c0e=[55.0_dp,52.0_dp,65.0_dp];c1e=100.0_dp-c0e
   eh=mcmc_hier_ei(r0,r1,c0e,c1e,0.0_dp,2.287656_dp,0.0_dp,2.287656_dp,0.825_dp,0.0105_dp,0.825_dp,0.0105_dp,5,10,1)
   if(eh%status/=0.or.size(eh%draws,1)/=10)error stop 1
   w=0.0_dp;w(1,2)=1;w(2,1)=1;w(2,3)=1;w(3,2)=1
   ed=mcmc_dynamic_ei(r0,r1,c0e,c1e,w,2.0_dp,0.1_dp,2.0_dp,0.1_dp,5,10,1)
   if(ed%status/=0.or.any(.not.(ed%draws==ed%draws)))error stop 2
   xx=reshape([1,0,1, 0,1,1, 1,1,0, 0,0,1],[3,4]);it=[1,1,2,2]
   ths=0.0_dp;as=0.0_dp;bs=1.0_dp;ts=0.5_dp;c0=2.0_dp;d0=1.0_dp;im=0.0_dp;iv=1.0_dp;eq=-999.0_dp;ineq=0.0_dp
   di=mcmc_dynamic_irt1d(xx,it,ths,as,bs,ts,0.0_dp,1.0_dp,0.0_dp,1.0_dp,c0,d0,im,iv,eq,ineq,5,10,1)
   if(di%status/=0.or.size(di%draws,1)/=10)error stop 3
   xr=reshape([1,0,1,0, 0,1,0,1, 1,1,0,0],[4,3]);lam=0.0_dp;lam(:,2)=1.0_dp;lam(:,3)=0.3_dp;thr=0.0_dp
   leq=-999.0_dp;lineq=0.0_dp;teq=-999.0_dp;tineq=0.0_dp;lpm=0.0_dp;lpp=1.0_dp
   ir=mcmc_irtkd_rob(xr,lam,thr,leq,lineq,teq,tineq,lpm,lpp,0.05_dp,0.05_dp,0.2_dp,0.2_dp,2.0_dp,2.0_dp,2.0_dp,2.0_dp, &
                      0.15_dp,0.1_dp,0.02_dp,0.02_dp,5,10,1,.true.,.true.)
   if(ir%status/=0.or.size(ir%draws,1)/=10.or.ir%accept_rate<0.0_dp.or.ir%accept_rate>1.0_dp)error stop 4
   print *,'test_special_ei_irt: PASS'
end program
