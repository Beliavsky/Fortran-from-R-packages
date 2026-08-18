program test_hdp_special
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_hdp_hmm, only : hdp_count_result,hdphmm_poisson,hdphmm_negbin,hdphsmm_negbin
   implicit none
   integer,parameter::n=18,k=3
   integer::y(n),i
   real(dp)::x(n,1),b(k,1),p(k,k),b0(1),bp(1,1),rho(k),rs(k),om(k),ph(k,k)
   type(hdp_count_result)::po,nb,hs
   call set_seed(9090)
   x=1.0_dp
   y=[0,1,0,1,0,1, 5,7,6,4,8,6, 1,0,2,1,0,1]
   b(:,1)=[-0.2_dp,1.5_dp,0.3_dp];p=0.1_dp
   do i=1,k;p(i,i)=0.8_dp;p(i,:)=p(i,:)/sum(p(i,:));end do
   b0=0.0_dp;bp=0.2_dp
   po=hdphmm_poisson(y,x,k,b,p,2.0_dp,5.0_dp,0.8_dp,b0,bp,1.0_dp,0.1_dp,1.0_dp,0.1_dp,10.0_dp,2.0_dp,5,10,1,0.1_dp)
   if(po%status/=0.or.size(po%state,1)/=10.or.any(po%state<1).or.any(po%state>k))error stop 1
   rho=[2.0_dp,2.0_dp,2.0_dp];rs=0.3_dp
   nb=hdphmm_negbin(y,x,k,b,p,rho,2.0_dp,5.0_dp,0.8_dp,b0,bp,1.0_dp,0.1_dp,1.0_dp,0.1_dp,10.0_dp,2.0_dp, &
                    2.0_dp,2.0_dp,10.0_dp,rs,5,10,1,0.08_dp)
   if(nb%status/=0.or.any(nb%rho<=0.0_dp))error stop 2
   ph=1.0_dp
   do i=1,k;ph(i,i)=0.0_dp;ph(i,:)=ph(i,:)/sum(ph(i,:));end do
   om=[0.4_dp,0.5_dp,0.6_dp]
   hs=hdphsmm_negbin(y,x,k,b,ph,rho,om,2.0_dp,5.0_dp,b0,bp,1.0_dp,0.1_dp,1.0_dp,0.1_dp,1.0_dp,5.0_dp, &
                     10.0_dp,2.0_dp,2.0_dp,1.0_dp,rs,3,6,1,0.08_dp)
   if(hs%status/=0.or.any(hs%rho<=0.0_dp).or.any(hs%omega<=0.0_dp).or.any(hs%omega>=1.0_dp))error stop 3
   print *,'test_hdp_special: PASS'
end program
