program test_panel_hmm_special
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : set_seed
   use mcmcpack_panel_hmm, only : panel_fe_result,panel_re_result,hmm_panel_fe,hmm_panel_re
   implicit none
   type(panel_fe_result)::fe
   type(panel_re_result)::re
   integer,parameter::ns=3,nt=8
   real(dp)::y(ns,nt),x(ns,nt,1),w(ns,nt,1),bst(1),dst(ns,2),sst(ns,2),b0(1),bp(1,1),ps(ns,2,2),p0(ns,2,2)
   integer::nst(ns),i,t
   real(dp)::br(1,2),sr(2),dr(1,1,2),rscale(1,1),pr(2,2),pr0(2,2)
   call set_seed(777)
   do i=1,ns;do t=1,nt
      x(i,t,1)=real(t-4,dp)/4.0_dp;w(i,t,1)=1.0_dp
      if(t<=4)then;y(i,t)=0.5_dp*x(i,t,1)+0.2_dp*real(i-2,dp)
      else;y(i,t)=0.5_dp*x(i,t,1)+1.5_dp+0.2_dp*real(i-2,dp);end if
   end do;end do
   bst=0.0_dp;dst=0.0_dp;dst(:,2)=1.0_dp;sst=1.0_dp;b0=0.0_dp;bp=0.2_dp;nst=2
   ps=0.0_dp;p0=0.0_dp
   do i=1,ns;ps(i,1,1)=0.8_dp;ps(i,1,2)=0.2_dp;ps(i,2,2)=1.0_dp;p0(i,1,1)=8.0_dp;p0(i,1,2)=2.0_dp;p0(i,2,2)=1.0_dp;end do
   fe=hmm_panel_fe(y,x,nst,bst,dst,sst,b0,bp,0.0_dp,0.2_dp,2.0_dp,1.0_dp,ps,p0,5,10,1)
   if(fe%status/=0.or.size(fe%beta,1)/=10.or.any(fe%state<1).or.any(fe%state>2))error stop 1
   br=0.0_dp;br(1,2)=1.0_dp;sr=1.0_dp;dr=0.5_dp;rscale=1.0_dp;pr=0.0_dp;pr0=0.0_dp
   pr(1,1)=0.8_dp;pr(1,2)=0.2_dp;pr(2,2)=1.0_dp;pr0(1,1)=8.0_dp;pr0(1,2)=2.0_dp;pr0(2,2)=1.0_dp
   re=hmm_panel_re(y,x,w,2,br,sr,dr,b0,bp,2.0_dp,1.0_dp,3.0_dp,rscale,pr,pr0,5,10,1)
   if(re%status/=0.or.size(re%beta,1)/=10.or.any(re%state<1).or.any(re%state>2))error stop 2
   if(any(re%sigma2<=0.0_dp).or.any(re%d<=0.0_dp))error stop 3
   print *,'test_panel_hmm_special: PASS'
end program
