! SPDX-License-Identifier: GPL-3.0-only
! Dynamic one-dimensional IRT sampler translated from cMCMCdynamicIRT1d.cc.
module mcmcpack_dynamic_irt
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : rnorm,rtruncnorm,rmvnorm,rinvgamma_rng
   use mcmcpack_linalg, only : inv_spd
   implicit none
   private
   public :: dynamic_irt_result,mcmc_dynamic_irt1d
   type :: dynamic_irt_result
      real(dp), allocatable :: draws(:,:)
      integer :: status=0
   end type dynamic_irt_result
contains
   real(dp) function constrained_normal(mu,var,eq,ineq) result(x)
      real(dp),intent(in)::mu,var,eq,ineq
      real(dp),parameter::big=huge(1.0_dp)/10.0_dp
      if(eq>-998.5_dp)then;x=eq;return;end if
      if(ineq>0.0_dp)then;x=rtruncnorm(mu,sqrt(max(var,tiny(1.0_dp))),0.0_dp,big)
      else if(ineq<0.0_dp)then;x=rtruncnorm(mu,sqrt(max(var,tiny(1.0_dp))),-big,0.0_dp)
      else;x=rnorm(mu,sqrt(max(var,tiny(1.0_dp))));end if
   end function constrained_normal

   function mcmc_dynamic_irt1d(x,item_time,theta_start,alpha_start,beta_start,tau2_start, &
                               a0_mean,a0_prec,b0_mean,b0_prec,c0,d0,init_mean,init_var,theta_eq,theta_ineq, &
                               burnin,mcmc,thin,store_ability,store_item) result(res)
      integer,intent(in)::x(:,:),item_time(:),burnin,mcmc,thin
      real(dp),intent(in)::theta_start(:,:),alpha_start(:),beta_start(:),tau2_start(:)
      real(dp),intent(in)::a0_mean,a0_prec,b0_mean,b0_prec,c0(:),d0(:),init_mean(:),init_var(:),theta_eq(:),theta_ineq(:)
      logical,intent(in),optional::store_ability,store_item
      type(dynamic_irt_result)::res
      integer::ns,ni,nt,nstore,iter,keep,s,i,t,nobs,info,ncol,col
      real(dp)::theta(size(theta_start,1),size(theta_start,2)),alpha(size(alpha_start)),beta(size(beta_start))
      real(dp)::tau2(size(tau2_start)),z(size(x,1),size(x,2)),mu
      real(dp)::prec2(2,2),cov2(2,2),rhs2(2),mean2(2),draw2(2)
      real(dp),allocatable::m(:),C(:),ap(:),R(:)
      real(dp)::obs_prec,obs_rhs,sse,Bgain,back_mean,back_var
      logical::sa,si
      ns=size(x,1);ni=size(x,2);nt=size(theta_start,2);nstore=mcmc/thin
      sa=.true.;if(present(store_ability))sa=store_ability
      si=.true.;if(present(store_item))si=store_item
      if(size(item_time)/=ni.or.any(item_time<1).or.any(item_time>nt).or.size(theta_start,1)/=ns.or. &
         size(alpha_start)/=ni.or.size(beta_start)/=ni.or.size(tau2_start)/=ns.or.size(c0)/=ns.or.size(d0)/=ns.or. &
         size(init_mean)/=ns.or.size(init_var)/=ns.or.size(theta_eq)/=ns.or. &
         size(theta_ineq)/=ns.or.a0_prec<0.0_dp.or.b0_prec<0.0_dp.or. &
         any(init_var<=0.0_dp).or.any(tau2_start<=0.0_dp).or.thin<=0.or.nstore<=0)then;res%status=1;return;end if
      theta=theta_start;alpha=alpha_start;beta=beta_start;tau2=tau2_start;z=-999.0_dp
      ncol=merge(ns*nt,0,sa)+merge(2*ni,0,si)+ns
      allocate(res%draws(nstore,ncol));res%draws=0.0_dp;keep=0
      allocate(m(nt),C(nt),ap(nt),R(nt))
      do iter=0,burnin+mcmc-1
         ! Albert-Chib latent utilities.
         do s=1,ns;do i=1,ni
            if(x(s,i)==-999)cycle
            t=item_time(i);mu=-alpha(i)+beta(i)*theta(s,t)
            if(x(s,i)==1)then;z(s,i)=rtruncnorm(mu,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else if(x(s,i)==0)then;z(s,i)=rtruncnorm(mu,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
            else;res%status=2;return;end if
         end do;end do

         ! Item difficulty/discrimination pairs.
         do i=1,ni
            prec2=0.0_dp;rhs2=0.0_dp;prec2(1,1)=a0_prec;prec2(2,2)=b0_prec
            rhs2(1)=a0_prec*a0_mean;rhs2(2)=b0_prec*b0_mean;nobs=0;t=item_time(i)
            do s=1,ns
               if(x(s,i)==-999)cycle;nobs=nobs+1
               prec2(1,1)=prec2(1,1)+1.0_dp
               prec2(1,2)=prec2(1,2)-theta(s,t);prec2(2,1)=prec2(1,2)
               prec2(2,2)=prec2(2,2)+theta(s,t)**2
               rhs2(1)=rhs2(1)-z(s,i);rhs2(2)=rhs2(2)+theta(s,t)*z(s,i)
            end do
            if(nobs>0)then
               call inv_spd(prec2,cov2,info);if(info/=0)then;res%status=10+info;return;end if
               mean2=matmul(cov2,rhs2);call rmvnorm(mean2,cov2,draw2,info);if(info/=0)then;res%status=20+info;return;end if
               alpha(i)=draw2(1);beta(i)=draw2(2)
            end if
         end do

         ! Scalar DLM forward-filter/backward-sample for every subject.
         do s=1,ns
            if(theta_eq(s)>-998.5_dp)then
               theta(s,:)=theta_eq(s)
            else
               do t=1,nt
                  if(t==1)then;ap(t)=init_mean(s);R(t)=init_var(s)
                  else;ap(t)=m(t-1);R(t)=C(t-1)+tau2(s);end if
                  obs_prec=0.0_dp;obs_rhs=0.0_dp
                  do i=1,ni
                     if(item_time(i)/=t.or.x(s,i)==-999)cycle
                     obs_prec=obs_prec+beta(i)**2
                     obs_rhs=obs_rhs+beta(i)*(z(s,i)+alpha(i))
                  end do
                  if(obs_prec>0.0_dp)then
                     C(t)=1.0_dp/(1.0_dp/R(t)+obs_prec)
                     m(t)=C(t)*(ap(t)/R(t)+obs_rhs)
                  else
                     C(t)=R(t);m(t)=ap(t)
                  end if
               end do
               theta(s,nt)=constrained_normal(m(nt),C(nt),theta_eq(s),theta_ineq(s))
               do t=nt-1,1,-1
                  Bgain=C(t)/max(C(t)+tau2(s),tiny(1.0_dp))
                  back_mean=m(t)+Bgain*(theta(s,t+1)-m(t))
                  back_var=C(t)-Bgain*C(t);back_var=max(back_var,1.0e-12_dp)
                  theta(s,t)=constrained_normal(back_mean,back_var,theta_eq(s),theta_ineq(s))
               end do
            end if
            if(c0(s)>0.0_dp.and.d0(s)>0.0_dp.and.nt>1)then
               sse=sum((theta(s,2:nt)-theta(s,1:nt-1))**2)
               tau2(s)=rinvgamma_rng(0.5_dp*(c0(s)+real(nt-1,dp)),0.5_dp*(d0(s)+sse))
            end if
         end do

         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            if(sa)then
               do s=1,ns;res%draws(keep,col+1:col+nt)=theta(s,:);col=col+nt;end do
            end if
            if(si)then
               res%draws(keep,col+1:col+ni)=alpha;col=col+ni
               res%draws(keep,col+1:col+ni)=beta;col=col+ni
            end if
            res%draws(keep,col+1:col+ns)=tau2
         end if
      end do
   end function mcmc_dynamic_irt1d
end module mcmcpack_dynamic_irt
