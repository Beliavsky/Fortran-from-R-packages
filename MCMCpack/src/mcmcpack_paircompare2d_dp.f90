! SPDX-License-Identifier: GPL-3.0-only
! Truncated Dirichlet-process two-dimensional paired-comparison sampler.
! Translation of cMCMCpaircompare2dDP.cc.
module mcmcpack_paircompare2d_dp
   use mcmcpack_kinds, only : dp,pi
   use mcmcpack_rng, only : runif,rnorm,rtruncnorm,rbeta,rgamma_mt
   use mcmcpack_special_utils, only : sample_categorical
   use mcmcpack_paircompare2d, only : update_theta_candidate
   implicit none
   private
   public :: paircompare2d_dp_result,mcmc_paircompare2d_dp
   type :: paircompare2d_dp_result
      real(dp),allocatable :: draws(:,:)
      real(dp),allocatable :: gamma_accept_rate(:)
      integer :: status=0
   end type paircompare2d_dp_result
contains
   real(dp) function respondent_ll(r,g,md,ystar,theta) result(v)
      integer,intent(in)::r,md(:,:)
      real(dp),intent(in)::g,ystar(:),theta(:,:)
      integer::i,c1,c2
      real(dp)::eta
      v=0.0_dp
      do i=1,size(md,1)
         if(md(i,1)/=r)cycle;c1=md(i,2);c2=md(i,3)
         eta=cos(g)*(theta(c1,1)-theta(c2,1))+sin(g)*(theta(c1,2)-theta(c2,2))
         v=v-0.5_dp*(ystar(i)-eta)**2
      end do
   end function respondent_ll

   real(dp) function cluster_ll(k,g,membership,md,ystar,theta) result(v)
      integer,intent(in)::k,membership(:),md(:,:)
      real(dp),intent(in)::g,ystar(:),theta(:,:)
      integer::r
      v=0.0_dp
      do r=1,size(membership)
         if(membership(r)==k)v=v+respondent_ll(r,g,md,ystar,theta)
      end do
   end function cluster_ll

   function mcmc_paircompare2d_dp(md,theta_start,cluster_gamma_start,membership_start,theta_eq,theta_ineq, &
                                  tune,cluster_mcmc,alpha_start,alpha_fixed,a0,b0,burnin,mcmc,thin, &
                                  store_theta,store_gamma) result(res)
      integer,intent(in)::md(:,:),membership_start(:),cluster_mcmc,burnin,mcmc,thin
      real(dp),intent(in)::theta_start(:,:),cluster_gamma_start(:),theta_eq(:,:),theta_ineq(:,:),tune,alpha_start,a0,b0
      logical,intent(in)::alpha_fixed
      logical,intent(in),optional::store_theta,store_gamma
      type(paircompare2d_dp_result)::res
      integer::n,ncand,nresp,kmax,nstore,iter,keep,i,r,k,h,c1,c2,info,ncol,col,nuniq,newk
      integer::membership(size(membership_start)),csize(size(cluster_gamma_start))
      real(dp)::theta(size(theta_start,1),2),cgamma(size(cluster_gamma_start)),gamma(size(membership_start))
      real(dp)::ystar(size(md,1)),weights(size(cluster_gamma_start)),vstick(size(cluster_gamma_start))
      real(dp)::logp(size(cluster_gamma_start)),prob(size(cluster_gamma_start))
      real(dp)::trial(size(membership_start)),accept(size(membership_start))
      real(dp)::mu,gold,gnew,llold,llnew,alpha,remaining,aa,bb,u
      logical::st,sg
      st=.true.;sg=.true.;if(present(store_theta))st=store_theta;if(present(store_gamma))sg=store_gamma
      n=size(md,1);ncand=size(theta_start,1);nresp=size(membership_start);kmax=size(cluster_gamma_start);nstore=mcmc/thin
      if(size(md,2)<4.or.size(theta_start,2)/=2.or.any(shape(theta_eq)/=[ncand,2]).or.any(shape(theta_ineq)/=[ncand,2]).or. &
         kmax<1.or.any(membership_start<1).or.any(membership_start>kmax).or.any(cluster_gamma_start<0.0_dp).or. &
         any(cluster_gamma_start>0.5_dp*pi).or.tune<=0.0_dp.or.cluster_mcmc<1.or. &
         alpha_start<=0.0_dp.or.a0<=0.0_dp.or.b0<=0.0_dp.or. &
         nstore<=0.or.(.not.st.and..not.sg))then;res%status=1;return;end if
      if(any(md(:,1)<1).or.any(md(:,1)>nresp).or.any(md(:,2)<1).or.any(md(:,2)>ncand).or.any(md(:,3)<1).or.any(md(:,3)>ncand))then
         res%status=2;return
      end if
      ncol=merge(2*ncand,0,st)+merge(2*nresp+1,0,sg)+merge(0,1,alpha_fixed)
      allocate(res%draws(nstore,ncol),res%gamma_accept_rate(nresp));res%draws=0.0_dp
      theta=theta_start;cgamma=cluster_gamma_start;membership=membership_start;alpha=alpha_start;trial=0.0_dp;accept=0.0_dp
      csize=0;do r=1,nresp;csize(membership(r))=csize(membership(r))+1;end do
      do r=1,nresp;gamma(r)=cgamma(membership(r));end do
      keep=0
      do iter=0,burnin+mcmc-1
         ! latent probit utilities
         do i=1,n
            r=md(i,1);c1=md(i,2);c2=md(i,3)
            mu=cos(gamma(r))*(theta(c1,1)-theta(c2,1))+sin(gamma(r))*(theta(c1,2)-theta(c2,2))
            if(md(i,4)==c1)then;ystar(i)=rtruncnorm(mu,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else if(md(i,4)==c2)then;ystar(i)=rtruncnorm(mu,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
            else;ystar(i)=rnorm(mu,1.0_dp);end if
         end do
         ! cluster-specific angle updates
         do k=1,kmax
            if(csize(k)==0)then
               cgamma(k)=0.5_dp*pi*runif()
            else
               gold=cgamma(k)
               do h=1,cluster_mcmc
                  do;gnew=gold+(1.0_dp-2.0_dp*runif())*tune;if(gnew>=0.0_dp.and.gnew<=0.5_dp*pi)exit;end do
                  llold=cluster_ll(k,gold,membership,md,ystar,theta);llnew=cluster_ll(k,gnew,membership,md,ystar,theta)
                  do r=1,nresp;if(membership(r)==k)trial(r)=trial(r)+1.0_dp;end do
                  if(log(runif())<min(0.0_dp,llnew-llold))then
                     gold=gnew;do r=1,nresp;if(membership(r)==k)accept(r)=accept(r)+1.0_dp;end do
                  end if
               end do
               cgamma(k)=gold
            end if
         end do
         ! blocked stick-breaking weights V_k | allocations.
         remaining=real(nresp,dp);weights=0.0_dp;vstick=0.0_dp;u=1.0_dp
         do k=1,kmax-1
            aa=1.0_dp+real(csize(k),dp);remaining=remaining-real(csize(k),dp);bb=alpha+remaining
            vstick(k)=min(0.9999_dp,rbeta(aa,bb));weights(k)=u*vstick(k);u=u*(1.0_dp-vstick(k))
         end do
         weights(kmax)=u;vstick(kmax)=1.0_dp
         ! respondent allocations
         do r=1,nresp
            do k=1,kmax;logp(k)=log(max(weights(k),tiny(1.0_dp)))+respondent_ll(r,cgamma(k),md,ystar,theta);end do
            prob=exp(logp-maxval(logp));newk=sample_categorical(prob)
            if(newk/=membership(r))then
               csize(membership(r))=csize(membership(r))-1
               csize(newk)=csize(newk)+1;membership(r)=newk
            end if
         end do
         do r=1,nresp;gamma(r)=cgamma(membership(r));end do
         do i=1,ncand
            call update_theta_candidate(i,md,ystar,gamma,theta,theta_eq,theta_ineq,info)
            if(info/=0)then;res%status=10+info;return;end if
         end do
         if(.not.alpha_fixed)then
            alpha=rgamma_mt(a0+real(kmax-1,dp),1.0_dp/max(b0-log(max(weights(kmax),tiny(1.0_dp))),tiny(1.0_dp)))
         end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            if(st)then;res%draws(keep,1:ncand)=theta(:,1);res%draws(keep,ncand+1:2*ncand)=theta(:,2);col=2*ncand;end if
            if(sg)then
               res%draws(keep,col+1:col+nresp)=gamma;col=col+nresp
               res%draws(keep,col+1:col+nresp)=real(membership,dp);col=col+nresp
               nuniq=count(csize>0);res%draws(keep,col+1)=real(nuniq,dp);col=col+1
            end if
            if(.not.alpha_fixed)res%draws(keep,col+1)=alpha
         end if
      end do
      do r=1,nresp
         if(trial(r)>0.0_dp)then
            res%gamma_accept_rate(r)=accept(r)/trial(r)
         else
            res%gamma_accept_rate(r)=0.0_dp
         end if
      end do
   end function mcmc_paircompare2d_dp
end module mcmcpack_paircompare2d_dp
