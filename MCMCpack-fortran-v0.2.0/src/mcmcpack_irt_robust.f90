! SPDX-License-Identifier: GPL-3.0-only
! Robust multidimensional IRT sampler translated from MCMCirtKdRob.cc.
! It keeps MCMCpack's four-parameter logistic response model and constraints.
! Coordinate random-walk Metropolis updates replace the original coordinate
! slice sampler; the posterior target and parameterization are unchanged.
module mcmcpack_irt_robust
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : runif,rnorm
   use mcmcpack_special_utils, only : logistic_safe,normal_logpdf
   implicit none
   private
   public :: irt_robust_result,mcmc_irtkd_rob
   type :: irt_robust_result
      real(dp), allocatable :: draws(:,:)
      real(dp) :: accept_rate=0.0_dp
      integer :: status=0
   end type irt_robust_result
contains
   real(dp) function robust_loglike(x,lambda,theta,delta0,delta1) result(v)
      integer,intent(in)::x(:,:)
      real(dp),intent(in)::lambda(:,:),theta(:,:),delta0,delta1
      integer::n,k,d,i,j,h
      real(dp)::eta,p
      n=size(x,1);k=size(x,2);d=size(theta,2);v=0.0_dp
      if(size(lambda,1)/=k.or.size(lambda,2)/=d+1.or.delta0<=0.0_dp.or.delta1<=0.0_dp.or.delta0+delta1>=1.0_dp)then
         v=-huge(1.0_dp);return
      end if
      do i=1,n;do j=1,k
         if(x(i,j)==-999)cycle
         eta=-lambda(j,1)
         do h=1,d;eta=eta+theta(i,h)*lambda(j,h+1);end do
         p=delta0+(1.0_dp-delta0-delta1)*logistic_safe(eta)
         p=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),p))
         if(x(i,j)==1)then;v=v+log(p)
         else if(x(i,j)==0)then;v=v+log(1.0_dp-p)
         else;v=-huge(1.0_dp);return;end if
      end do;end do
   end function robust_loglike

   real(dp) function lambda_logpost(x,lambda,theta,delta0,delta1,pm,pp) result(v)
      integer,intent(in)::x(:,:)
      real(dp),intent(in)::lambda(:,:),theta(:,:),delta0,delta1,pm(:,:),pp(:,:)
      integer::i,j
      v=robust_loglike(x,lambda,theta,delta0,delta1)
      if(v<=-0.5_dp*huge(1.0_dp))return
      do i=1,size(lambda,1);do j=1,size(lambda,2)
         if(pp(i,j)>0.0_dp)v=v-0.5_dp*pp(i,j)*(lambda(i,j)-pm(i,j))**2
      end do;end do
   end function lambda_logpost

   real(dp) function full_logpost(x,lambda,theta,delta0,delta1,pm,pp,k0,k1,c0,d0,c1,d1) result(v)
      integer,intent(in)::x(:,:)
      real(dp),intent(in)::lambda(:,:),theta(:,:),delta0,delta1,pm(:,:),pp(:,:),k0,k1,c0,d0,c1,d1
      real(dp)::q0,q1
      if(delta0<=0.0_dp.or.delta0>=k0.or.delta1<=0.0_dp.or.delta1>=k1.or.delta0+delta1>=1.0_dp)then
         v=-huge(1.0_dp);return
      end if
      v=lambda_logpost(x,lambda,theta,delta0,delta1,pm,pp)
      if(v<=-0.5_dp*huge(1.0_dp))return
      v=v-0.5_dp*sum(theta**2)
      q0=delta0/k0;q1=delta1/k1
      v=v+(c0-1.0_dp)*log(q0)+(d0-1.0_dp)*log(1.0_dp-q0) &
         +(c1-1.0_dp)*log(q1)+(d1-1.0_dp)*log(1.0_dp-q1)
   end function full_logpost

   function mcmc_irtkd_rob(x,lambda_start,theta_start,lambda_eq,lambda_ineq,theta_eq,theta_ineq, &
                           lambda_prior_mean,lambda_prior_prec,delta0_start,delta1_start,k0,k1,c0,d0,c1,d1, &
                           theta_tune,lambda_tune,delta0_tune,delta1_tune,burnin,mcmc,thin,store_item,store_ability) result(res)
      integer,intent(in)::x(:,:),burnin,mcmc,thin
      real(dp),intent(in)::lambda_start(:,:),theta_start(:,:),lambda_eq(:,:),lambda_ineq(:,:),theta_eq(:,:),theta_ineq(:,:)
      real(dp),intent(in)::lambda_prior_mean(:,:),lambda_prior_prec(:,:),delta0_start,delta1_start,k0,k1,c0,d0,c1,d1
      real(dp),intent(in)::theta_tune,lambda_tune,delta0_tune,delta1_tune
      logical,intent(in),optional::store_item,store_ability
      type(irt_robust_result)::res
      integer::n,ki,d,nstore,iter,keep,i,j,ncol,col,accepts,trials
      real(dp)::lambda(size(lambda_start,1),size(lambda_start,2)),theta(size(theta_start,1),size(theta_start,2))
      real(dp)::delta0,delta1,cur,canlp,old,can
      logical::sti,sta
      n=size(x,1);ki=size(x,2);d=size(theta_start,2);nstore=mcmc/thin
      sti=.true.;if(present(store_item))sti=store_item
      sta=.false.;if(present(store_ability))sta=store_ability
      if(size(lambda_start,1)/=ki.or.size(lambda_start,2)/=d+1.or.size(theta_start,1)/=n.or. &
         any(shape(lambda_eq)/=shape(lambda_start)).or.any(shape(lambda_ineq)/=shape(lambda_start)).or. &
         any(shape(theta_eq)/=shape(theta_start)).or.any(shape(theta_ineq)/=shape(theta_start)).or. &
         any(shape(lambda_prior_mean)/=shape(lambda_start)).or.any(shape(lambda_prior_prec)/=shape(lambda_start)).or. &
         k0<=0.0_dp.or.k0>0.5_dp.or.k1<=0.0_dp.or.k1>0.5_dp.or.c0<=0.0_dp.or.d0<=0.0_dp.or.c1<=0.0_dp.or.d1<=0.0_dp.or. &
         theta_tune<=0.0_dp.or.lambda_tune<=0.0_dp.or.delta0_tune<=0.0_dp.or. &
         delta1_tune<=0.0_dp.or.nstore<=0.or.(.not.sti.and..not.sta))then
         res%status=1;return
      end if
      lambda=lambda_start;theta=theta_start;delta0=delta0_start;delta1=delta1_start
      do i=1,ki;do j=1,d+1;if(lambda_eq(i,j)>-998.5_dp)lambda(i,j)=lambda_eq(i,j);end do;end do
      do i=1,n;do j=1,d;if(theta_eq(i,j)>-998.5_dp)theta(i,j)=theta_eq(i,j);end do;end do
      ncol=merge(ki*(d+1),0,sti)+merge(n*d,0,sta)+2
      allocate(res%draws(nstore,ncol));res%draws=0.0_dp;keep=0;accepts=0;trials=0
      cur=full_logpost(x,lambda,theta,delta0,delta1,lambda_prior_mean,lambda_prior_prec,k0,k1,c0,d0,c1,d1)
      do iter=0,burnin+mcmc-1
         do i=1,n;do j=1,d
            if(theta_eq(i,j)>-998.5_dp)cycle
            old=theta(i,j);can=old+theta_tune*rnorm();trials=trials+1
            if(theta_ineq(i,j)*can<0.0_dp)cycle
            theta(i,j)=can;canlp=full_logpost(x,lambda,theta,delta0,delta1,lambda_prior_mean,lambda_prior_prec,k0,k1,c0,d0,c1,d1)
            if(log(runif())<min(0.0_dp,canlp-cur))then;cur=canlp;accepts=accepts+1
            else;theta(i,j)=old;end if
         end do;end do
         do i=1,ki;do j=1,d+1
            if(lambda_eq(i,j)>-998.5_dp)cycle
            old=lambda(i,j);can=old+lambda_tune*rnorm();trials=trials+1
            if(lambda_ineq(i,j)*can<0.0_dp)cycle
            lambda(i,j)=can;canlp=full_logpost(x,lambda,theta,delta0,delta1,lambda_prior_mean,lambda_prior_prec,k0,k1,c0,d0,c1,d1)
            if(log(runif())<min(0.0_dp,canlp-cur))then;cur=canlp;accepts=accepts+1
            else;lambda(i,j)=old;end if
         end do;end do
         old=delta0;can=old+delta0_tune*rnorm();trials=trials+1;delta0=can
         canlp=full_logpost(x,lambda,theta,delta0,delta1,lambda_prior_mean,lambda_prior_prec,k0,k1,c0,d0,c1,d1)
         if(log(runif())<min(0.0_dp,canlp-cur))then;cur=canlp;accepts=accepts+1
         else;delta0=old;end if
         old=delta1;can=old+delta1_tune*rnorm();trials=trials+1;delta1=can
         canlp=full_logpost(x,lambda,theta,delta0,delta1,lambda_prior_mean,lambda_prior_prec,k0,k1,c0,d0,c1,d1)
         if(log(runif())<min(0.0_dp,canlp-cur))then;cur=canlp;accepts=accepts+1
         else;delta1=old;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            if(sti)then
               do i=1,ki;res%draws(keep,col+1:col+d+1)=lambda(i,:);col=col+d+1;end do
            end if
            if(sta)then
               do i=1,n;res%draws(keep,col+1:col+d)=theta(i,:);col=col+d;end do
            end if
            res%draws(keep,col+1:col+2)=[delta0,delta1]
         end if
      end do
      if(trials>0)res%accept_rate=real(accepts,dp)/real(trials,dp)
   end function mcmc_irtkd_rob
end module mcmcpack_irt_robust
