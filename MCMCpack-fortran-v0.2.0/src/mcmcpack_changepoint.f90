! SPDX-License-Identifier: GPL-3.0-only
! Ordered changepoint samplers translated from MCMCbinaryChange,
! MCMCprobitChange, and MCMCregressChange.
module mcmcpack_changepoint
   use mcmcpack_kinds, only : dp, pi
   use mcmcpack_math, only : normal_cdf
   use mcmcpack_rng, only : runif, rbeta, rtruncnorm, rmvnorm, rinvgamma_rng
   use mcmcpack_negbin, only : negbin_logpost, rho_nb_logcond
   use mcmcpack_linalg, only : inv_spd
   implicit none
   private
   public :: change_result, mcmc_binary_change, mcmc_probit_change, mcmc_regress_change, mcmc_poisson_change
   public :: mcmc_negbin_change, mcmc_residual_break_analysis
   public :: ordered_state_sample

   type :: change_result
      real(dp), allocatable :: draws(:,:)
      integer, allocatable :: states(:,:)
      real(dp), allocatable :: prob_state(:,:)
      integer :: status = 0
   end type change_result
contains
   subroutine sample_discrete_prob(prob,draw)
      real(dp),intent(in)::prob(:)
      integer,intent(out)::draw
      real(dp)::u,c
      integer::j
      u=runif();c=0.0_dp;draw=size(prob)
      do j=1,size(prob);c=c+prob(j);if(u<c)then;draw=j;return;end if;end do
   end subroutine

   subroutine ordered_state_sample(log_emit,p,s,ps,status)
      real(dp),intent(in)::log_emit(:,:),p(:,:)
      integer,intent(out)::s(size(log_emit,1))
      real(dp),intent(out)::ps(size(log_emit,1),size(log_emit,2))
      integer,intent(out)::status
      integer::n,ns,t,st
      real(dp)::f(size(log_emit,1),size(log_emit,2)),pred(size(log_emit,2)),w(size(log_emit,2)),mx,z
      n=size(log_emit,1);ns=size(log_emit,2);status=0
      if(n<1.or.ns<1.or.any(shape(p)/=[ns,ns]))then;status=1;return;end if
      pred=0.0_dp;pred(1)=1.0_dp
      do t=1,n
         if(t>1)pred=matmul(f(t-1,:),p)
         mx=maxval(log_emit(t,:));w=pred*exp(log_emit(t,:)-mx);z=sum(w)
         if(z<=tiny(1.0_dp))then;status=2;return;end if
         f(t,:)=w/z
      end do
      ps(n,:)=f(n,:);s(n)=ns
      do t=n-1,1,-1
         st=s(t+1);w=f(t,:)*p(:,st);z=sum(w)
         if(z<=tiny(1.0_dp))then;status=3;return;end if
         ps(t,:)=w/z;call sample_discrete_prob(ps(t,:),s(t))
      end do
   end subroutine ordered_state_sample

   subroutine update_ordered_p(p,a0,nstate,status)
      real(dp),intent(inout)::p(:,:)
      real(dp),intent(in)::a0(:,:)
      integer,intent(in)::nstate(:)
      integer,intent(out)::status
      integer::ns,j
      real(dp)::aa,bb
      ns=size(p,1);status=0
      if(any(shape(p)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]).or.size(nstate)/=ns)then;status=1;return;end if
      p=0.0_dp;p(ns,ns)=1.0_dp
      do j=1,ns-1
         aa=a0(j,j)+real(nstate(j)-1,dp);bb=a0(j,j+1)+1.0_dp
         if(aa<=0.0_dp.or.bb<=0.0_dp)then;status=2;return;end if
         p(j,j)=rbeta(aa,bb);p(j,j+1)=1.0_dp-p(j,j)
      end do
   end subroutine update_ordered_p

   subroutine flatten_p(row,offset,p)
      real(dp),intent(inout)::row(:)
      integer,intent(in)::offset
      real(dp),intent(in)::p(:,:)
      integer::i,j,k
      k=offset
      do i=1,size(p,1);do j=1,size(p,2);k=k+1;row(k)=p(i,j);end do;end do
   end subroutine

   function mcmc_binary_change(y,phi_start,p_start,a0,c0,d0,burnin,mcmc,thin) result(res)
      integer,intent(in)::y(:),burnin,mcmc,thin
      real(dp),intent(in)::phi_start(:),p_start(:,:),a0(:,:),c0,d0
      type(change_result)::res
      integer::n,ns,nstore,iter,keep,j,info
      integer::s(size(y)),nstate(size(phi_start)),ysum(size(phi_start))
      real(dp)::phi(size(phi_start)),p(size(phi_start),size(phi_start)),loge(size(y),size(phi_start)),ps(size(y),size(phi_start))
      n=size(y);ns=size(phi_start);nstore=mcmc/thin
      if(ns<1.or.nstore<=0.or.any(y<0).or.any(y>1).or.any(shape(p_start)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]) &
         .or.c0<=0.0_dp.or.d0<=0.0_dp)then;res%status=1;return;end if
      allocate(res%draws(nstore,ns+ns*ns),res%states(nstore,n),res%prob_state(n,ns));res%prob_state=0.0_dp
      phi=phi_start;p=p_start;keep=0
      do iter=0,burnin+mcmc-1
         do j=1,ns
            loge(:,j)=real(y,dp)*log(max(phi(j),tiny(1.0_dp)))+real(1-y,dp)*log(max(1.0_dp-phi(j),tiny(1.0_dp)))
         end do
         call ordered_state_sample(loge,p,s,ps,info);if(info/=0)then;res%status=10+info;return;end if
         do j=1,ns
            nstate(j)=count(s==j);ysum(j)=sum(pack(y,s==j));phi(j)=rbeta(c0+real(ysum(j),dp),d0+real(nstate(j)-ysum(j),dp))
         end do
         call update_ordered_p(p,a0,nstate,info);if(info/=0)then;res%status=20+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;res%draws(keep,1:ns)=phi;call flatten_p(res%draws(keep,:),ns,p)
            res%states(keep,:)=s;res%prob_state=res%prob_state+ps/real(nstore,dp)
         end if
      end do
   end function mcmc_binary_change

   function mcmc_probit_change(y,x,beta_start,p_start,b0,b0prec,a0,burnin,mcmc,thin) result(res)
      integer,intent(in)::y(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),p_start(:,:),b0(:),b0prec(:,:),a0(:,:)
      type(change_result)::res
      integer::n,ns,k,nstore,iter,keep,j,i,h,info,nj
      integer::s(size(y)),nstate(size(beta_start,1)),idx(size(y))
      real(dp)::beta(size(beta_start,1),size(beta_start,2)),p(size(beta_start,1),size(beta_start,1))
      real(dp)::loge(size(y),size(beta_start,1)),ps(size(y),size(beta_start,1)),z(size(y)),eta,prob
      real(dp),allocatable::xj(:,:),zj(:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:)
      n=size(y);ns=size(beta_start,1);k=size(beta_start,2);nstore=mcmc/thin
      if(size(x,1)/=n.or.size(x,2)/=k.or.size(b0)/=k.or.any(shape(b0prec)/=[k,k]).or. &
         any(shape(p_start)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]).or. &
         any(y<0).or.any(y>1).or.nstore<=0)then
         res%status=1;return
      end if
      allocate(res%draws(nstore,ns*k+ns*ns),res%states(nstore,n),res%prob_state(n,ns));res%prob_state=0.0_dp
      beta=beta_start;p=p_start;keep=0
      do iter=0,burnin+mcmc-1
         do j=1,ns;do i=1,n
            eta=dot_product(x(i,:),beta(j,:));prob=min(1.0_dp-epsilon(1.0_dp),max(epsilon(1.0_dp),normal_cdf(eta)))
            loge(i,j)=merge(log(prob),log(1.0_dp-prob),y(i)==1)
         end do;end do
         call ordered_state_sample(loge,p,s,ps,info);if(info/=0)then;res%status=10+info;return;end if
         do i=1,n
            eta=max(-200.0_dp,min(200.0_dp,dot_product(x(i,:),beta(s(i),:))))
            if(y(i)==1)then;z(i)=rtruncnorm(eta,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else;z(i)=rtruncnorm(eta,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp);end if
         end do
         do j=1,ns
            nj=count(s==j);nstate(j)=nj;h=0
            do i=1,n;if(s(i)==j)then;h=h+1;idx(h)=i;end if;end do
            allocate(xj(nj,k),zj(nj),prec(k,k),cov(k,k),rhs(k),mu(k),draw(k))
            do h=1,nj;xj(h,:)=x(idx(h),:);zj(h)=z(idx(h));end do
            prec=b0prec+matmul(transpose(xj),xj);rhs=matmul(b0prec,b0)+matmul(transpose(xj),zj)
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=20+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=30+info;return;end if;beta(j,:)=draw
            deallocate(xj,zj,prec,cov,rhs,mu,draw)
         end do
         call update_ordered_p(p,a0,nstate,info);if(info/=0)then;res%status=40+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;h=0;do j=1,ns;res%draws(keep,h+1:h+k)=beta(j,:);h=h+k;end do;call flatten_p(res%draws(keep,:),h,p)
            res%states(keep,:)=s;res%prob_state=res%prob_state+ps/real(nstore,dp)
         end if
      end do
   end function mcmc_probit_change

   function mcmc_regress_change(y,x,beta_start,sigma2_start,p_start,state_start,b0,b0prec,c0,d0,a0, &
                                burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:,:),sigma2_start(:),p_start(:,:),b0(:),b0prec(:,:),c0,d0,a0(:,:)
      integer,intent(in)::state_start(:),burnin,mcmc,thin
      type(change_result)::res
      integer::n,ns,k,nstore,iter,keep,j,i,h,info,nj
      integer::s(size(y)),nstate(size(beta_start,1)),idx(size(y))
      real(dp)::beta(size(beta_start,1),size(beta_start,2)),sig2(size(beta_start,1)),p(size(beta_start,1),size(beta_start,1))
      real(dp)::loge(size(y),size(beta_start,1)),ps(size(y),size(beta_start,1)),eta,e,sse
      real(dp),allocatable::xj(:,:),yj(:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:)
      n=size(y);ns=size(beta_start,1);k=size(beta_start,2);nstore=mcmc/thin
      if(size(x,1)/=n.or.size(x,2)/=k.or.size(sigma2_start)/=ns.or.size(state_start)/=n.or.size(b0)/=k.or. &
         any(shape(b0prec)/=[k,k]).or.any(shape(p_start)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]).or. &
         any(sigma2_start<=0.0_dp).or.nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,ns*k+ns+ns*ns),res%states(nstore,n),res%prob_state(n,ns));res%prob_state=0.0_dp
      beta=beta_start;sig2=sigma2_start;p=p_start;s=state_start;keep=0
      do iter=0,burnin+mcmc-1
         do j=1,ns
            nj=count(s==j);nstate(j)=nj;h=0
            do i=1,n;if(s(i)==j)then;h=h+1;idx(h)=i;end if;end do
            if(nj<=0)then;res%status=2;return;end if
            allocate(xj(nj,k),yj(nj),prec(k,k),cov(k,k),rhs(k),mu(k),draw(k))
            do h=1,nj;xj(h,:)=x(idx(h),:);yj(h)=y(idx(h));end do
            prec=b0prec+matmul(transpose(xj),xj)/sig2(j);rhs=matmul(b0prec,b0)+matmul(transpose(xj),yj)/sig2(j)
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=10+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=20+info;return;end if;beta(j,:)=draw
            sse=0.0_dp;do h=1,nj;e=yj(h)-dot_product(xj(h,:),beta(j,:));sse=sse+e*e;end do
            sig2(j)=rinvgamma_rng(0.5_dp*(c0+real(nj,dp)),0.5_dp*(d0+sse))
            deallocate(xj,yj,prec,cov,rhs,mu,draw)
         end do
         call update_ordered_p(p,a0,nstate,info);if(info/=0)then;res%status=30+info;return;end if
         do i=1,n;do j=1,ns
            eta=dot_product(x(i,:),beta(j,:));loge(i,j)=-0.5_dp*(log(2.0_dp*pi*sig2(j))+(y(i)-eta)**2/sig2(j))
         end do;end do
         call ordered_state_sample(loge,p,s,ps,info);if(info/=0)then;res%status=40+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;h=0;do j=1,ns;res%draws(keep,h+1:h+k)=beta(j,:);h=h+k;end do
            res%draws(keep,h+1:h+ns)=sig2;h=h+ns;call flatten_p(res%draws(keep,:),h,p)
            res%states(keep,:)=s;res%prob_state=res%prob_state+ps/real(nstore,dp)
         end if
      end do
   end function mcmc_regress_change

   real(dp) function poisson_state_logpost(y,x,beta,b0,b0prec) result(v)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:),beta(:),b0(:),b0prec(:,:)
      real(dp)::eta(size(y)),d(size(beta))
      eta=matmul(x,beta)
      if(maxval(eta)>700.0_dp)then;v=-huge(1.0_dp);return;end if
      d=beta-b0;v=sum(real(y,dp)*eta-exp(eta))-0.5_dp*dot_product(d,matmul(b0prec,d))
   end function poisson_state_logpost

   function mcmc_poisson_change(y,x,beta_start,p_start,b0,b0prec,a0,burnin,mcmc,thin,beta_tune) result(res)
      ! Same posterior target as MCMCpoissonChange.  The source's auxiliary
      ! normal-mixture beta step is replaced by a Gaussian MH update.
      integer,intent(in)::y(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),p_start(:,:),b0(:),b0prec(:,:),a0(:,:)
      real(dp),intent(in),optional::beta_tune
      type(change_result)::res
      integer::n,ns,k,nstore,iter,keep,j,i,h,info,nj,accepts,tries
      integer::s(size(y)),nstate(size(beta_start,1)),idx(size(y))
      real(dp)::beta(size(beta_start,1),size(beta_start,2)),p(size(beta_start,1),size(beta_start,1))
      real(dp)::loge(size(y),size(beta_start,1)),ps(size(y),size(beta_start,1)),eta,tune,cur,canlp
      real(dp),allocatable::xj(:,:),yjreal(:),prec(:,:),cov(:,:),zero(:),z(:),can(:)
      integer,allocatable::yj(:)
      n=size(y);ns=size(beta_start,1);k=size(beta_start,2);nstore=mcmc/thin;tune=0.20_dp;if(present(beta_tune))tune=beta_tune
      if(size(x,1)/=n.or.size(x,2)/=k.or.size(b0)/=k.or.any(shape(b0prec)/=[k,k]).or. &
         any(shape(p_start)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]).or. &
         any(y<0).or.nstore<=0.or.tune<=0.0_dp)then
         res%status=1;return
      end if
      allocate(res%draws(nstore,ns*k+ns*ns),res%states(nstore,n),res%prob_state(n,ns));res%prob_state=0.0_dp
      beta=beta_start;p=p_start;keep=0;accepts=0;tries=0
      do iter=0,burnin+mcmc-1
         do j=1,ns;do i=1,n
            eta=dot_product(x(i,:),beta(j,:));if(eta>700.0_dp)then;loge(i,j)=-huge(1.0_dp);else
               loge(i,j)=real(y(i),dp)*eta-exp(eta)-log_gamma(real(y(i)+1,dp));end if
         end do;end do
         call ordered_state_sample(loge,p,s,ps,info);if(info/=0)then;res%status=10+info;return;end if
         do j=1,ns
            nj=count(s==j);nstate(j)=nj;h=0
            do i=1,n;if(s(i)==j)then;h=h+1;idx(h)=i;end if;end do
            allocate(xj(nj,k),yj(nj),yjreal(nj),prec(k,k),cov(k,k),zero(k),z(k),can(k))
            do h=1,nj;xj(h,:)=x(idx(h),:);yj(h)=y(idx(h));end do;yjreal=real(yj,dp)
            prec=b0prec+matmul(transpose(xj),xj);call inv_spd(prec,cov,info);if(info/=0)then;res%status=20+info;return;end if
            cov=tune*tune*cov;zero=0.0_dp;call rmvnorm(zero,cov,z,info);if(info/=0)then;res%status=30+info;return;end if
            can=beta(j,:)+z;cur=poisson_state_logpost(yj,xj,beta(j,:),b0,b0prec);canlp=poisson_state_logpost(yj,xj,can,b0,b0prec)
            tries=tries+1;if(log(max(runif(),tiny(1.0_dp)))<min(0.0_dp,canlp-cur))then;beta(j,:)=can;accepts=accepts+1;end if
            deallocate(xj,yj,yjreal,prec,cov,zero,z,can)
         end do
         call update_ordered_p(p,a0,nstate,info);if(info/=0)then;res%status=40+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;h=0;do j=1,ns;res%draws(keep,h+1:h+k)=beta(j,:);h=h+k;end do;call flatten_p(res%draws(keep,:),h,p)
            res%states(keep,:)=s;res%prob_state=res%prob_state+ps/real(nstore,dp)
         end if
      end do
   end function mcmc_poisson_change

   function mcmc_residual_break_analysis(resid,nbreak,b0,b0prec,c0,d0,a0,burnin,mcmc,thin) result(res)
      ! Intercept-only Gaussian changepoint model used by MCMCresidualBreakAnalysis.
      ! R-specific WAIC/Chib95 attributes and random perturbation are not part of
      ! this computational entry point.
      real(dp),intent(in)::resid(:),b0,b0prec,c0,d0,a0(:,:)
      integer,intent(in)::nbreak,burnin,mcmc,thin
      type(change_result)::res
      integer::n,ns,i,j
      real(dp),allocatable::x(:,:),beta0(:,:),sig20(:),p0(:,:),bvec(:),bprec(:,:)
      integer,allocatable::s0(:)
      real(dp)::mu,var0
      n=size(resid);ns=nbreak+1
      if(nbreak<0.or.n<ns.or.any(shape(a0)/=[ns,ns]).or.b0prec<0.0_dp.or.c0<=0.0_dp.or.d0<=0.0_dp)then
         res%status=1;return
      end if
      allocate(x(n,1),beta0(ns,1),sig20(ns),p0(ns,ns),bvec(1),bprec(1,1),s0(n))
      x=1.0_dp;mu=sum(resid)/real(n,dp)
      if(n>1)then
         var0=sum((resid-mu)**2)/real(n-1,dp)
      else
         var0=1.0_dp
      end if
      var0=max(var0,1.0e-8_dp);beta0(:,1)=mu;sig20=var0;bvec(1)=b0;bprec(1,1)=b0prec
      p0=0.0_dp
      do j=1,ns-1;p0(j,j)=0.9_dp;p0(j,j+1)=0.1_dp;end do
      p0(ns,ns)=1.0_dp
      do i=1,n;s0(i)=min(ns,1+(i-1)*ns/n);end do
      s0(n)=ns
      res=mcmc_regress_change(resid,x,beta0,sig20,p0,s0,bvec,bprec,c0,d0,a0,burnin,mcmc,thin)
   end function mcmc_residual_break_analysis

   pure real(dp) function nb_log_emit(yi,eta,rho) result(v)
      integer,intent(in)::yi
      real(dp),intent(in)::eta,rho
      real(dp)::lr
      if(rho<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      if(eta>log(rho))then;lr=eta+log(1.0_dp+rho*exp(-eta))
      else;lr=log(rho)+log(1.0_dp+exp(eta)/rho);end if
      v=log_gamma(rho+real(yi,dp))-log_gamma(rho)-log_gamma(real(yi+1,dp))+ &
        rho*log(rho)+real(yi,dp)*eta-(rho+real(yi,dp))*lr
   end function nb_log_emit

   real(dp) function slice_rho_state(rho,y,lambda,step,e,f,g) result(newrho)
      real(dp),intent(in)::rho,lambda(:),step,e,f,g
      integer,intent(in)::y(:)
      real(dp)::level,left,right,u,val
      integer::j,k,tries
      level=rho_nb_logcond(rho,y,lambda,e,f,g)+log(max(runif(),tiny(1.0_dp)))
      u=runif();left=max(tiny(1.0_dp),rho-step*u);right=left+step
      j=int(100.0_dp*runif());k=99-j
      do while(j>0)
         val=rho_nb_logcond(left,y,lambda,e,f,g);if(val<=level)exit
         left=max(tiny(1.0_dp),left-step);j=j-1
      end do
      do while(k>0)
         val=rho_nb_logcond(right,y,lambda,e,f,g);if(val<=level)exit
         right=right+step;k=k-1
      end do
      tries=0
      do
         newrho=left+runif()*(right-left)
         val=rho_nb_logcond(newrho,y,lambda,e,f,g)
         if(val>level)exit
         if(newrho>rho)then;right=newrho;else;left=newrho;end if
         tries=tries+1;if(tries>100000)then;newrho=rho;exit;end if
      end do
   end function slice_rho_state

   function mcmc_negbin_change(y,x,beta_start,rho_start,p_start,b0,b0prec,a0,e,f,g,rho_step, &
                               burnin,mcmc,thin,beta_tune) result(res)
      ! Same posterior target as MCMCnegbinChange. As in mcmc_negbin, the
      ! original auxiliary normal-mixture beta block is replaced by Gaussian MH.
      integer,intent(in)::y(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),rho_start(:),p_start(:,:),b0(:),b0prec(:,:),a0(:,:),e,f,g,rho_step(:)
      real(dp),intent(in),optional::beta_tune
      type(change_result)::res
      integer::n,ns,npar,nstore,iter,keep,j,i,h,nj,info
      integer::s(size(y)),nstate(size(rho_start)),idx(size(y))
      integer,allocatable::yj(:)
      real(dp)::beta(size(beta_start,1),size(beta_start,2)),rho(size(rho_start)),p(size(rho_start),size(rho_start))
      real(dp)::loge(size(y),size(rho_start)),ps(size(y),size(rho_start)),eta,tune,cur,canlp
      real(dp),allocatable::xj(:,:),lambda(:),prec(:,:),cov(:,:),zero(:),z(:),can(:)
      n=size(y);ns=size(rho_start);npar=size(beta_start,2);nstore=mcmc/thin;tune=0.20_dp;if(present(beta_tune))tune=beta_tune
      if(n<ns.or.size(x,1)/=n.or.size(x,2)/=npar.or.size(beta_start,1)/=ns.or.size(b0)/=npar.or. &
         any(shape(b0prec)/=[npar,npar]).or.any(shape(p_start)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]).or. &
         size(rho_step)/=ns.or.any(rho_start<=0.0_dp).or.any(rho_step<=0.0_dp).or.any(y<0).or. &
         nstore<=0.or.tune<=0.0_dp)then;res%status=1;return;end if
      allocate(res%draws(nstore,ns*npar+ns+ns*ns),res%states(nstore,n),res%prob_state(n,ns));res%prob_state=0.0_dp
      beta=beta_start;rho=rho_start;p=p_start;keep=0
      do iter=0,burnin+mcmc-1
         ! State path conditional on current regression and dispersion parameters.
         do j=1,ns;do i=1,n
            eta=dot_product(x(i,:),beta(j,:));loge(i,j)=nb_log_emit(y(i),eta,rho(j))
         end do;end do
         call ordered_state_sample(loge,p,s,ps,info);if(info/=0)then;res%status=10+info;return;end if

         do j=1,ns
            nj=count(s==j);nstate(j)=nj;if(nj<=0)then;res%status=2;return;end if
            h=0;do i=1,n;if(s(i)==j)then;h=h+1;idx(h)=i;end if;end do
            allocate(xj(nj,npar),yj(nj),lambda(nj),prec(npar,npar),cov(npar,npar),zero(npar),z(npar),can(npar))
            do h=1,nj;xj(h,:)=x(idx(h),:);yj(h)=y(idx(h));end do
            lambda=exp(min(700.0_dp,matmul(xj,beta(j,:))))
            rho(j)=slice_rho_state(rho(j),yj,lambda,rho_step(j),e,f,g)
            prec=b0prec+matmul(transpose(xj),xj);call inv_spd(prec,cov,info);if(info/=0)then;res%status=20+info;return;end if
            cov=tune*tune*cov;zero=0.0_dp;call rmvnorm(zero,cov,z,info);if(info/=0)then;res%status=30+info;return;end if
            can=beta(j,:)+z;cur=negbin_logpost(yj,xj,beta(j,:),rho(j),b0,b0prec,e,f,g)
            canlp=negbin_logpost(yj,xj,can,rho(j),b0,b0prec,e,f,g)
            if(log(max(runif(),tiny(1.0_dp)))<min(0.0_dp,canlp-cur))beta(j,:)=can
            deallocate(xj,yj,lambda,prec,cov,zero,z,can)
         end do
         call update_ordered_p(p,a0,nstate,info);if(info/=0)then;res%status=40+info;return;end if

         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;h=0
            do j=1,ns;res%draws(keep,h+1:h+npar)=beta(j,:);h=h+npar;end do
            res%draws(keep,h+1:h+ns)=rho;h=h+ns;call flatten_p(res%draws(keep,:),h,p)
            res%states(keep,:)=s;res%prob_state=res%prob_state+ps/real(nstore,dp)
         end if
      end do
   end function mcmc_negbin_change
end module mcmcpack_changepoint
