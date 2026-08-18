! SPDX-License-Identifier: GPL-3.0-only
! Hierarchical logit and Poisson models translated from cMCMChlogit.cc and
! cMCMChpoisson.cc.  The latent linear-predictor MH step follows the original;
! Gaussian fixed/random effects are updated with algebraically equivalent
! conditional Gibbs steps.
module mcmcpack_hglm
   use mcmcpack_kinds, only : dp,pi
   use mcmcpack_math, only : logistic,log1pexp,log_factorial
   use mcmcpack_rng, only : runif,rnorm,rmvnorm,rinvgamma_rng
   use mcmcpack_distributions, only : riwish
   use mcmcpack_linalg, only : inv_spd
   implicit none
   private
   type, public :: hglm_result
      real(dp), allocatable :: draws(:,:)
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: latent_accept_rate(:)
      integer :: status=0
   end type hglm_result
   public :: mcmc_hlogit,mcmc_hpoisson
contains
   function mcmc_hlogit(y,x,w,group,beta_start,b_start,vb_start,v_start,latent_start, &
                        mubeta,vbeta,r,rmat,s1,s2,burnin,mcmc,thin,fix_od) result(res)
      integer,intent(in)::y(:),group(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),w(:,:),beta_start(:),b_start(:,:),vb_start(:,:),v_start,latent_start(:)
      real(dp),intent(in)::mubeta(:),vbeta(:,:),r,rmat(:,:),s1,s2
      logical,intent(in),optional::fix_od
      type(hglm_result)::res
      res=hglm_engine(1,y,x,w,group,beta_start,b_start,vb_start,v_start,latent_start, &
                      mubeta,vbeta,r,rmat,s1,s2,burnin,mcmc,thin,fix_od)
   end function mcmc_hlogit

   function mcmc_hpoisson(y,x,w,group,beta_start,b_start,vb_start,v_start,latent_start, &
                          mubeta,vbeta,r,rmat,s1,s2,burnin,mcmc,thin,fix_od) result(res)
      integer,intent(in)::y(:),group(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),w(:,:),beta_start(:),b_start(:,:),vb_start(:,:),v_start,latent_start(:)
      real(dp),intent(in)::mubeta(:),vbeta(:,:),r,rmat(:,:),s1,s2
      logical,intent(in),optional::fix_od
      type(hglm_result)::res
      res=hglm_engine(2,y,x,w,group,beta_start,b_start,vb_start,v_start,latent_start, &
                      mubeta,vbeta,r,rmat,s1,s2,burnin,mcmc,thin,fix_od)
   end function mcmc_hpoisson

   function hglm_engine(family,y,x,w,group,beta_start,b_start,vb_start,v_start,latent_start, &
                        mubeta,vbeta,r,rmat,s1,s2,burnin,mcmc,thin,fix_od) result(res)
      integer,intent(in)::family,y(:),group(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),w(:,:),beta_start(:),b_start(:,:),vb_start(:,:),v_start,latent_start(:)
      real(dp),intent(in)::mubeta(:),vbeta(:,:),r,rmat(:,:),s1,s2
      logical,intent(in),optional::fix_od
      type(hglm_result)::res
      integer::n,p,q,ng,nsamp,iter,keep,i,j,g,k,info,ncol,adapt_every
      real(dp),allocatable::beta(:),b(:,:),vb(:,:),ivbeta(:,:),ivb(:,:),z(:),proposal_sd(:),accepted(:),trials(:)
      real(dp),allocatable::fit(:),offset(:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:),ssb(:,:)
      real(dp)::v,zp,eta,lpold,lpnew,sse,dev,cconst,rate
      logical::fixedv
      n=size(y);p=size(x,2);q=size(w,2);ng=size(b_start,1);nsamp=mcmc/thin
      fixedv=.false.;if(present(fix_od))fixedv=fix_od
      if(thin<=0.or.mcmc<=0.or.nsamp<=0.or.size(x,1)/=n.or.size(w,1)/=n.or.size(group)/=n.or. &
         size(beta_start)/=p.or.size(mubeta)/=p.or.any(shape(vbeta)/=[p,p]).or.size(b_start,2)/=q.or. &
         any(shape(vb_start)/=[q,q]).or.any(shape(rmat)/=[q,q]).or.size(latent_start)/=n.or.v_start<=0.0_dp)then
         res%status=1;return
      end if
      if(any(group<1).or.any(group>ng).or.any(y<0))then;res%status=2;return;end if
      if(family==1.and.any(y>1))then;res%status=3;return;end if
      allocate(beta(p),b(ng,q),vb(q,q),ivbeta(p,p),ivb(q,q),z(n),proposal_sd(n),accepted(n),trials(n), &
               fit(n),offset(n),ssb(q,q))
      beta=beta_start;b=b_start;vb=vb_start;v=v_start;z=latent_start;proposal_sd=1.0_dp;accepted=0.0_dp;trials=0.0_dp
      call inv_spd(vbeta,ivbeta,info);if(info/=0)then;res%status=10+info;return;end if
      ncol=p+ng*q+q*q+2
      allocate(res%draws(nsamp,ncol),res%prediction(n),res%latent_accept_rate(n));res%prediction=0.0_dp;keep=0
      adapt_every=merge(100,max(1,(burnin+mcmc)/10),burnin+mcmc>=1000)
      cconst=(16.0_dp*sqrt(3.0_dp)/(15.0_dp*pi))**2

      do iter=1,burnin+mcmc
         ! Latent transformed response, random-walk MH.
         do i=1,n
            eta=dot_product(x(i,:),beta)+dot_product(w(i,:),b(group(i),:))
            zp=rnorm(z(i),proposal_sd(i))
            if(family==1)then
               lpold=real(y(i),dp)*z(i)-log1pexp(z(i))-0.5_dp*(z(i)-eta)**2/v
               lpnew=real(y(i),dp)*zp-log1pexp(zp)-0.5_dp*(zp-eta)**2/v
            else
               lpold=real(y(i),dp)*z(i)-exp(min(z(i),log(huge(1.0_dp))))-0.5_dp*(z(i)-eta)**2/v
               lpnew=real(y(i),dp)*zp-exp(min(zp,log(huge(1.0_dp))))-0.5_dp*(zp-eta)**2/v
            end if
            trials(i)=trials(i)+1.0_dp
            if(log(runif())<min(0.0_dp,lpnew-lpold))then;z(i)=zp;accepted(i)=accepted(i)+1.0_dp;end if
         end do

         ! beta | z,b,V
         allocate(prec(p,p),cov(p,p),rhs(p),mu(p),draw(p))
         do i=1,n;offset(i)=dot_product(w(i,:),b(group(i),:));end do
         prec=ivbeta+matmul(transpose(x),x)/v
         rhs=matmul(ivbeta,mubeta)+matmul(transpose(x),z-offset)/v
         call inv_spd(prec,cov,info);if(info/=0)then;res%status=20+info;return;end if
         mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=30+info;return;end if
         beta=draw;deallocate(prec,cov,rhs,mu,draw)

         ! b_g | z,beta,V,Vb
         call inv_spd(vb,ivb,info);if(info/=0)then;res%status=40+info;return;end if
         do g=1,ng
            allocate(prec(q,q),cov(q,q),rhs(q),mu(q),draw(q));prec=ivb;rhs=0.0_dp
            do i=1,n
               if(group(i)/=g)cycle
               do j=1,q;do k=1,q;prec(j,k)=prec(j,k)+w(i,j)*w(i,k)/v;end do;end do
               rhs=rhs+w(i,:)*(z(i)-dot_product(x(i,:),beta))/v
            end do
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=50+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=60+info;return;end if
            b(g,:)=draw;deallocate(prec,cov,rhs,mu,draw)
         end do

         ssb=0.0_dp
         do g=1,ng;do j=1,q;do k=1,q;ssb(j,k)=ssb(j,k)+b(g,j)*b(g,k);end do;end do;end do
         call riwish(r+real(ng,dp),ssb+r*rmat,vb,info);if(info/=0)then;res%status=70+info;return;end if

         do i=1,n;fit(i)=dot_product(x(i,:),beta)+dot_product(w(i,:),b(group(i),:));end do
         if(.not.fixedv)then
            sse=dot_product(z-fit,z-fit);v=rinvgamma_rng(s1+0.5_dp*real(n,dp),s2+0.5_dp*sse)
         end if
         dev=0.0_dp
         if(family==1)then
            do i=1,n;dev=dev-2.0_dp*(real(y(i),dp)*z(i)-log1pexp(z(i)));end do
         else
            do i=1,n;dev=dev-2.0_dp*(real(y(i),dp)*z(i)-exp(min(z(i),log(huge(1.0_dp))))-log_factorial(y(i)));end do
         end if

         if(iter>burnin.and.mod(iter-burnin,thin)==0)then
            keep=keep+1;k=0;res%draws(keep,1:p)=beta;k=p
            do g=1,ng;res%draws(keep,k+1:k+q)=b(g,:);k=k+q;end do
            do j=1,q;do g=1,q;k=k+1;res%draws(keep,k)=vb(j,g);end do;end do
            res%draws(keep,k+1)=v;res%draws(keep,k+2)=dev
            if(family==1)then
               do i=1,n;res%prediction(i)=res%prediction(i)+logistic(fit(i)/sqrt(1.0_dp+cconst*v))/real(nsamp,dp);end do
            else
               do i=1,n;res%prediction(i)=res%prediction(i)+exp(min(fit(i)+0.5_dp*v,log(huge(1.0_dp))))/real(nsamp,dp);end do
            end if
         end if

         if(iter<=burnin.and.mod(iter,adapt_every)==0)then
            do i=1,n
               if(trials(i)<=0.0_dp)cycle
               rate=accepted(i)/trials(i)
               if(rate>=0.44_dp)then;proposal_sd(i)=proposal_sd(i)*(2.0_dp-(1.0_dp-rate)/(1.0_dp-0.44_dp))
               else;proposal_sd(i)=proposal_sd(i)/(2.0_dp-rate/0.44_dp);end if
               accepted(i)=0.0_dp;trials(i)=0.0_dp
            end do
         end if
      end do
      do i=1,n
         ! If the last adaptation reset counters, this is the post-adaptation rate;
         ! otherwise it is the rate over all iterations since the last reset.
         if(trials(i)>0.0_dp)then;res%latent_accept_rate(i)=accepted(i)/trials(i);else;res%latent_accept_rate(i)=0.0_dp;end if
      end do
   end function hglm_engine
end module mcmcpack_hglm
