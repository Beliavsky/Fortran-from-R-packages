! SPDX-License-Identifier: GPL-3.0-only
! Gaussian factor-analysis Gibbs sampler translated from cMCMCfactanal.cc/MCMCfcds.h.
module mcmcpack_factor
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : rmvnorm,rinvgamma_rng
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_samplers, only : mcmc_result
   implicit none
   private
   public :: mcmc_factanal
contains
   function mcmc_factanal(x,lambda_start,psi_start,lambda_eq,lambda_ineq,lambda_prior_mean,lambda_prior_prec, &
                          a0,b0,burnin,mcmc,thin,store_scores) result(res)
      real(dp),intent(in)::x(:,:),lambda_start(:,:),psi_start(:),lambda_eq(:,:),lambda_ineq(:,:)
      real(dp),intent(in)::lambda_prior_mean(:,:),lambda_prior_prec(:,:),a0(:),b0(:)
      integer,intent(in)::burnin,mcmc,thin
      logical,intent(in),optional::store_scores
      type(mcmc_result)::res
      integer::n,k,d,iter,store_count,nsamp,info,i,j,h,nfree,col,tries
      logical::ss,ok
      integer,allocatable::idx(:)
      real(dp)::lambda(size(lambda_start,1),size(lambda_start,2)),psi(size(psi_start))
      real(dp)::phi(size(x,1),size(lambda_start,2)),psiinv(size(psi_start))
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mean(:),draw(:),pf(:,:),r(:)
      real(dp)::phiprec(size(lambda_start,2),size(lambda_start,2)),phicov(size(lambda_start,2),size(lambda_start,2))
      real(dp)::phimean(size(lambda_start,2)),tmp(size(lambda_start,2)),e(size(x,1)),sse
      n=size(x,1);k=size(x,2);d=size(lambda_start,2);nsamp=mcmc/thin;ss=.false.;if(present(store_scores))ss=store_scores
      if(size(lambda_start,1)/=k.or.size(psi_start)/=k.or.any(shape(lambda_eq)/=[k,d]).or.any(shape(lambda_ineq)/=[k,d]).or. &
         any(shape(lambda_prior_mean)/=[k,d]).or.any(shape(lambda_prior_prec)/=[k,d]).or.size(a0)/=k.or.size(b0)/=k.or.nsamp<=0)then
         res%status=1;return
      end if
      allocate(res%draws(nsamp,k*d+k+merge(n*d,0,ss)));lambda=lambda_start;psi=psi_start;phi=0.0_dp;store_count=0
      do iter=0,burnin+mcmc-1
         psiinv=1.0_dp/psi
         ! Factor scores: precision I + Lambda' Psi^-1 Lambda.
         phiprec=0.0_dp;do j=1,d;phiprec(j,j)=1.0_dp;end do
         do i=1,k
            do j=1,d;do h=1,d;phiprec(j,h)=phiprec(j,h)+psiinv(i)*lambda(i,j)*lambda(i,h);end do;end do
         end do
         call inv_spd(phiprec,phicov,info);if(info/=0)then;res%status=10+info;return;end if
         do i=1,n
            tmp=0.0_dp
            do j=1,d;do h=1,k;tmp(j)=tmp(j)+lambda(h,j)*psiinv(h)*x(i,h);end do;end do
            phimean=matmul(phicov,tmp);call rmvnorm(phimean,phicov,phi(i,:),info);if(info/=0)then;res%status=20+info;return;end if
         end do
         ! Loadings row by row, conditioning on fixed loadings.
         do i=1,k
            nfree=count(lambda_eq(i,:)<-998.5_dp);if(nfree==0)cycle
            allocate(idx(nfree),prec(nfree,nfree),cov(nfree,nfree),rhs(nfree),mean(nfree),draw(nfree),pf(n,nfree),r(n))
            h=0;do j=1,d;if(lambda_eq(i,j)<-998.5_dp)then;h=h+1;idx(h)=j;end if;end do
            r=x(:,i)
            do j=1,d
               if(lambda_eq(i,j)>=-998.5_dp)r=r-phi(:,j)*lambda(i,j)
            end do
            do h=1,nfree;pf(:,h)=phi(:,idx(h));end do
            prec=psiinv(i)*matmul(transpose(pf),pf);rhs=psiinv(i)*matmul(transpose(pf),r)
            do h=1,nfree
               prec(h,h)=prec(h,h)+lambda_prior_prec(i,idx(h));rhs(h)=rhs(h)+lambda_prior_prec(i,idx(h))*lambda_prior_mean(i,idx(h))
            end do
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=30+info;return;end if;mean=matmul(cov,rhs)
            tries=0
            do
               call rmvnorm(mean,cov,draw,info);if(info/=0)then;res%status=40+info;return;end if
               ok=.true.;do h=1,nfree
                  if(lambda_ineq(i,idx(h))*draw(h)<0.0_dp)ok=.false.
               end do
               tries=tries+1;if(ok.or.tries>100000)exit
            end do
            if(.not.ok)then;res%status=50;return;end if
            do h=1,nfree;lambda(i,idx(h))=draw(h);end do
            deallocate(idx,prec,cov,rhs,mean,draw,pf,r)
         end do
         ! Uniqueness variances.
         do i=1,k
            e=x(:,i)-matmul(phi,lambda(i,:));sse=dot_product(e,e)
            psi(i)=rinvgamma_rng(0.5_dp*(a0(i)+real(n,dp)),0.5_dp*(b0(i)+sse))
         end do
         if(iter>=burnin.and.mod(iter,thin)==0)then
            store_count=store_count+1;col=0
            do i=1,k;do j=1,d;col=col+1;res%draws(store_count,col)=lambda(i,j);end do;end do
            res%draws(store_count,col+1:col+k)=psi;col=col+k
            if(ss)then;do i=1,n;do j=1,d;col=col+1;res%draws(store_count,col)=phi(i,j);end do;end do;end if
         end if
      end do
   end function mcmc_factanal
end module mcmcpack_factor
