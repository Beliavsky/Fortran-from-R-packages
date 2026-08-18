! SPDX-License-Identifier: GPL-3.0-only
! Mixed continuous/ordinal factor-analysis sampler translated from
! MCMCmixfactanal.cc.
module mcmcpack_mixfactor
   use mcmcpack_kinds, only : dp
   use mcmcpack_math, only : normal_cdf
   use mcmcpack_rng, only : runif,rnorm,rtruncnorm,rmvnorm,rinvgamma_rng
   use mcmcpack_linalg, only : inv_spd
   implicit none
   private
   type, public :: mixfactor_result
      real(dp),allocatable::draws(:,:)
      real(dp),allocatable::threshold_accept_rate(:)
      integer::status=0
   end type
   public :: mcmc_mixfactanal
contains
   function mcmc_mixfactanal(x,ncat,lambda_start,gamma_start,psi_start,lambda_eq,lambda_ineq, &
                             lambda_prior_mean,lambda_prior_prec,a0,b0,tune,burnin,mcmc,thin, &
                             store_lambda,store_scores) result(res)
      real(dp),intent(in)::x(:,:),lambda_start(:,:),gamma_start(:,:),psi_start(:),lambda_eq(:,:),lambda_ineq(:,:), &
                           lambda_prior_mean(:,:),lambda_prior_prec(:,:),a0(:),b0(:),tune(:)
      integer,intent(in)::ncat(:),burnin,mcmc,thin
      logical,intent(in),optional::store_lambda,store_scores
      type(mixfactor_result)::res
      integer::n,k,d,gdim,nsamp,iter,keep,i,j,h,c,nfree,info,col,tries,cat
      integer,allocatable::idx(:),accepts(:)
      real(dp)::lambda(size(lambda_start,1),size(lambda_start,2)),gamma(size(gamma_start,1),size(gamma_start,2))
      real(dp)::psi(size(psi_start)),phi(size(x,1),size(lambda_start,2)),xstar(size(x,1),size(x,2)),xmean(size(x,2))
      real(dp)::precphi(max(1,size(lambda_start,2)-1),max(1,size(lambda_start,2)-1))
      real(dp)::covphi(max(1,size(lambda_start,2)-1),max(1,size(lambda_start,2)-1))
      real(dp)::rhsphi(max(1,size(lambda_start,2)-1)),meanphi(max(1,size(lambda_start,2)-1)),drawphi(max(1,size(lambda_start,2)-1))
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mean(:),draw(:),pf(:,:),rr(:),gp(:),eta(:)
      real(dp)::sse,llr,qratio,pnew,pold,nf,nr
      logical::sl,ss,ok
      sl=.true.;ss=.false.;if(present(store_lambda))sl=store_lambda;if(present(store_scores))ss=store_scores
      n=size(x,1);k=size(x,2);d=size(lambda_start,2);gdim=size(gamma_start,1);nsamp=mcmc/thin
      if(thin<=0.or.nsamp<=0.or.size(ncat)/=k.or.size(lambda_start,1)/=k.or.size(gamma_start,2)/=k.or. &
         size(psi_start)/=k.or.any(psi_start<=0.0_dp).or.any(shape(lambda_eq)/=[k,d]).or. &
         any(shape(lambda_ineq)/=[k,d]).or.any(shape(lambda_prior_mean)/=[k,d]).or.any(shape(lambda_prior_prec)/=[k,d]).or. &
         size(a0)/=k.or.size(b0)/=k.or.size(tune)/=k.or.maxval(ncat)+1>gdim.or.d<2)then;res%status=1;return;end if
      allocate(accepts(k));accepts=0
      col=gdim*k+k+merge(k*d,0,sl)+merge(n*d,0,ss)
      allocate(res%draws(nsamp,col),res%threshold_accept_rate(k))
      lambda=lambda_start;gamma=gamma_start;psi=psi_start;phi=0.0_dp;phi(:,1)=1.0_dp;xstar=x;keep=0
      do j=1,k;if(ncat(j)>=2)psi(j)=1.0_dp;end do

      do iter=0,burnin+mcmc-1
         do i=1,n
            xmean=matmul(lambda,phi(i,:))
            do j=1,k
               if(ncat(j)>=2)then
                  if(nint(x(i,j))==-999)then;xstar(i,j)=rnorm(xmean(j),1.0_dp)
                  else;cat=nint(x(i,j));if(cat<1.or.cat>ncat(j))then;res%status=2;return;end if
                     xstar(i,j)=rtruncnorm(xmean(j),1.0_dp,gamma(cat,j),gamma(cat+1,j));end if
               else if(nint(x(i,j))==-999)then
                  xstar(i,j)=rnorm(xmean(j),sqrt(psi(j)))
               else
                  xstar(i,j)=x(i,j)
               end if
            end do
         end do

         precphi=0.0_dp;do h=1,d-1;precphi(h,h)=1.0_dp;end do
         do j=1,k;do h=1,d-1;do c=1,d-1
            precphi(h,c)=precphi(h,c)+lambda(j,h+1)*lambda(j,c+1)/psi(j)
         end do;end do;end do
         call inv_spd(precphi(1:d-1,1:d-1),covphi(1:d-1,1:d-1),info);if(info/=0)then;res%status=10+info;return;end if
         do i=1,n
            rhsphi(1:d-1)=0.0_dp
            do j=1,k;rhsphi(1:d-1)=rhsphi(1:d-1)+lambda(j,2:d)*(xstar(i,j)-lambda(j,1))/psi(j);end do
            meanphi(1:d-1)=matmul(covphi(1:d-1,1:d-1),rhsphi(1:d-1))
            call rmvnorm(meanphi(1:d-1),covphi(1:d-1,1:d-1),drawphi(1:d-1),info);if(info/=0)then;res%status=20+info;return;end if
            phi(i,2:d)=drawphi(1:d-1)
         end do

         do j=1,k
            nfree=count(lambda_eq(j,:)<-998.5_dp);if(nfree==0)then;lambda(j,:)=lambda_eq(j,:);cycle;end if
            allocate(idx(nfree),prec(nfree,nfree),cov(nfree,nfree),rhs(nfree),mean(nfree),draw(nfree),pf(n,nfree),rr(n))
            h=0;do c=1,d;if(lambda_eq(j,c)<-998.5_dp)then;h=h+1;idx(h)=c;else;lambda(j,c)=lambda_eq(j,c);end if;end do
            rr=xstar(:,j);do c=1,d;if(lambda_eq(j,c)>=-998.5_dp)rr=rr-phi(:,c)*lambda(j,c);end do
            do h=1,nfree;pf(:,h)=phi(:,idx(h));end do
            prec=matmul(transpose(pf),pf)/psi(j);rhs=matmul(transpose(pf),rr)/psi(j)
            do h=1,nfree;prec(h,h)=prec(h,h)+lambda_prior_prec(j,idx(h)); &
               rhs(h)=rhs(h)+lambda_prior_prec(j,idx(h))*lambda_prior_mean(j,idx(h));end do
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=30+info;return;end if;mean=matmul(cov,rhs);tries=0
            do
               call rmvnorm(mean,cov,draw,info);if(info/=0)then;res%status=40+info;return;end if;ok=.true.
               do h=1,nfree;if(lambda_ineq(j,idx(h))*draw(h)<0.0_dp)ok=.false.;end do
               tries=tries+1;if(ok.or.tries>100000)exit
            end do
            if(.not.ok)then;res%status=50;return;end if
            do h=1,nfree;lambda(j,idx(h))=draw(h);end do
            deallocate(idx,prec,cov,rhs,mean,draw,pf,rr)
         end do

         do j=1,k
            if(ncat(j)<2)then
               rr_alloc:block
                  real(dp)::evec(n)
                  evec=xstar(:,j)-matmul(phi,lambda(j,:));sse=dot_product(evec,evec)
                  psi(j)=rinvgamma_rng(0.5_dp*(a0(j)+real(n,dp)),0.5_dp*(b0(j)+sse))
               end block rr_alloc
            else;psi(j)=1.0_dp;end if
         end do

         do j=1,k
            if(ncat(j)<=2)then;if(ncat(j)==2)accepts(j)=accepts(j)+1;cycle;end if
            if(tune(j)<=0.0_dp)cycle
            allocate(gp(gdim),eta(n));gp=gamma(:,j);eta=matmul(phi,lambda(j,:))
            do c=3,ncat(j)
               if(c==ncat(j))then;gp(c)=rtruncnorm(gamma(c,j),tune(j),gp(c-1),huge(1.0_dp)/10.0_dp)
               else;gp(c)=rtruncnorm(gamma(c,j),tune(j),gp(c-1),gamma(c+1,j));end if
            end do
            llr=0.0_dp
            do i=1,n
               if(nint(x(i,j))==-999)cycle;cat=nint(x(i,j))
               pnew=max(tiny(1.0_dp),ord_prob(cat,ncat(j),eta(i),gp));pold=max(tiny(1.0_dp),ord_prob(cat,ncat(j),eta(i),gamma(:,j)))
               llr=llr+log(pnew)-log(pold)
            end do
            qratio=0.0_dp
            do c=3,ncat(j)
               nf=normal_cdf((gamma(c+1,j)-gamma(c,j))/tune(j))-normal_cdf((gp(c-1)-gamma(c,j))/tune(j))
               nr=normal_cdf((gp(c+1)-gp(c))/tune(j))-normal_cdf((gamma(c-1,j)-gp(c))/tune(j))
               qratio=qratio+log(max(nf,tiny(1.0_dp)))-log(max(nr,tiny(1.0_dp)))
            end do
            if(log(runif())<min(0.0_dp,llr+qratio))then;gamma(1:ncat(j)+1,j)=gp(1:ncat(j)+1);accepts(j)=accepts(j)+1;end if
            deallocate(gp,eta)
         end do

         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            if(sl)then;do j=1,k;do c=1,d;col=col+1;res%draws(keep,col)=lambda(j,c);end do;end do;end if
            do j=1,k;do c=1,gdim;col=col+1;res%draws(keep,col)=gamma(c,j);end do;end do
            if(ss)then;do i=1,n;do c=1,d;col=col+1;res%draws(keep,col)=phi(i,c);end do;end do;end if
            res%draws(keep,col+1:col+k)=psi
         end if
      end do
      res%threshold_accept_rate=real(accepts,dp)/real(max(1,burnin+mcmc),dp)
   end function

   pure real(dp) function ord_prob(category,ncategory,eta,g) result(p)
      integer,intent(in)::category,ncategory
      real(dp),intent(in)::eta,g(:)
      if(category==1)then;p=normal_cdf(g(2)-eta)
      else if(category==ncategory)then;p=1.0_dp-normal_cdf(g(ncategory)-eta)
      else;p=normal_cdf(g(category+1)-eta)-normal_cdf(g(category)-eta);end if
      p=max(0.0_dp,p)
   end function
end module mcmcpack_mixfactor
