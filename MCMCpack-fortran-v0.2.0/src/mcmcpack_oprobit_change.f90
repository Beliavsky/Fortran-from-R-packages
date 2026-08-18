! SPDX-License-Identifier: GPL-3.0-only
! Ordered-probit changepoint sampler translated from cMCMCoprobitChange.cc.
module mcmcpack_oprobit_change
   use mcmcpack_kinds, only : dp,pi
   use mcmcpack_math, only : normal_cdf
   use mcmcpack_rng, only : runif,rtruncnorm,rmvnorm,rbeta
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_changepoint, only : ordered_state_sample
   implicit none
   private
   type, public :: oprobit_change_result
      real(dp),allocatable::draws(:,:)
      integer,allocatable::states(:,:)
      real(dp),allocatable::prob_state(:,:)
      real(dp),allocatable::gamma_accept_rate(:)
      integer::status=0
   end type
   public :: mcmc_oprobit_change
contains
   pure real(dp) function ordprob(cat,ncat,eta,gam) result(p)
      integer,intent(in)::cat,ncat
      real(dp),intent(in)::eta,gam(:)
      if(cat==1)then;p=normal_cdf(gam(2)-eta)
      else if(cat==ncat)then;p=1.0_dp-normal_cdf(gam(ncat)-eta)
      else;p=normal_cdf(gam(cat+1)-eta)-normal_cdf(gam(cat)-eta);end if
      p=max(tiny(1.0_dp),p)
   end function

   subroutine sample_transition(p,a0,nstate,status)
      real(dp),intent(inout)::p(:,:)
      real(dp),intent(in)::a0(:,:)
      integer,intent(in)::nstate(:)
      integer,intent(out)::status
      integer::j,ns
      real(dp)::aa,bb
      ns=size(p,1);status=0;p=0.0_dp;p(ns,ns)=1.0_dp
      do j=1,ns-1
         aa=a0(j,j)+real(nstate(j)-1,dp);bb=a0(j,j+1)+1.0_dp
         if(aa<=0.0_dp.or.bb<=0.0_dp)then;status=1;return;end if
         p(j,j)=rbeta(aa,bb);p(j,j+1)=1.0_dp-p(j,j)
      end do
   end subroutine

   function mcmc_oprobit_change(y,x,beta_start,beta_linear_start,gamma_start,p_start,linear_variance, &
                                b0,b0prec,a0,tune,burnin,mcmc,thin,gamma_fixed) result(res)
      integer,intent(in)::y(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:,:),beta_linear_start(:,:),gamma_start(:,:),p_start(:,:), &
                           linear_variance,b0(:),b0prec(:,:),a0(:,:),tune(:)
      logical,intent(in),optional::gamma_fixed
      type(oprobit_change_result)::res
      integer::n,ns,k,ncat,gk,nstore,iter,keep,i,j,h,nj,info,col,c
      integer::s(size(y)),nstate(size(beta_start,1)),idx(size(y))
      real(dp)::beta(size(beta_start,1),size(beta_start,2)),blin(size(beta_linear_start,1),size(beta_linear_start,2))
      real(dp),allocatable::gamma(:,:),xj(:,:),zj(:),yjreal(:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:)
      real(dp)::p(size(beta_start,1),size(beta_start,1)),loge(size(y),size(beta_start,1))
      real(dp)::ps(size(y),size(beta_start,1)),z(size(y)),eta
      integer,allocatable::accepts(:)
      logical::gf
      gf=.false.;if(present(gamma_fixed))gf=gamma_fixed
      n=size(y);ns=size(beta_start,1);k=size(beta_start,2);ncat=maxval(y);gk=ncat+1;nstore=mcmc/thin
      if(thin<=0.or.nstore<=0.or.n<ns.or.size(x,1)/=n.or.size(x,2)/=k.or.any(y<1).or. &
         size(beta_linear_start,1)/=ns.or.size(beta_linear_start,2)/=k.or.size(b0)/=k.or. &
         any(shape(b0prec)/=[k,k]).or.any(shape(p_start)/=[ns,ns]).or.any(shape(a0)/=[ns,ns]).or. &
         size(tune)/=ns.or.any(tune<0.0_dp).or.linear_variance<=0.0_dp)then;res%status=1;return;end if
      if(gf)then
         if(size(gamma_start,1)/=1.or.size(gamma_start,2)<gk)then;res%status=2;return;end if
         allocate(gamma(1,gk));gamma=gamma_start(:,1:gk)
      else
         if(size(gamma_start,1)/=ns.or.size(gamma_start,2)<gk)then;res%status=2;return;end if
         allocate(gamma(ns,gk));gamma=gamma_start(:,1:gk)
      end if
      if(any(abs(gamma(:,1)+300.0_dp)>1e-8_dp).or.any(abs(gamma(:,2))>1e-8_dp))then;res%status=3;return;end if
      allocate(accepts(ns));accepts=0
      col=2*ns*k+merge(gk,ns*gk,gf)+ns*ns
      allocate(res%draws(nstore,col),res%states(nstore,n),res%prob_state(n,ns),res%gamma_accept_rate(ns));res%prob_state=0.0_dp
      beta=beta_start;blin=beta_linear_start;p=p_start;keep=0
      do iter=0,burnin+mcmc-1
         ! State draw uses MCMCpack's Gaussian linear approximation to the ordinal response.
         do j=1,ns;do i=1,n
            eta=dot_product(x(i,:),blin(j,:));loge(i,j)=-0.5_dp*(log(2.0_dp*pi*linear_variance)+ &
                         (real(y(i),dp)-eta)**2/linear_variance)
         end do;end do
         call ordered_state_sample(loge,p,s,ps,info);if(info/=0)then;res%status=10+info;return;end if

         ! Latent ordinal utilities.
         do i=1,n
            eta=dot_product(x(i,:),beta(s(i),:));c=y(i)
            if(gf)then;z(i)=rtruncnorm(eta,1.0_dp,gamma(1,c),gamma(1,c+1))
            else;z(i)=rtruncnorm(eta,1.0_dp,gamma(s(i),c),gamma(s(i),c+1));end if
         end do

         ! State-specific exact latent-probit beta and Gaussian-approximation beta.
         do j=1,ns
            nj=count(s==j);nstate(j)=nj;if(nj<=0)then;res%status=4;return;end if
            h=0;do i=1,n;if(s(i)==j)then;h=h+1;idx(h)=i;end if;end do
            allocate(xj(nj,k),zj(nj),yjreal(nj),prec(k,k),cov(k,k),rhs(k),mu(k),draw(k))
            do h=1,nj;xj(h,:)=x(idx(h),:);zj(h)=z(idx(h));yjreal(h)=real(y(idx(h)),dp);end do
            prec=b0prec+matmul(transpose(xj),xj);rhs=matmul(b0prec,b0)+matmul(transpose(xj),zj)
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=20+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=30+info;return;end if;beta(j,:)=draw
            prec=b0prec+matmul(transpose(xj),xj)/linear_variance
            rhs=matmul(b0prec,b0)+matmul(transpose(xj),yjreal)/linear_variance
            call inv_spd(prec,cov,info);if(info/=0)then;res%status=40+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=50+info;return;end if;blin(j,:)=draw
            deallocate(xj,zj,yjreal,prec,cov,rhs,mu,draw)
         end do

         ! Threshold Metropolis updates. For gamma_fixed a common threshold vector
         ! is updated against all observations; otherwise each regime has its own.
         if(gf)then
            call update_gamma_row(1,.true.,y,x,beta,s,gamma,tune(1),ncat,accepts(1))
         else
            do j=1,ns;call update_gamma_row(j,.false.,y,x,beta,s,gamma,tune(j),ncat,accepts(j));end do
         end if

         call sample_transition(p,a0,nstate,info);if(info/=0)then;res%status=60+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            do j=1,ns;res%draws(keep,col+1:col+k)=beta(j,:);col=col+k;end do
            do j=1,ns;res%draws(keep,col+1:col+k)=blin(j,:);col=col+k;end do
            if(gf)then;res%draws(keep,col+1:col+gk)=gamma(1,:);col=col+gk
            else;do j=1,ns;res%draws(keep,col+1:col+gk)=gamma(j,:);col=col+gk;end do;end if
            do i=1,ns;do j=1,ns;col=col+1;res%draws(keep,col)=p(i,j);end do;end do
            res%states(keep,:)=s;res%prob_state=res%prob_state+ps/real(nstore,dp)
         end if
      end do
      res%gamma_accept_rate=real(accepts,dp)/real(max(1,burnin+mcmc),dp)
   end function

   subroutine update_gamma_row(regime,shared,y,x,beta,s,gamma,tune,ncat,accept_count)
      integer,intent(in)::regime,y(:),s(:),ncat
      real(dp),intent(in)::x(:,:),beta(:,:),tune
      logical,intent(in)::shared
      real(dp),intent(inout)::gamma(:,:)
      integer,intent(inout)::accept_count
      real(dp),allocatable::gp(:)
      integer::c,i,rrow
      real(dp)::eta,llr,qratio,pnew,pold,nf,nr
      if(ncat<=2.or.tune<=0.0_dp)return
      rrow=merge(1,regime,shared);allocate(gp(ncat+1));gp=gamma(rrow,1:ncat+1)
      do c=3,ncat
         if(c==ncat)then;gp(c)=rtruncnorm(gamma(rrow,c),tune,gp(c-1),huge(1.0_dp)/10.0_dp)
         else;gp(c)=rtruncnorm(gamma(rrow,c),tune,gp(c-1),gamma(rrow,c+1));end if
      end do
      llr=0.0_dp
      do i=1,size(y)
         if((.not.shared).and.s(i)/=regime)cycle
         eta=dot_product(x(i,:),beta(s(i),:))
         pnew=ordprob(y(i),ncat,eta,gp);pold=ordprob(y(i),ncat,eta,gamma(rrow,:));llr=llr+log(pnew)-log(pold)
      end do
      qratio=0.0_dp
      do c=3,ncat
         nf=normal_cdf((gamma(rrow,c+1)-gamma(rrow,c))/tune)-normal_cdf((gp(c-1)-gamma(rrow,c))/tune)
         nr=normal_cdf((gp(c+1)-gp(c))/tune)-normal_cdf((gamma(rrow,c-1)-gp(c))/tune)
         qratio=qratio+log(max(nf,tiny(1.0_dp)))-log(max(nr,tiny(1.0_dp)))
      end do
      if(log(runif())<min(0.0_dp,llr+qratio))then;gamma(rrow,1:ncat+1)=gp;accept_count=accept_count+1;end if
   end subroutine
end module mcmcpack_oprobit_change
