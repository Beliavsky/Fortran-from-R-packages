! SPDX-License-Identifier: GPL-3.0-only
! Hierarchical one-dimensional IRT sampler translated from cMCMCirtHier1d.cc.
! Includes the Liu-Wu parameter-expanded latent-data kernel and the Chib-style
! level-2 marginal-likelihood ordinate used by MCMCpack.
module mcmcpack_irthier
   use mcmcpack_kinds, only : dp,pi
   use mcmcpack_rng, only : rnorm,rtruncnorm,rmvnorm,rinvgamma_rng
   use mcmcpack_linalg, only : inv_spd,logdet_spd
   implicit none
   private
   public :: irthier_result,mcmc_irt_hier1d

   type :: irthier_result
      real(dp),allocatable :: draws(:,:)
      real(dp) :: log_marginal=-huge(1.0_dp)
      integer :: status=0
   end type irthier_result
contains
   pure real(dp) function log_ig(x,a,b) result(v)
      real(dp),intent(in)::x,a,b
      if(x<=0.0_dp.or.a<=0.0_dp.or.b<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      v=a*log(b)-log_gamma(a)-(a+1.0_dp)*log(x)-b/x
   end function log_ig

   real(dp) function log_mvn_prec(x,m,prec,info) result(v)
      real(dp),intent(in)::x(:),m(:),prec(:,:)
      integer,intent(out)::info
      real(dp)::d(size(x)),ld
      d=x-m;call logdet_spd(prec,ld,info)
      if(info/=0)then;v=-huge(1.0_dp);return;end if
      v=0.5_dp*ld-0.5_dp*real(size(x),dp)*log(2.0_dp*pi)-0.5_dp*dot_product(d,matmul(prec,d))
   end function log_mvn_prec

   function mcmc_irt_hier1d(x,theta_start,eta_start,ab0,ab0prec,xj,beta_start,b0,b0prec,c0,d0, &
                            burnin,mcmc,thin,store_theta,store_items,px,px_a0,px_b0,chib) result(res)
      integer,intent(in)::x(:,:),burnin,mcmc,thin
      real(dp),intent(in)::theta_start(:),eta_start(:,:),ab0(:),ab0prec(:,:),xj(:,:),beta_start(:),b0(:),b0prec(:,:),c0,d0
      logical,intent(in),optional::store_theta,store_items,px,chib
      real(dp),intent(in),optional::px_a0,px_b0
      type(irthier_result)::res
      integer::nj,nitem,l,nstore,iter,keep,i,k,info,ncol,col,df,totiter
      real(dp)::theta(size(theta_start)),eta(size(eta_start,1),2),beta(size(beta_start)),z(size(x,1),size(x,2))
      real(dp)::theta_hat(size(theta_start)),eta_hat(size(eta_start,1),2)
      real(dp)::sigma2,mean_z,var_theta,sd_theta,rhs_theta,sse,rss,alpha_px,alpha1,pa,pb
      real(dp)::tpt(2,2),eta_cov(2,2),eta_rhs(2),eta_mu(2),eta_draw(2),ab0rhs(2)
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mu(:),draw(:)
      real(dp)::theta_acc(size(theta_start)),beta_acc(size(beta_start)),sigma_acc
      logical::st,si,use_px,use_chib
      nj=size(x,1);nitem=size(x,2);l=size(xj,2);nstore=mcmc/thin;totiter=burnin+mcmc
      st=.true.;if(present(store_theta))st=store_theta
      si=.true.;if(present(store_items))si=store_items
      use_px=.false.;if(present(px))use_px=px
      use_chib=.false.;if(present(chib))use_chib=chib
      pa=1.0_dp;pb=1.0_dp;if(present(px_a0))pa=px_a0;if(present(px_b0))pb=px_b0
      if(size(theta_start)/=nj.or.size(eta_start,1)/=nitem.or.size(eta_start,2)/=2.or.size(ab0)/=2.or. &
         any(shape(ab0prec)/=[2,2]).or.size(xj,1)/=nj.or.size(beta_start)/=l.or.size(b0)/=l.or. &
         any(shape(b0prec)/=[l,l]).or.c0<=0.0_dp.or.d0<=0.0_dp.or.thin<=0.or.nstore<=0.or. &
         (use_px.and.(pa<=0.0_dp.or.pb<=0.0_dp)).or.(use_chib.and..not.st))then;res%status=1;return;end if
      theta=theta_start;eta=eta_start;theta_hat=theta_start;eta_hat=eta_start;beta=beta_start
      sse=sum((theta-matmul(xj,beta))**2)
      sigma2=rinvgamma_rng(0.5_dp*(c0+real(nj,dp)),0.5_dp*(d0+sse));if(sigma2<=0.0_dp)then;res%status=2;return;end if
      if(use_px)then;alpha_px=rinvgamma_rng(pa,pb);else;alpha_px=1.0_dp;end if
      ncol=merge(nj,0,st)+merge(2*nitem,0,si)+l+1
      allocate(res%draws(nstore,ncol));res%draws=0.0_dp;keep=0;ab0rhs=matmul(ab0prec,ab0)
      theta_acc=0.0_dp;beta_acc=0.0_dp;sigma_acc=0.0_dp
      do iter=0,totiter-1
         ! Albert-Chib latent utilities or the MCMCpack parameter-expanded W update.
         rss=0.0_dp;df=0
         do i=1,nj;do k=1,nitem
            if(use_px)then
               mean_z=alpha_px*(-eta(k,1)+theta(i)*eta(k,2))
               select case(x(i,k))
               case(1);z(i,k)=rtruncnorm(mean_z,alpha_px,0.0_dp,huge(1.0_dp)/10.0_dp);df=df+1
               case(0);z(i,k)=rtruncnorm(mean_z,alpha_px,-huge(1.0_dp)/10.0_dp,0.0_dp);df=df+1
               case default;z(i,k)=rnorm(mean_z,alpha_px)
               end select
               z(i,k)=z(i,k)/alpha_px
               rss=rss+(z(i,k)-(-eta_hat(k,1)+theta_hat(i)*eta_hat(k,2)))**2
            else
               mean_z=-eta(k,1)+theta(i)*eta(k,2)
               select case(x(i,k))
               case(1);z(i,k)=rtruncnorm(mean_z,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
               case(0);z(i,k)=rtruncnorm(mean_z,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
               case default;z(i,k)=rnorm(mean_z,1.0_dp)
               end select
            end if
         end do;end do
         if(use_px)then
            alpha1=rinvgamma_rng(0.5_dp*(pa+real(df,dp)),0.5_dp*(pb+rss))
            alpha_px=sqrt(max(alpha1/max(alpha_px,1.0e-12_dp),1.0e-12_dp))
         end if

         ! Subject abilities conditional on the level-2 regression.
         var_theta=1.0_dp/(sum(eta(:,2)**2)+1.0_dp/sigma2);sd_theta=sqrt(var_theta)
         do i=1,nj
            rhs_theta=sum(eta(:,2)*(z(i,:)+eta(:,1)))+dot_product(xj(i,:),beta)/sigma2
            theta_hat(i)=var_theta*rhs_theta
            theta(i)=rnorm(theta_hat(i)/alpha_px,sd_theta)
         end do

         ! Item intercept/discrimination pairs.
         tpt=0.0_dp;tpt(1,1)=real(nj,dp);tpt(1,2)=-sum(theta);tpt(2,1)=tpt(1,2);tpt(2,2)=sum(theta**2)
         call inv_spd(tpt+ab0prec,eta_cov,info);if(info/=0)then;res%status=10+info;return;end if
         do k=1,nitem
            eta_rhs(1)=-sum(z(:,k));eta_rhs(2)=dot_product(z(:,k),theta)
            eta_hat(k,:)=matmul(eta_cov,eta_rhs+ab0rhs);eta_mu=eta_hat(k,:)/alpha_px
            call rmvnorm(eta_mu,eta_cov,eta_draw,info);if(info/=0)then;res%status=20+info;return;end if;eta(k,:)=eta_draw
         end do

         ! Level-2 Gaussian regression beta and variance.
         allocate(prec(l,l),cov(l,l),rhs(l),mu(l),draw(l))
         prec=b0prec+matmul(transpose(xj),xj)/sigma2;rhs=matmul(b0prec,b0)+matmul(transpose(xj),theta)/sigma2
         call inv_spd(prec,cov,info);if(info/=0)then;res%status=30+info;return;end if
         mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=40+info;return;end if;beta=draw
         deallocate(prec,cov,rhs,mu,draw)
         sse=sum((theta-matmul(xj,beta))**2);sigma2=rinvgamma_rng(0.5_dp*(c0+real(nj,dp)),0.5_dp*(d0+sse))

         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            if(st)then;res%draws(keep,col+1:col+nj)=theta;col=col+nj;end if
            if(si)then;do k=1,nitem;res%draws(keep,col+1:col+2)=eta(k,:);col=col+2;end do;end if
            res%draws(keep,col+1:col+l)=beta;col=col+l;res%draws(keep,col+1)=sigma2
            theta_acc=theta_acc+theta;beta_acc=beta_acc+beta;sigma_acc=sigma_acc+sigma2
         end if
      end do
      if(use_chib)call chib_level2(theta_acc/real(nstore,dp),beta_acc/real(nstore,dp),sigma_acc/real(nstore,dp), &
                                    xj,b0,b0prec,c0,d0,totiter,res%log_marginal,res%status)
   end function mcmc_irt_hier1d

   subroutine chib_level2(theta_star,beta_star,sigma_star,xj,b0,b0prec,c0,d0,nscan,logm,status)
      real(dp),intent(in)::theta_star(:),beta_star(:),sigma_star,xj(:,:),b0(:),b0prec(:,:),c0,d0
      integer,intent(in)::nscan
      real(dp),intent(out)::logm
      integer,intent(inout)::status
      integer::l,n,iter,info
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mu(:),beta(:),draw(:)
      real(dp)::sig,sse,shape,scale,log_sig_ord,log_beta_ord,loglike,logprior
      l=size(beta_star);n=size(theta_star);allocate(prec(l,l),cov(l,l),rhs(l),mu(l),beta(l),draw(l));beta=beta_star
      log_sig_ord=0.0_dp
      do iter=1,max(nscan,1)
         sse=sum((theta_star-matmul(xj,beta))**2);shape=0.5_dp*(c0+real(n,dp));scale=0.5_dp*(d0+sse)
         log_sig_ord=log_sig_ord+log_ig(sigma_star,shape,scale);sig=rinvgamma_rng(shape,scale)
         prec=b0prec+matmul(transpose(xj),xj)/sig;rhs=matmul(b0prec,b0)+matmul(transpose(xj),theta_star)/sig
         call inv_spd(prec,cov,info);if(info/=0)then;status=80+info;return;end if
         mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;status=90+info;return;end if;beta=draw
      end do
      log_sig_ord=log_sig_ord/real(max(nscan,1),dp)
      prec=b0prec+matmul(transpose(xj),xj)/sigma_star;rhs=matmul(b0prec,b0)+matmul(transpose(xj),theta_star)/sigma_star
      call inv_spd(prec,cov,info);if(info/=0)then;status=100+info;return;end if;mu=matmul(cov,rhs)
      log_beta_ord=log_mvn_prec(beta_star,mu,prec,info);if(info/=0)then;status=110+info;return;end if
      loglike=-0.5_dp*real(n,dp)*log(2.0_dp*pi*sigma_star)-0.5_dp*sum((theta_star-matmul(xj,beta_star))**2)/sigma_star
      logprior=log_mvn_prec(beta_star,b0,b0prec,info)+log_ig(sigma_star,0.5_dp*c0,0.5_dp*d0)
      if(info/=0)then;status=120+info;return;end if
      logm=loglike+logprior-log_beta_ord-log_sig_ord
   end subroutine chib_level2
end module mcmcpack_irthier
