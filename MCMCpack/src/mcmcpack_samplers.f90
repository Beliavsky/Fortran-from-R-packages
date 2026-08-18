! SPDX-License-Identifier: GPL-3.0-only
! Major model samplers translated from MCMCpack compiled C++ routines.
module mcmcpack_samplers
   use mcmcpack_kinds, only : dp
   use mcmcpack_math, only : logistic,log1pexp
   use mcmcpack_rng, only : runif,rnorm,rmvnorm,rinvgamma_rng,rtruncnorm,rinvgauss,rchisq
   use mcmcpack_linalg, only : inv_spd,chol_lower,solve_spd
   implicit none
   private
   type, public :: mcmc_result
      real(dp), allocatable :: draws(:,:)
      real(dp) :: accept_rate = 1.0_dp
      integer :: status = 0
   end type mcmc_result
   abstract interface
      function logpost_fn(theta) result(value)
         import dp
         real(dp), intent(in) :: theta(:)
         real(dp) :: value
      end function logpost_fn
   end interface
   public :: mcmc_regress,mcmc_probit,mcmc_logit,mcmc_poisson,mcmc_tobit,mcmc_quantreg
   public :: mcmc_metrop1r, logpost_fn
contains
   subroutine beta_normal_draw(xtx,xty,b0,b0prec,sigma2,beta,info)
      real(dp),intent(in)::xtx(:,:),xty(:),b0(:),b0prec(:,:),sigma2
      real(dp),intent(out)::beta(size(b0))
      integer,intent(out)::info
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mean(:)
      integer::k
      k=size(b0); allocate(prec(k,k),cov(k,k),rhs(k),mean(k))
      prec=b0prec+xtx/sigma2; rhs=matmul(b0prec,b0)+xty/sigma2
      call inv_spd(prec,cov,info); if(info/=0)return
      mean=matmul(cov,rhs); call rmvnorm(mean,cov,beta,info)
   end subroutine beta_normal_draw

   real(dp) function sigma2_draw(x,y,beta,c0,d0) result(s2)
      real(dp),intent(in)::x(:,:),y(:),beta(:),c0,d0
      real(dp)::e(size(y)),sse
      e=y-matmul(x,beta); sse=dot_product(e,e)
      s2=rinvgamma_rng(0.5_dp*(c0+real(size(y),dp)),0.5_dp*(d0+sse))
   end function sigma2_draw

   function mcmc_regress(y,x,beta_start,b0,b0prec,c0,d0,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:),c0,d0
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      real(dp)::beta(size(beta_start)),s2
      real(dp),allocatable::xtx(:,:),xty(:)
      integer::iter,count,nstore,info,k
      k=size(beta); nstore=mcmc/thin
      if(size(x,1)/=size(y).or.size(x,2)/=k.or.nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k+1),xtx(k,k),xty(k)); xtx=matmul(transpose(x),x); xty=matmul(transpose(x),y)
      beta=beta_start; count=0
      do iter=0,burnin+mcmc-1
         s2=sigma2_draw(x,y,beta,c0,d0)
         call beta_normal_draw(xtx,xty,b0,b0prec,s2,beta,info); if(info/=0)then;res%status=10+info;return;end if
         if(iter>=burnin .and. mod(iter,thin)==0)then
            count=count+1; res%draws(count,1:k)=beta; res%draws(count,k+1)=s2
         end if
      end do
   end function mcmc_regress

   function mcmc_probit(y,x,beta_start,b0,b0prec,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:)
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      real(dp)::beta(size(beta_start)),z(size(y)),mu(size(y))
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mean(:),xtx(:,:)
      integer::iter,i,k,n,nstore,count,info
      n=size(y);k=size(beta);nstore=mcmc/thin
      if(size(x,1)/=n.or.size(x,2)/=k.or.nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k),prec(k,k),cov(k,k),rhs(k),mean(k),xtx(k,k)); xtx=matmul(transpose(x),x)
      prec=b0prec+xtx; call inv_spd(prec,cov,info); if(info/=0)then;res%status=10+info;return;end if
      beta=beta_start; count=0
      do iter=0,burnin+mcmc-1
         mu=matmul(x,beta)
         do i=1,n
            if(y(i)>=0.5_dp)then
               z(i)=rtruncnorm(mu(i),1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else
               z(i)=rtruncnorm(mu(i),1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
            end if
         end do
         rhs=matmul(b0prec,b0)+matmul(transpose(x),z); mean=matmul(cov,rhs)
         call rmvnorm(mean,cov,beta,info); if(info/=0)then;res%status=20+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,:)=beta;end if
      end do
   end function mcmc_probit

   real(dp) function logpost_logit(y,x,beta,b0,b0prec) result(v)
      real(dp),intent(in)::y(:),x(:,:),beta(:),b0(:),b0prec(:,:)
      real(dp)::eta(size(y)),d(size(beta))
      eta=matmul(x,beta); d=beta-b0
      v=sum(y*eta-log1pexp(eta))-0.5_dp*dot_product(d,matmul(b0prec,d))
   end function logpost_logit

   real(dp) function logpost_poisson(y,x,beta,b0,b0prec) result(v)
      real(dp),intent(in)::y(:),x(:,:),beta(:),b0(:),b0prec(:,:)
      real(dp)::eta(size(y)),d(size(beta))
      eta=matmul(x,beta)
      if(maxval(eta)>700.0_dp)then;v=-huge(1.0_dp);return;end if
      d=beta-b0; v=sum(y*eta-exp(eta))-0.5_dp*dot_product(d,matmul(b0prec,d))
   end function logpost_poisson

   function mh_glm(y,x,beta_start,b0,b0prec,v,tune,burnin,mcmc,thin,is_logit) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:),v(:,:),tune(:,:)
      integer,intent(in)::burnin,mcmc,thin
      logical,intent(in)::is_logit
      type(mcmc_result)::res
      real(dp)::beta(size(beta_start)),can(size(beta_start)),z(size(beta_start))
      real(dp),allocatable::vinv(:,:),mid(:,:),propv(:,:),l(:,:)
      real(dp)::cur,lp,loga
      integer::k,nstore,iter,count,accepts,info,i
      k=size(beta);nstore=mcmc/thin
      if(nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k),vinv(k,k),mid(k,k),propv(k,k),l(k,k))
      call inv_spd(v,vinv,info);if(info/=0)then;res%status=11;return;end if
      call inv_spd(b0prec+vinv,mid,info);if(info/=0)then;res%status=12;return;end if
      propv=matmul(tune,matmul(mid,transpose(tune)))
      call chol_lower(propv,l,info);if(info/=0)then;res%status=13;return;end if
      beta=beta_start
      if(is_logit)then;cur=logpost_logit(y,x,beta,b0,b0prec);else;cur=logpost_poisson(y,x,beta,b0,b0prec);end if
      count=0;accepts=0
      do iter=0,burnin+mcmc-1
         do i=1,k;z(i)=rnorm();end do;can=beta+matmul(l,z)
         if(is_logit)then;lp=logpost_logit(y,x,can,b0,b0prec);else;lp=logpost_poisson(y,x,can,b0,b0prec);end if
         loga=lp-cur
         if(log(runif())<min(0.0_dp,loga))then;beta=can;cur=lp;accepts=accepts+1;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,:)=beta;end if
      end do
      res%accept_rate=real(accepts,dp)/real(burnin+mcmc,dp)
   end function mh_glm

   function mcmc_logit(y,x,beta_start,b0,b0prec,v,tune,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:),v(:,:),tune(:,:)
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      res=mh_glm(y,x,beta_start,b0,b0prec,v,tune,burnin,mcmc,thin,.true.)
   end function mcmc_logit

   function mcmc_poisson(y,x,beta_start,b0,b0prec,v,tune,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:),v(:,:),tune(:,:)
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      res=mh_glm(y,x,beta_start,b0,b0prec,v,tune,burnin,mcmc,thin,.false.)
   end function mcmc_poisson

   function mcmc_tobit(y,x,beta_start,b0,b0prec,c0,d0,below,above,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:),c0,d0,below,above
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      real(dp)::beta(size(beta_start)),z(size(y)),mu(size(y)),s2
      real(dp),allocatable::xtx(:,:),xtz(:)
      integer::n,k,i,iter,nstore,count,info
      n=size(y);k=size(beta);nstore=mcmc/thin
      if(nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k+1),xtx(k,k),xtz(k));xtx=matmul(transpose(x),x);beta=beta_start;z=y;count=0
      do iter=0,burnin+mcmc-1
         s2=sigma2_draw(x,z,beta,c0,d0);mu=matmul(x,beta)
         do i=1,n
            if(y(i)<=below)z(i)=rtruncnorm(mu(i),sqrt(s2),-huge(1.0_dp)/10.0_dp,below)
            if(y(i)>=above)z(i)=rtruncnorm(mu(i),sqrt(s2),above,huge(1.0_dp)/10.0_dp)
         end do
         xtz=matmul(transpose(x),z);call beta_normal_draw(xtx,xtz,b0,b0prec,s2,beta,info)
         if(info/=0)then;res%status=10+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,1:k)=beta;res%draws(count,k+1)=s2;end if
      end do
   end function mcmc_tobit

   subroutine quant_weights(e,w)
      real(dp),intent(in)::e(:)
      real(dp),intent(out)::w(size(e))
      real(dp)::nu,chi,small
      integer::i
      do i=1,size(e)
         nu=1.0_dp/max(abs(e(i)),sqrt(epsilon(1.0_dp)));chi=rchisq(1.0_dp)
         small=nu*(nu*chi+1.0_dp-sqrt(max(0.0_dp,nu*nu*chi*chi+2.0_dp*nu*chi)))
         small=max(small,tiny(1.0_dp))
         if(runif()<nu/(nu+small))then;w(i)=1.0_dp/small;else;w(i)=small/(nu*nu);end if
      end do
   end subroutine quant_weights

   subroutine quant_beta_draw(tau,x,y,w,b0,b0prec,beta,info)
      real(dp),intent(in)::tau,x(:,:),y(:),w(:),b0(:),b0prec(:,:)
      real(dp),intent(out)::beta(size(b0))
      integer,intent(out)::info
      integer::n,k,i,j,m
      real(dp)::u(size(y))
      real(dp),allocatable::xtwx(:,:),xtwu(:),prec(:,:),cov(:,:),mean(:),rhs(:)
      n=size(y);k=size(b0);allocate(xtwx(k,k),xtwu(k),prec(k,k),cov(k,k),mean(k),rhs(k))
      u=y-(1.0_dp-2.0_dp*tau)*w;xtwx=0.0_dp;xtwu=0.0_dp
      do i=1,k
         do m=1,n;xtwu(i)=xtwu(i)+x(m,i)*u(m)/w(m);end do
         do j=1,k
            do m=1,n;xtwx(i,j)=xtwx(i,j)+x(m,i)*x(m,j)/w(m);end do
         end do
      end do
      prec=b0prec+0.5_dp*xtwx;rhs=matmul(b0prec,b0)+0.5_dp*xtwu
      call inv_spd(prec,cov,info);if(info/=0)return;mean=matmul(cov,rhs);call rmvnorm(mean,cov,beta,info)
   end subroutine quant_beta_draw

   function mcmc_quantreg(tau,y,x,beta_start,b0,b0prec,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::tau,y(:),x(:,:),beta_start(:),b0(:),b0prec(:,:)
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      real(dp)::beta(size(beta_start)),e(size(y)),w(size(y))
      integer::nstore,k,iter,count,info
      k=size(beta);nstore=mcmc/thin
      if(tau<=0.0_dp.or.tau>=1.0_dp.or.nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k));beta=beta_start;count=0
      do iter=0,burnin+mcmc-1
         e=y-matmul(x,beta);call quant_weights(e,w);call quant_beta_draw(tau,x,y,w,b0,b0prec,beta,info)
         if(info/=0)then;res%status=10+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,:)=beta;end if
      end do
   end function mcmc_quantreg

   function mcmc_metrop1r(logpost,theta_start,proposal_cov,burnin,mcmc,thin) result(res)
      procedure(logpost_fn)::logpost
      real(dp),intent(in)::theta_start(:),proposal_cov(:,:)
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      real(dp)::theta(size(theta_start)),can(size(theta_start)),z(size(theta_start)),cur,lp
      real(dp),allocatable::l(:,:)
      integer::k,nstore,info,iter,count,accepts,i
      k=size(theta);nstore=mcmc/thin;if(nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k),l(k,k));call chol_lower(proposal_cov,l,info);if(info/=0)then;res%status=10+info;return;end if
      theta=theta_start;cur=logpost(theta);count=0;accepts=0
      do iter=0,burnin+mcmc-1
         do i=1,k;z(i)=rnorm();end do;can=theta+matmul(l,z);lp=logpost(can)
         if(log(runif())<min(0.0_dp,lp-cur))then;theta=can;cur=lp;accepts=accepts+1;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,:)=theta;end if
      end do
      res%accept_rate=real(accepts,dp)/real(burnin+mcmc,dp)
   end function mcmc_metrop1r
end module mcmcpack_samplers
