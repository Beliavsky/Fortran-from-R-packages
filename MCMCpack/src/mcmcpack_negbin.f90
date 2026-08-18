! SPDX-License-Identifier: GPL-3.0-only
! Negative-binomial regression target from cMCMCnegbin.cc/MCMCnbutil.h.
! The original auxiliary-mixture Gibbs beta update is replaced here by an
! exact-target Gaussian random-walk Metropolis step; rho retains the original
! stepping-out slice update and the same MCMCpack prior parameterization.
module mcmcpack_negbin
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : runif, rmvnorm
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_samplers, only : mcmc_result
   implicit none
   private
   public :: mcmc_negbin, negbin_logpost, rho_nb_logcond
contains
   pure real(dp) function log_rho_plus_expeta(rho,eta) result(v)
      real(dp),intent(in)::rho,eta
      if(eta>log(rho))then
         v=eta+log(1.0_dp+rho*exp(-eta))
      else
         v=log(rho)+log(1.0_dp+exp(eta)/rho)
      end if
   end function

   real(dp) function negbin_loglike(y,x,beta,rho) result(v)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:),beta(:),rho
      integer::i
      real(dp)::eta
      if(rho<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      v=0.0_dp
      do i=1,size(y)
         eta=dot_product(x(i,:),beta)
         v=v+log_gamma(rho+real(y(i),dp))-log_gamma(rho)-log_gamma(real(y(i)+1,dp)) &
             +rho*log(rho)+real(y(i),dp)*eta-(rho+real(y(i),dp))*log_rho_plus_expeta(rho,eta)
      end do
   end function negbin_loglike

   real(dp) function negbin_logpost(y,x,beta,rho,b0,b0prec,e,f,g) result(v)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:),beta(:),rho,b0(:),b0prec(:,:),e,f,g
      real(dp)::d(size(beta))
      if(rho<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      d=beta-b0
      v=negbin_loglike(y,x,beta,rho)-0.5_dp*dot_product(d,matmul(b0prec,d)) &
        +(e-1.0_dp)*log(rho)-(e+f)*log(rho+g)
   end function negbin_logpost

   real(dp) function rho_nb_logcond(rho,y,lambda,e,f,g) result(v)
      real(dp),intent(in)::rho,lambda(:),e,f,g
      integer,intent(in)::y(:)
      integer::i
      if(rho<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      v=(e-1.0_dp)*log(rho)-(e+f)*log(rho+g)
      do i=1,size(y)
         v=v+log_gamma(rho+real(y(i),dp))-log_gamma(rho)-log_gamma(real(y(i)+1,dp)) &
            +rho*log(rho)+real(y(i),dp)*log(max(lambda(i),tiny(1.0_dp))) &
            -(rho+real(y(i),dp))*log(rho+lambda(i))
      end do
   end function rho_nb_logcond

   real(dp) function slice_rho(rho,y,lambda,step,e,f,g) result(newrho)
      real(dp),intent(in)::rho,lambda(:),step,e,f,g
      integer,intent(in)::y(:)
      real(dp)::level,l,r,u,val
      integer::j,k,tries
      level=rho_nb_logcond(rho,y,lambda,e,f,g)+log(max(runif(),tiny(1.0_dp)))
      u=runif();l=max(0.0_dp,rho-step*u);r=l+step
      j=int(100.0_dp*runif());k=99-j
      do while(j>0)
         val=rho_nb_logcond(max(l,tiny(1.0_dp)),y,lambda,e,f,g);if(val<=level)exit
         l=max(0.0_dp,l-step);j=j-1
      end do
      do while(k>0)
         val=rho_nb_logcond(r,y,lambda,e,f,g);if(val<=level)exit
         r=r+step;k=k-1
      end do
      tries=0
      do
         newrho=l+runif()*(r-l);val=rho_nb_logcond(max(newrho,tiny(1.0_dp)),y,lambda,e,f,g)
         if(val>level.and.newrho>0.0_dp)exit
         if(newrho>rho)then;r=newrho;else;l=newrho;end if
         tries=tries+1;if(tries>100000)then;newrho=rho;exit;end if
      end do
   end function slice_rho

   function mcmc_negbin(y,x,beta_start,rho_start,b0,b0prec,e,f,g,rho_step,beta_tune,burnin,mcmc,thin) result(res)
      integer,intent(in)::y(:),burnin,mcmc,thin
      real(dp),intent(in)::x(:,:),beta_start(:),rho_start,b0(:),b0prec(:,:),e,f,g,rho_step
      real(dp),intent(in),optional::beta_tune
      type(mcmc_result)::res
      integer::n,k,nstore,iter,count,accepts,info
      real(dp)::beta(size(beta_start)),can(size(beta_start)),zero(size(beta_start)),rho,lp,lpc,tune
      real(dp),allocatable::basecov(:,:),propcov(:,:),z(:),lambda(:),xtx(:,:)
      n=size(y);k=size(beta_start);nstore=mcmc/thin;tune=0.15_dp;if(present(beta_tune))tune=beta_tune
      if(size(x,1)/=n.or.size(x,2)/=k.or.size(b0)/=k.or.any(shape(b0prec)/=[k,k]).or. &
         any(y<0).or.rho_start<=0.0_dp.or.rho_step<=0.0_dp.or.tune<=0.0_dp.or.nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,k+1),basecov(k,k),propcov(k,k),z(k),lambda(n),xtx(k,k));xtx=matmul(transpose(x),x)
      call inv_spd(b0prec+xtx,basecov,info);if(info/=0)then;res%status=10+info;return;end if
      propcov=tune*tune*basecov;zero=0.0_dp;beta=beta_start;rho=rho_start;lp=negbin_logpost(y,x,beta,rho,b0,b0prec,e,f,g)
      count=0;accepts=0
      do iter=0,burnin+mcmc-1
         lambda=exp(min(700.0_dp,matmul(x,beta)));rho=slice_rho(rho,y,lambda,rho_step,e,f,g)
         lp=negbin_logpost(y,x,beta,rho,b0,b0prec,e,f,g)
         call rmvnorm(zero,propcov,z,info);if(info/=0)then;res%status=20+info;return;end if;can=beta+z
         lpc=negbin_logpost(y,x,can,rho,b0,b0prec,e,f,g)
         if(log(max(runif(),tiny(1.0_dp)))<min(0.0_dp,lpc-lp))then;beta=can;lp=lpc;accepts=accepts+1;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,1:k)=beta;res%draws(count,k+1)=rho;end if
      end do
      res%accept_rate=real(accepts,dp)/real(burnin+mcmc,dp)
   end function mcmc_negbin
end module mcmcpack_negbin
