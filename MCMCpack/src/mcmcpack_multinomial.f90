! SPDX-License-Identifier: GPL-3.0-only
! Multinomial-logit Metropolis sampler translated from MCMCmnlMH.cc/MCMCmnl.h.
module mcmcpack_multinomial
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : runif, rmvnorm, rchisq
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_samplers, only : mcmc_result
   implicit none
   private
   public :: mcmc_mnl, mnl_logpost
contains
   real(dp) function mnl_logpost(y,x,beta,b0,b0prec) result(lp)
      real(dp),intent(in)::y(:,:),x(:,:),beta(:),b0(:),b0prec(:,:)
      integer::n,c,k,i,j,row
      real(dp)::eta(size(y,2)),mx,den,d(size(beta))
      n=size(y,1);c=size(y,2);k=size(beta);lp=0.0_dp
      if(size(x,1)/=n*c .or. size(x,2)/=k)then;lp=-huge(1.0_dp);return;end if
      do i=1,n
         mx=-huge(1.0_dp)
         do j=1,c
            row=(i-1)*c+j;eta(j)=dot_product(x(row,:),beta)
            if(y(i,j)/=-999.0_dp)mx=max(mx,eta(j))
         end do
         den=0.0_dp
         do j=1,c;if(y(i,j)/=-999.0_dp)den=den+exp(eta(j)-mx);end do
         if(den<=0.0_dp)then;lp=-huge(1.0_dp);return;end if
         do j=1,c;if(y(i,j)==1.0_dp)lp=lp+eta(j)-mx-log(den);end do
      end do
      d=beta-b0;lp=lp-0.5_dp*dot_product(d,matmul(b0prec,d))
   end function mnl_logpost

   real(dp) function log_multit_kernel(theta,mu,vinv,df) result(v)
      real(dp),intent(in)::theta(:),mu(:),vinv(:,:),df
      real(dp)::d(size(theta)),q
      d=theta-mu;q=dot_product(d,matmul(vinv,d))
      v=-0.5_dp*(df+real(size(theta),dp))*log(1.0_dp+q/df)
   end function log_multit_kernel

   function mcmc_mnl(y,x,beta_start,beta_mode,b0,b0prec,v,tune,burnin,mcmc,thin,rw,tdf) result(res)
      real(dp),intent(in)::y(:,:),x(:,:),beta_start(:),beta_mode(:),b0(:),b0prec(:,:),v(:,:),tune(:,:)
      integer,intent(in)::burnin,mcmc,thin
      logical,intent(in),optional::rw
      real(dp),intent(in),optional::tdf
      type(mcmc_result)::res
      integer::k,nstore,iter,count,accepts,info
      real(dp)::beta(size(beta_start)),can(size(beta_start)),z(size(beta_start)),zero(size(beta_start))
      real(dp)::propv(size(v,1),size(v,2)),pinv(size(v,1),size(v,2)),cur,lpcan,jcur,jcan,ratio,df
      logical::dorw
      k=size(beta_start);nstore=mcmc/thin;dorw=.false.;if(present(rw))dorw=rw;df=6.0_dp;if(present(tdf))df=tdf
      if(nstore<=0.or.df<=0.0_dp.or.size(beta_mode)/=k.or.size(b0)/=k.or.any(shape(b0prec)/=[k,k]).or. &
         any(shape(v)/=[k,k]).or.any(shape(tune)/=[k,k]))then;res%status=1;return;end if
      allocate(res%draws(nstore,k));propv=matmul(tune,matmul(v,transpose(tune)))
      call inv_spd(propv,pinv,info);if(info/=0)then;res%status=10+info;return;end if
      beta=beta_start;zero=0.0_dp;cur=mnl_logpost(y,x,beta,b0,b0prec);jcur=log_multit_kernel(beta,beta_mode,pinv,df)
      count=0;accepts=0
      do iter=0,burnin+mcmc-1
         if(dorw)then
            call rmvnorm(zero,propv,z,info);if(info/=0)then;res%status=20+info;return;end if;can=beta+z
            lpcan=mnl_logpost(y,x,can,b0,b0prec);ratio=lpcan-cur
            if(log(runif())<min(0.0_dp,ratio))then;beta=can;cur=lpcan;accepts=accepts+1;end if
         else if(runif()<0.75_dp)then
            call rmvnorm(zero,propv,z,info);if(info/=0)then;res%status=20+info;return;end if
            can=beta_mode+z*sqrt(df/rchisq(df));lpcan=mnl_logpost(y,x,can,b0,b0prec);jcan=log_multit_kernel(can,beta_mode,pinv,df)
            ratio=lpcan-jcan-cur+jcur
            if(log(runif())<min(0.0_dp,ratio))then;beta=can;cur=lpcan;jcur=jcan;accepts=accepts+1;end if
         else
            can=2.0_dp*beta_mode-beta;lpcan=mnl_logpost(y,x,can,b0,b0prec);ratio=lpcan-cur
            if(log(runif())<min(0.0_dp,ratio))then
               beta=can;cur=lpcan;jcur=log_multit_kernel(beta,beta_mode,pinv,df)
               accepts=accepts+1
            end if
         end if
         if(iter>=burnin.and.mod(iter,thin)==0)then;count=count+1;res%draws(count,:)=beta;end if
      end do
      res%accept_rate=real(accepts,dp)/real(burnin+mcmc,dp)
   end function mcmc_mnl
end module mcmcpack_multinomial
