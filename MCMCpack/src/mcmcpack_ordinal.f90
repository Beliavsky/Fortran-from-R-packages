! SPDX-License-Identifier: GPL-3.0-only
! Ordinal probit sampler, Cowles threshold update, translated from cMCMCoprobit.cc.
module mcmcpack_ordinal
   use mcmcpack_kinds, only : dp
   use mcmcpack_math, only : normal_cdf
   use mcmcpack_rng, only : runif,rtruncnorm,rmvnorm
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_samplers, only : mcmc_result
   implicit none
   private
   public :: mcmc_oprobit
contains
   pure real(dp) function ncdf_ms(x,mu,sd) result(p)
      real(dp),intent(in)::x,mu,sd
      p=normal_cdf((x-mu)/sd)
   end function ncdf_ms

   real(dp) function ordinal_loglik(y,mu,gamma) result(v)
      integer,intent(in)::y(:)
      real(dp),intent(in)::mu(:),gamma(:)
      integer::i,c,ncat
      real(dp)::p
      ncat=size(gamma)-1;v=0.0_dp
      do i=1,size(y)
         c=y(i)
         if(c<1.or.c>ncat)then;v=-huge(1.0_dp);return;end if
         p=normal_cdf(gamma(c+1)-mu(i))-normal_cdf(gamma(c)-mu(i))
         if(p<=tiny(1.0_dp))then;v=-huge(1.0_dp);return;end if
         v=v+log(p)
      end do
   end function ordinal_loglik

   function mcmc_oprobit(y,x,beta_start,gamma_start,b0,b0prec,tune,burnin,mcmc,thin) result(res)
      integer,intent(in)::y(:)
      real(dp),intent(in)::x(:,:),beta_start(:),gamma_start(:),b0(:),b0prec(:,:),tune
      integer,intent(in)::burnin,mcmc,thin
      type(mcmc_result)::res
      integer::n,k,ncat,nsamp,iter,count,j,i,info,accepts
      real(dp)::beta(size(beta_start)),gamma(size(gamma_start)),gp(size(gamma_start)),z(size(y)),mu(size(y))
      real(dp),allocatable::xtx(:,:),cov(:,:),rhs(:),mean(:)
      real(dp)::llcur,llnew,logqrat,lo,hi,num,den
      n=size(y);k=size(beta);ncat=maxval(y);nsamp=mcmc/thin
      if(size(x,1)/=n.or.size(x,2)/=k.or.size(gamma_start)/=ncat+1.or.tune<=0.0_dp.or.nsamp<=0)then;res%status=1;return;end if
      allocate(res%draws(nsamp,k+max(0,ncat-2)),xtx(k,k),cov(k,k),rhs(k),mean(k));xtx=matmul(transpose(x),x)
      call inv_spd(b0prec+xtx,cov,info);if(info/=0)then;res%status=10+info;return;end if
      beta=beta_start;gamma=gamma_start;gp=gamma;count=0;accepts=0
      do iter=0,burnin+mcmc-1
         mu=matmul(x,beta);gp=gamma
         ! Propose all unidentified interior cutpoints except fixed gamma(2)=0.
         do j=3,ncat
            lo=gp(j-1)
            if(j==ncat)then;hi=huge(1.0_dp)/10.0_dp;else;hi=gamma(j+1);end if
            gp(j)=rtruncnorm(gamma(j),tune,lo,hi)
         end do
         llcur=ordinal_loglik(y,mu,gamma);llnew=ordinal_loglik(y,mu,gp);logqrat=0.0_dp
         do j=3,ncat
            num=ncdf_ms(gamma(j+1),gamma(j),tune)-ncdf_ms(gp(j-1),gamma(j),tune)
            den=ncdf_ms(gp(j+1),gp(j),tune)-ncdf_ms(gamma(j-1),gp(j),tune)
            logqrat=logqrat+log(max(num,tiny(1.0_dp)))-log(max(den,tiny(1.0_dp)))
         end do
         if(log(runif())<min(0.0_dp,llnew-llcur+logqrat))then;gamma=gp;accepts=accepts+1;end if
         mu=matmul(x,beta)
         do i=1,n;z(i)=rtruncnorm(mu(i),1.0_dp,gamma(y(i)),gamma(y(i)+1));end do
         rhs=matmul(b0prec,b0)+matmul(transpose(x),z);mean=matmul(cov,rhs);call rmvnorm(mean,cov,beta,info)
         if(info/=0)then;res%status=20+info;return;end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            count=count+1;res%draws(count,1:k)=beta
            if(ncat>2)res%draws(count,k+1:)=gamma(3:ncat)
         end if
      end do
      res%accept_rate=real(accepts,dp)/real(burnin+mcmc,dp)
   end function mcmc_oprobit
end module mcmcpack_ordinal
