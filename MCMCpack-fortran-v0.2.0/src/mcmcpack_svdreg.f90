! SPDX-License-Identifier: GPL-3.0-only
! MCMCSVDreg translated from cMCMCSVDreg.cc (West large-p/small-n model).
module mcmcpack_svdreg
   use mcmcpack_kinds, only : dp, pi
   use mcmcpack_rng, only : runif, rnorm, rinvgamma_rng
   use mcmcpack_linalg, only : inv_spd, logdet_spd
   use mcmcpack_samplers, only : mcmc_result
   implicit none
   private
   public :: mcmc_svdreg
contains
   real(dp) function logdnorm(x,mu,sd) result(v)
      real(dp),intent(in)::x,mu,sd
      v=-0.5_dp*(log(2.0_dp*pi)+2.0_dp*log(sd)+((x-mu)/sd)**2)
   end function

   real(dp) function logdmvnorm(x,mu,cov,info) result(v)
      real(dp),intent(in)::x(:),mu(:),cov(:,:)
      integer,intent(out)::info
      real(dp),allocatable::cinv(:,:)
      real(dp)::d(size(x)),ld
      allocate(cinv(size(x),size(x)));call inv_spd(cov,cinv,info);if(info/=0)then;v=-huge(1.0_dp);return;end if
      call logdet_spd(cov,ld,info);if(info/=0)then;v=-huge(1.0_dp);return;end if
      d=x-mu;v=-0.5_dp*(real(size(x),dp)*log(2.0_dp*pi)+ld+dot_product(d,matmul(cinv,d)))
   end function

   function mcmc_svdreg(y_start,ymiss,a,dmat,fmat,tau_start,g0,a0,b0,c0,d0,w0,burnin,mcmc,thin,betasamp) result(res)
      real(dp),intent(in)::y_start(:),a(:,:),dmat(:,:),fmat(:,:),tau_start(:),g0(:),a0,b0,c0(:),d0(:),w0(:)
      logical,intent(in)::ymiss(:)
      integer,intent(in)::burnin,mcmc,thin
      logical,intent(in),optional::betasamp
      type(mcmc_result)::res
      integer::n,k,nmiss,nstore,iter,keep,i,j,col,info
      logical::bs
      real(dp)::y(size(y_start)),tau2(size(y_start)),gamma(size(y_start)),fy(size(y_start)),dg0(size(y_start))
      real(dp)::gammahat(size(y_start)),eta(size(y_start)),beta(size(a,1)),d2(size(y_start))
      real(dp)::ftd(size(y_start),size(y_start)),dinvf(size(y_start),size(y_start))
      real(dp)::sigma2,q,mstar,vstar,logw0,logw1,logf0,logf1,ldenom,wstar,gamg2
      real(dp)::resid(size(y_start)),margmean(size(y_start)),gammanot(size(y_start))
      real(dp),allocatable::margvar(:,:)
      n=size(y_start);k=size(a,1);nmiss=count(ymiss);nstore=mcmc/thin;bs=.false.;if(present(betasamp))bs=betasamp
      if(any(shape(a)/=[k,n]).or.any(shape(dmat)/=[n,n]).or.any(shape(fmat)/=[n,n]).or.size(tau_start)/=n.or. &
         size(g0)/=n.or.size(c0)/=n.or.size(d0)/=n.or.size(w0)/=n.or.size(ymiss)/=n.or.any(tau_start<=0.0_dp).or. &
         any(w0<=0.0_dp).or.any(w0>=1.0_dp).or.nstore<=0)then;res%status=1;return;end if
      do i=1,n;if(abs(dmat(i,i))<=tiny(1.0_dp))then;res%status=2;return;end if;end do
      ftd=matmul(transpose(fmat),dmat);dinvf=0.0_dp
      do i=1,n;dinvf(i,:)=fmat(i,:)/dmat(i,i);d2(i)=dmat(i,i)**2;end do
      dg0=matmul(dmat,g0);y=y_start;tau2=tau_start;gamma=0.0_dp;sigma2=1.0_dp
      allocate(res%draws(nstore,nmiss+n+n+1+merge(k,0,bs)),margvar(n,n));keep=0
      do iter=0,burnin+mcmc-1
         fy=matmul(fmat,y)-dg0;q=sum(fy*fy/(1.0_dp+tau2))
         sigma2=rinvgamma_rng(0.5_dp*(a0+real(n,dp)),0.5_dp*(b0+q))
         gammahat=matmul(dinvf,y)
         do i=1,n
            mstar=(g0(i)+tau2(i)*gammahat(i))/(1.0_dp+tau2(i))
            vstar=sigma2*tau2(i)/(d2(i)*(1.0_dp+tau2(i)))
            gammanot=gamma;gammanot(i)=0.0_dp;resid=y-matmul(ftd,gammanot);margmean=ftd(:,i)*g0(i)
            margvar=0.0_dp;do j=1,n;margvar(j,j)=sigma2;end do
            do j=1,n;margvar(j,:)=margvar(j,:)+ftd(j,i)*ftd(:,i)*(sigma2*tau2(i)/d2(i));end do
            logw0=log(w0(i));logw1=log(1.0_dp-w0(i));logf0=0.0_dp
            do j=1,n;logf0=logf0+logdnorm(resid(j),0.0_dp,sqrt(sigma2));end do
            logf1=logdmvnorm(resid,margmean,margvar,info);if(info/=0)then;res%status=10+info;return;end if
            if(logw0+logf0>logw1+logf1)then
               ldenom=logw0+logf0+log(1.0_dp+exp(logw1+logf1-logw0-logf0))
            else
               ldenom=logw1+logf1+log(1.0_dp+exp(logw0+logf0-logw1-logf1))
            end if
            wstar=exp(logw0+logf0-ldenom)
            if(runif()<wstar)then;gamma(i)=0.0_dp;else;gamma(i)=rnorm(mstar,sqrt(vstar));end if
         end do
         do i=1,n
            gamg2=(gamma(i)-g0(i))**2;tau2(i)=rinvgamma_rng(0.5_dp*(1.0_dp+c0(i)),0.5_dp*(gamg2*d2(i)/sigma2+d0(i)))
         end do
         eta=matmul(ftd,gamma);do i=1,n;if(ymiss(i))y(i)=rnorm(eta(i),sqrt(sigma2));end do
         if(bs)beta=matmul(a,gamma)
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;col=0
            do i=1,n;if(ymiss(i))then;col=col+1;res%draws(keep,col)=y(i);end if;end do
            res%draws(keep,col+1:col+n)=gamma;col=col+n;res%draws(keep,col+1:col+n)=tau2;col=col+n
            res%draws(keep,col+1)=sigma2;col=col+1;if(bs)res%draws(keep,col+1:col+k)=beta
         end if
      end do
   end function mcmc_svdreg
end module mcmcpack_svdreg
