! SPDX-License-Identifier: GPL-3.0-only
! Stochastic-search variable-selection quantile regression.
! Translated from cSSVSquantreg.cc and the QR_SSVS/ALaplace helpers in MCMCfcds.h.
module mcmcpack_ssvs_quantreg
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : runif, rbeta, rgamma_mt, rinvgauss, rmvnorm
   use mcmcpack_linalg, only : inv_spd, logdet_spd
   implicit none
   private
   public :: ssvs_quantreg_result, ssvs_quantreg

   type :: ssvs_quantreg_result
      integer, allocatable :: gamma(:,:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: pi0(:)
      integer :: status = 0
   end type ssvs_quantreg_result
contains
   subroutine draw_al_weights(e,w)
      real(dp), intent(in) :: e(:)
      real(dp), intent(out) :: w(size(e))
      integer :: i
      real(dp) :: nu, z
      do i=1,size(e)
         nu = 1.0_dp/max(abs(e(i)),1.0e-10_dp)
         z = rinvgauss(nu,1.0_dp)
         w(i) = 1.0_dp/max(z,1.0e-12_dp)
      end do
   end subroutine draw_al_weights

   subroutine active_system(x,y,tau,w,gamma,lambda,q,prec,rhs,cov,mu,active,info)
      real(dp),intent(in)::x(:,:),y(:),tau,w(:),lambda(:)
      integer,intent(in)::gamma(:),q
      real(dp),allocatable,intent(out)::prec(:,:),rhs(:),cov(:,:),mu(:)
      integer,allocatable,intent(out)::active(:)
      integer,intent(out)::info
      integer::n,k,na,i,j,a,b,idx
      real(dp)::u(size(y))
      n=size(x,1);k=size(x,2);na=count(gamma==1)
      allocate(active(na),prec(na,na),rhs(na),cov(na,na),mu(na))
      idx=0
      do j=1,k
         if(gamma(j)==1)then;idx=idx+1;active(idx)=j;end if
      end do
      u=y-(1.0_dp-2.0_dp*tau)*w
      prec=0.0_dp;rhs=0.0_dp
      do a=1,na
         j=active(a)
         do i=1,n
            rhs(a)=rhs(a)+0.5_dp*x(i,j)*u(i)/max(w(i),1.0e-12_dp)
         end do
         do b=1,a
            do i=1,n
               prec(a,b)=prec(a,b)+0.5_dp*x(i,j)*x(i,active(b))/max(w(i),1.0e-12_dp)
            end do
            prec(b,a)=prec(a,b)
         end do
         if(j>q)prec(a,a)=prec(a,a)+lambda(j-q)
      end do
      call inv_spd(prec,cov,info)
      if(info/=0)return
      mu=matmul(cov,rhs)
   end subroutine active_system

   real(dp) function log_model_weight(x,y,tau,w,gamma,lambda,q,pi0,ok) result(val)
      real(dp),intent(in)::x(:,:),y(:),tau,w(:),lambda(:),pi0
      integer,intent(in)::gamma(:),q
      logical,intent(out)::ok
      real(dp),allocatable::prec(:,:),rhs(:),cov(:,:),mu(:)
      integer,allocatable::active(:)
      integer::info,a,j,nopt,nin
      real(dp)::ld
      ok=.false.;val=-huge(1.0_dp)
      call active_system(x,y,tau,w,gamma,lambda,q,prec,rhs,cov,mu,active,info)
      if(info/=0)return
      call logdet_spd(prec,ld,info);if(info/=0)return
      val=0.5_dp*dot_product(rhs,mu)-0.5_dp*ld
      do a=1,size(active)
         j=active(a);if(j>q)val=val+0.5_dp*log(max(lambda(j-q),tiny(1.0_dp)))
      end do
      nopt=size(gamma)-q;nin=sum(gamma(q+1:))
      if(nopt>0)val=val+real(nin,dp)*log(max(pi0,tiny(1.0_dp)))+ &
         real(nopt-nin,dp)*log(max(1.0_dp-pi0,tiny(1.0_dp)))
      ok=.true.
   end function log_model_weight

   function ssvs_quantreg(y,x,tau,q,burnin,mcmc,thin,pi0a0,pi0b0) result(res)
      real(dp),intent(in)::y(:),x(:,:),tau,pi0a0,pi0b0
      integer,intent(in)::q,burnin,mcmc,thin
      type(ssvs_quantreg_result)::res
      integer::n,k,nstore,iter,keep,j,na,info,nin
      integer::gamma(size(x,2)),trial(size(x,2))
      real(dp)::beta(size(x,2)),w(size(y)),e(size(y)),pi0,l0,l1,pinc
      real(dp),allocatable::lambda(:),prec(:,:),rhs(:),cov(:,:),mu(:),draw(:)
      integer,allocatable::active(:)
      logical::ok0,ok1
      n=size(x,1);k=size(x,2);nstore=mcmc/thin
      if(size(y)/=n.or.q<0.or.q>k.or.tau<=0.0_dp.or.tau>=1.0_dp.or.thin<=0.or.nstore<=0.or. &
         pi0a0<=0.0_dp.or.pi0b0<=0.0_dp)then;res%status=1;return;end if
      gamma=1;beta=0.0_dp;allocate(lambda(max(k-q,0)))
      if(k>q)then
         do j=1,k-q;lambda(j)=rgamma_mt(1.0_dp,2.0_dp);end do ! Exp(rate=.5)
      end if
      pi0=rbeta(pi0a0,pi0b0)
      do j=1,n;w(j)=rgamma_mt(1.0_dp,1.0_dp/(tau*(1.0_dp-tau)));end do
      allocate(res%gamma(nstore,k),res%beta(nstore,k),res%pi0(nstore));keep=0
      do iter=0,burnin+mcmc-1
         ! Exact integrated Gaussian indicator conditionals under the AL augmentation.
         do j=q+1,k
            trial=gamma;trial(j)=0;l0=log_model_weight(x,y,tau,w,trial,lambda,q,pi0,ok0)
            trial(j)=1;l1=log_model_weight(x,y,tau,w,trial,lambda,q,pi0,ok1)
            if(.not.ok0.and..not.ok1)then;res%status=10;return;end if
            if(.not.ok0)then;pinc=1.0_dp
            else if(.not.ok1)then;pinc=0.0_dp
            else if(l1-l0>40.0_dp)then;pinc=1.0_dp
            else if(l1-l0< -40.0_dp)then;pinc=0.0_dp
            else;pinc=1.0_dp/(1.0_dp+exp(l0-l1));end if
            gamma(j)=merge(1,0,runif()<pinc)
         end do
         beta=0.0_dp;na=count(gamma==1)
         if(na>0)then
            call active_system(x,y,tau,w,gamma,lambda,q,prec,rhs,cov,mu,active,info)
            if(info/=0)then;res%status=20+info;return;end if
            allocate(draw(na));call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=30+info;return;end if
            do j=1,na;beta(active(j))=draw(j);end do
            deallocate(draw,prec,rhs,cov,mu,active)
         end if
         e=y-matmul(x,beta);call draw_al_weights(e,w)
         do j=q+1,k
            if(gamma(j)==1)then
               lambda(j-q)=rgamma_mt(1.0_dp,1.0_dp/(0.5_dp*(1.0_dp+beta(j)**2)))
            else
               lambda(j-q)=rgamma_mt(1.0_dp,2.0_dp)
            end if
         end do
         nin=sum(gamma(q+1:));pi0=rbeta(pi0a0+real(nin,dp),pi0b0+real(k-q-nin,dp))
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1;res%gamma(keep,:)=gamma;res%beta(keep,:)=beta;res%pi0(keep)=pi0
         end if
      end do
   end function ssvs_quantreg
end module mcmcpack_ssvs_quantreg
