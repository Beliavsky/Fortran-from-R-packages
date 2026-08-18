! SPDX-License-Identifier: GPL-3.0-only
! Wakefield hierarchical ecological-inference and Quinn dynamic-EI samplers.
! Translated from MCMChierEI.cc and MCMCdynamicEI.cc.  The original adaptive
! doubling slice sampler is represented by a standard stepping-out/shrinkage
! slice transition targeting the same one-dimensional full conditionals.
module mcmcpack_ei
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : runif,rnorm,rinvgamma_rng
   use mcmcpack_special_utils, only : logistic_safe,normal_logpdf
   implicit none
   private
   public :: ei_result,mcmc_hier_ei,mcmc_dynamic_ei
   type :: ei_result
      real(dp), allocatable :: draws(:,:)
      integer :: status=0
   end type ei_result
contains
   real(dp) function ei_logpost_pair(t0,t1,r0,r1,c0,mu0,mu1,sig0,sig1) result(v)
      real(dp),intent(in)::t0,t1,r0,r1,c0,mu0,mu1,sig0,sig1
      real(dp)::p0,p1,m,var
      if(sig0<=0.0_dp.or.sig1<=0.0_dp)then;v=-huge(1.0_dp);return;end if
      p0=logistic_safe(t0);p1=logistic_safe(t1)
      m=r0*p0+r1*p1
      var=r0*p0*(1.0_dp-p0)+r1*p1*(1.0_dp-p1)
      var=max(var,1.0e-12_dp)
      v=normal_logpdf(t0,mu0,sig0)+normal_logpdf(t1,mu1,sig1)+normal_logpdf(c0,m,var)
   end function ei_logpost_pair

   subroutine slice_ei_coord(theta,index,r0,r1,c0,mu0,mu1,sig0,sig1,width)
      real(dp),intent(inout)::theta(2)
      integer,intent(in)::index
      real(dp),intent(in)::r0,r1,c0,mu0,mu1,sig0,sig1,width
      real(dp)::x0,level,l,r,x,lp
      integer::j,k,tries
      x0=theta(index)
      level=ei_logpost_pair(theta(1),theta(2),r0,r1,c0,mu0,mu1,sig0,sig1)+log(runif())
      l=x0-width*runif();r=l+width
      j=20;k=20
      do while(j>0)
         if(index==1)then;lp=ei_logpost_pair(l,theta(2),r0,r1,c0,mu0,mu1,sig0,sig1)
         else;lp=ei_logpost_pair(theta(1),l,r0,r1,c0,mu0,mu1,sig0,sig1);end if
         if(lp<=level)exit;l=l-width;j=j-1
      end do
      do while(k>0)
         if(index==1)then;lp=ei_logpost_pair(r,theta(2),r0,r1,c0,mu0,mu1,sig0,sig1)
         else;lp=ei_logpost_pair(theta(1),r,r0,r1,c0,mu0,mu1,sig0,sig1);end if
         if(lp<=level)exit;r=r+width;k=k-1
      end do
      tries=0
      do
         x=l+runif()*(r-l)
         if(index==1)then;lp=ei_logpost_pair(x,theta(2),r0,r1,c0,mu0,mu1,sig0,sig1)
         else;lp=ei_logpost_pair(theta(1),x,r0,r1,c0,mu0,mu1,sig0,sig1);end if
         if(lp>=level)then;theta(index)=x;return;end if
         if(x<x0)then;l=x;else;r=x;end if
         tries=tries+1;if(tries>10000)return
      end do
   end subroutine slice_ei_coord

   function mcmc_hier_ei(r0,r1,c0,c1,mu0_prior,var0_prior,mu1_prior,var1_prior,a0,b0,a1,b1,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::r0(:),r1(:),c0(:),c1(:),mu0_prior,var0_prior,mu1_prior,var1_prior,a0,b0,a1,b1
      integer,intent(in)::burnin,mcmc,thin
      type(ei_result)::res
      integer::n,nstore,iter,keep,i
      real(dp)::theta0(size(r0)),theta1(size(r0)),mu0,mu1,sig0,sig1,postv,postm,sse
      real(dp)::pair(2),width
      n=size(r0);nstore=mcmc/thin
      if(n<1.or.size(r1)/=n.or.size(c0)/=n.or.size(c1)/=n.or.any(abs((r0+r1)-(c0+c1))>1.0e-8_dp).or. &
         var0_prior<=0.0_dp.or.var1_prior<=0.0_dp.or.a0<=0.0_dp.or.b0<=0.0_dp.or.a1<=0.0_dp.or.b1<=0.0_dp.or.nstore<=0)then
         res%status=1;return
      end if
      allocate(res%draws(nstore,2*n+4));theta0=0.0_dp;theta1=0.0_dp;mu0=0.0_dp;mu1=0.0_dp;sig0=1.0_dp;sig1=1.0_dp
      width=1.0_dp;keep=0
      do iter=0,burnin+mcmc-1
         do i=1,n
            pair=[theta0(i),theta1(i)]
            call slice_ei_coord(pair,1,r0(i),r1(i),c0(i),mu0,mu1,sig0,sig1,width)
            call slice_ei_coord(pair,2,r0(i),r1(i),c0(i),mu0,mu1,sig0,sig1,width)
            theta0(i)=pair(1);theta1(i)=pair(2)
         end do
         postv=1.0_dp/(1.0_dp/var0_prior+real(n,dp)/sig0)
         postm=postv*(sum(theta0)/sig0+mu0_prior/var0_prior)
         mu0=rnorm(postm,sqrt(postv))
         postv=1.0_dp/(1.0_dp/var1_prior+real(n,dp)/sig1)
         postm=postv*(sum(theta1)/sig1+mu1_prior/var1_prior)
         mu1=rnorm(postm,sqrt(postv))
         sse=sum((theta0-mu0)**2);sig0=rinvgamma_rng(0.5_dp*(a0+real(n,dp)),0.5_dp*(b0+sse))
         sse=sum((theta1-mu1)**2);sig1=rinvgamma_rng(0.5_dp*(a1+real(n,dp)),0.5_dp*(b1+sse))
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1
            do i=1,n;res%draws(keep,i)=logistic_safe(theta0(i));res%draws(keep,n+i)=logistic_safe(theta1(i));end do
            res%draws(keep,2*n+1:2*n+4)=[mu0,mu1,sig0,sig1]
         end if
      end do
   end function mcmc_hier_ei

   function mcmc_dynamic_ei(r0,r1,c0,c1,w,a0,b0,a1,b1,burnin,mcmc,thin) result(res)
      real(dp),intent(in)::r0(:),r1(:),c0(:),c1(:),w(:,:),a0,b0,a1,b1
      integer,intent(in)::burnin,mcmc,thin
      type(ei_result)::res
      integer::n,nstore,iter,keep,i,j
      real(dp)::theta0(size(r0)),theta1(size(r0)),sig0,sig1,wsum(size(r0)),mu0,mu1,var0,var1,sse0,sse1
      real(dp)::pair(2),width,d0,d1
      n=size(r0);nstore=mcmc/thin
      if(n<2.or.any(shape(w)/=[n,n]).or.size(r1)/=n.or.size(c0)/=n.or.size(c1)/=n.or. &
         any(abs((r0+r1)-(c0+c1))>1.0e-8_dp).or.any(w<0.0_dp).or.a0<=0.0_dp.or.b0<=0.0_dp.or. &
         a1<=0.0_dp.or.b1<=0.0_dp.or.nstore<=0)then;res%status=1;return;end if
      allocate(res%draws(nstore,2*n+2));theta0=0.0_dp;theta1=0.0_dp;sig0=0.0625_dp;sig1=0.0625_dp
      wsum=sum(w,dim=2);width=1.0_dp;keep=0
      do iter=0,burnin+mcmc-1
         do i=1,n
            if(wsum(i)>0.0_dp)then
               mu0=dot_product(w(i,:),theta0)/wsum(i);mu1=dot_product(w(i,:),theta1)/wsum(i)
               var0=sig0/wsum(i);var1=sig1/wsum(i)
            else
               mu0=sum(theta0)/real(n,dp);mu1=sum(theta1)/real(n,dp);var0=sig0;var1=sig1
            end if
            pair=[theta0(i),theta1(i)]
            call slice_ei_coord(pair,1,r0(i),r1(i),c0(i),mu0,mu1,var0,var1,width)
            call slice_ei_coord(pair,2,r0(i),r1(i),c0(i),mu0,mu1,var0,var1,width)
            theta0(i)=pair(1);theta1(i)=pair(2)
         end do
         sse0=0.0_dp;sse1=0.0_dp
         do i=1,n;do j=i+1,n
            if(w(i,j)>0.0_dp)then
               d0=theta0(i)-theta0(j);d1=theta1(i)-theta1(j)
               sse0=sse0+w(i,j)*d0*d0;sse1=sse1+w(i,j)*d1*d1
            end if
         end do;end do
         sig0=rinvgamma_rng(0.5_dp*(a0+real(n,dp)),0.5_dp*(b0+sse0))
         sig1=rinvgamma_rng(0.5_dp*(a1+real(n,dp)),0.5_dp*(b1+sse1))
         if(iter>=burnin.and.mod(iter,thin)==0)then
            keep=keep+1
            do i=1,n;res%draws(keep,i)=logistic_safe(theta0(i));res%draws(keep,n+i)=logistic_safe(theta1(i));end do
            res%draws(keep,2*n+1:2*n+2)=[sig0,sig1]
         end if
      end do
   end function mcmc_dynamic_ei
end module mcmcpack_ei
