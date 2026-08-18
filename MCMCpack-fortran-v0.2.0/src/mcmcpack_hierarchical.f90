! SPDX-License-Identifier: GPL-3.0-only
! Hierarchical Gaussian regression sampler corresponding to MCMChregress.
! The fixed and random effects are updated by conditional Gibbs steps; this has
! the same posterior target as MCMCpack's Chib-Carlin blocked implementation.
module mcmcpack_hierarchical
   use mcmcpack_kinds, only : dp, pi
   use mcmcpack_rng, only : rmvnorm, rinvgamma_rng
   use mcmcpack_distributions, only : riwish
   use mcmcpack_linalg, only : inv_spd
   implicit none
   private
   public :: hregress_result, mcmc_hregress

   type :: hregress_result
      real(dp), allocatable :: draws(:,:)
      real(dp), allocatable :: y_pred(:)
      integer :: status = 0
   end type hregress_result
contains
   function mcmc_hregress(y,x,w,group,beta_start,b_start,vb_start,v_start,mubeta,vbeta,r,rmat,s1,s2, &
                          burnin,mcmc,thin) result(res)
      real(dp), intent(in) :: y(:),x(:,:),w(:,:),beta_start(:),b_start(:,:),vb_start(:,:),v_start
      integer, intent(in) :: group(:),burnin,mcmc,thin
      real(dp), intent(in) :: mubeta(:),vbeta(:,:),r,rmat(:,:),s1,s2
      type(hregress_result) :: res
      integer :: n,p,q,g,ng,nsamp,iter,keep,info,i,j,k,ncol
      integer, allocatable :: idx(:)
      real(dp), allocatable :: beta(:),b(:,:),vb(:,:),ivbeta(:,:),ivb(:,:),prec(:,:),cov(:,:),rhs(:),mu(:),draw(:)
      real(dp), allocatable :: yg(:),xg(:,:),wg(:,:),fit(:),resid(:),ssb(:,:)
      real(dp) :: v, sse, dev

      n=size(y); p=size(x,2); q=size(w,2); nsamp=mcmc/thin
      if(size(x,1)/=n .or. size(w,1)/=n .or. size(group)/=n .or. size(beta_start)/=p .or. &
         any(shape(vbeta)/=[p,p]) .or. size(mubeta)/=p .or. size(b_start,2)/=q .or. &
         any(shape(vb_start)/=[q,q]) .or. any(shape(rmat)/=[q,q]) .or. v_start<=0.0_dp .or. &
         nsamp<=0 .or. thin<=0 .or. burnin<0) then
         res%status=1; return
      end if
      ng=size(b_start,1)
      if(any(group<1).or.any(group>ng)) then; res%status=2; return; end if
      allocate(beta(p),b(ng,q),vb(q,q),ivbeta(p,p),ivb(q,q),fit(n),resid(n),ssb(q,q))
      beta=beta_start; b=b_start; vb=vb_start; v=v_start
      call inv_spd(vbeta,ivbeta,info); if(info/=0)then;res%status=10+info;return;end if
      ncol=p+ng*q+q*q+2
      allocate(res%draws(nsamp,ncol),res%y_pred(n));res%y_pred=0.0_dp;keep=0

      do iter=1,burnin+mcmc
         ! beta | b, V
         allocate(prec(p,p),cov(p,p),rhs(p),mu(p),draw(p))
         prec=ivbeta+matmul(transpose(x),x)/v
         fit=0.0_dp
         do i=1,n; fit(i)=dot_product(w(i,:),b(group(i),:)); end do
         rhs=matmul(ivbeta,mubeta)+matmul(transpose(x),y-fit)/v
         call inv_spd(prec,cov,info); if(info/=0)then;res%status=20+info;return;end if
         mu=matmul(cov,rhs); call rmvnorm(mu,cov,draw,info); if(info/=0)then;res%status=30+info;return;end if
         beta=draw; deallocate(prec,cov,rhs,mu,draw)

         call inv_spd(vb,ivb,info); if(info/=0)then;res%status=40+info;return;end if
         do g=1,ng
            k=count(group==g); allocate(idx(k)); j=0
            do i=1,n;if(group(i)==g)then;j=j+1;idx(j)=i;end if;end do
            allocate(xg(k,p),wg(k,q),yg(k),prec(q,q),cov(q,q),rhs(q),mu(q),draw(q))
            do j=1,k;xg(j,:)=x(idx(j),:);wg(j,:)=w(idx(j),:);yg(j)=y(idx(j));end do
            prec=ivb+matmul(transpose(wg),wg)/v
            rhs=matmul(transpose(wg),yg-matmul(xg,beta))/v
            call inv_spd(prec,cov,info); if(info/=0)then;res%status=50+info;return;end if
            mu=matmul(cov,rhs);call rmvnorm(mu,cov,draw,info);if(info/=0)then;res%status=60+info;return;end if
            b(g,:)=draw
            deallocate(idx,xg,wg,yg,prec,cov,rhs,mu,draw)
         end do

         ssb=0.0_dp
         do g=1,ng
            do i=1,q;do j=1,q;ssb(i,j)=ssb(i,j)+b(g,i)*b(g,j);end do;end do
         end do
         call riwish(r+real(ng,dp),ssb+r*rmat,vb,info);if(info/=0)then;res%status=70+info;return;end if

         do i=1,n;fit(i)=dot_product(x(i,:),beta)+dot_product(w(i,:),b(group(i),:));end do
         resid=y-fit;sse=dot_product(resid,resid)
         v=rinvgamma_rng(s1+0.5_dp*real(n,dp),s2+0.5_dp*sse)
         dev=real(n,dp)*log(2.0_dp*pi*v)+sse/v

         if(iter>burnin.and.mod(iter-burnin,thin)==0)then
            keep=keep+1;k=0
            res%draws(keep,1:p)=beta;k=p
            do g=1,ng;res%draws(keep,k+1:k+q)=b(g,:);k=k+q;end do
            do i=1,q;do j=1,q;k=k+1;res%draws(keep,k)=vb(i,j);end do;end do
            res%draws(keep,k+1)=v;res%draws(keep,k+2)=dev
            res%y_pred=res%y_pred+fit/real(nsamp,dp)
         end if
      end do
   end function mcmc_hregress
end module mcmcpack_hierarchical
