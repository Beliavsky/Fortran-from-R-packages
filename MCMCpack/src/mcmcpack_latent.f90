! SPDX-License-Identifier: GPL-3.0-only
! Latent Gaussian samplers translated from MCMCpack's MCMCfcds.h and C++ drivers.
module mcmcpack_latent
   use mcmcpack_kinds, only : dp
   use mcmcpack_rng, only : rnorm,rtruncnorm,rmvnorm
   use mcmcpack_linalg, only : inv_spd
   use mcmcpack_samplers, only : mcmc_result
   implicit none
   private
   public :: mcmc_irt1d,mcmc_paircompare
contains
   function mcmc_irt1d(response,theta_start,alpha_start,beta_start,t0,t0prec,ab0,ab0prec, &
                       theta_eq,theta_ineq,burnin,mcmc,thin,store_item,store_ability) result(res)
      integer,intent(in)::response(:,:)
      real(dp),intent(in)::theta_start(:),alpha_start(:),beta_start(:),t0,t0prec,ab0(:),ab0prec(:,:)
      real(dp),intent(in)::theta_eq(:),theta_ineq(:)
      integer,intent(in)::burnin,mcmc,thin
      logical,intent(in),optional::store_item,store_ability
      type(mcmc_result)::res
      integer::j,k,iter,count,nsamp,nj,nk,info,ncol,col
      logical::si,sa
      real(dp)::theta(size(theta_start)),eta(size(alpha_start),2),z(size(response,1),size(response,2))
      real(dp)::tpt(2,2),postv(2,2),rhs(2),postm(2),draw(2),v,m,bsum
      si=.false.;sa=.true.;if(present(store_item))si=store_item;if(present(store_ability))sa=store_ability
      nj=size(response,1);nk=size(response,2);nsamp=mcmc/thin
      if(size(theta_start)/=nj.or.size(alpha_start)/=nk.or.size(beta_start)/=nk.or.size(ab0)/=2.or. &
         any(shape(ab0prec)/=[2,2]).or.size(theta_eq)/=nj.or.size(theta_ineq)/=nj.or.nsamp<=0.or.(.not.si.and..not.sa))then
         res%status=1;return
      end if
      ncol=merge(nj,0,sa)+merge(2*nk,0,si)
      allocate(res%draws(nsamp,ncol));theta=theta_start
      eta(:,1)=alpha_start;eta(:,2)=beta_start;count=0
      do iter=0,burnin+mcmc-1
         ! Z | theta, eta, response. Values other than 0/1 are treated as missing.
         do j=1,nj;do k=1,nk
            m=-eta(k,1)+theta(j)*eta(k,2)
            if(response(j,k)==1)then;z(j,k)=rtruncnorm(m,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else if(response(j,k)==0)then;z(j,k)=rtruncnorm(m,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
            else;z(j,k)=rnorm(m,1.0_dp);end if
         end do;end do
         ! Item parameters eta=(alpha,beta), with design (-1,theta).
         tpt=0.0_dp;tpt(1,1)=real(nj,dp);tpt(1,2)=-sum(theta);tpt(2,1)=tpt(1,2);tpt(2,2)=dot_product(theta,theta)
         call inv_spd(ab0prec+tpt,postv,info);if(info/=0)then;res%status=10+info;return;end if
         do k=1,nk
            rhs=matmul(ab0prec,ab0);rhs(1)=rhs(1)-sum(z(:,k));rhs(2)=rhs(2)+dot_product(theta,z(:,k));postm=matmul(postv,rhs)
            call rmvnorm(postm,postv,draw,info);if(info/=0)then;res%status=20+info;return;end if;eta(k,:)=draw
         end do
         ! Ability parameters.
         v=1.0_dp/(t0prec+dot_product(eta(:,2),eta(:,2)))
         do j=1,nj
            if(theta_eq(j)>-998.5_dp)then;theta(j)=theta_eq(j)
            else
               bsum=dot_product(eta(:,2),z(j,:)+eta(:,1));m=v*(t0prec*t0+bsum)
               if(theta_ineq(j)>0.0_dp)then;theta(j)=rtruncnorm(m,sqrt(v),0.0_dp,huge(1.0_dp)/10.0_dp)
               else if(theta_ineq(j)<0.0_dp)then;theta(j)=rtruncnorm(m,sqrt(v),-huge(1.0_dp)/10.0_dp,0.0_dp)
               else;theta(j)=rnorm(m,sqrt(v));end if
            end if
         end do
         if(iter>=burnin.and.mod(iter,thin)==0)then
            count=count+1;col=0
            if(sa)then;res%draws(count,col+1:col+nj)=theta;col=col+nj;end if
            if(si)then
               do k=1,nk;res%draws(count,col+2*k-1)=eta(k,1);res%draws(count,col+2*k)=eta(k,2);end do
            end if
         end if
      end do
   end function mcmc_irt1d

   function mcmc_paircompare(md,theta_start,alpha_start,theta_eq,theta_ineq,a0,a0prec,alpha_fixed, &
                             burnin,mcmc,thin,store_theta,store_alpha) result(res)
      ! md(:,1)=respondent/judge, md(:,2)=candidate 1, md(:,3)=candidate 2, md(:,4)=chosen candidate.
      ! Indices are Fortran 1-based in this API.
      integer,intent(in)::md(:,:)
      real(dp),intent(in)::theta_start(:),alpha_start(:),theta_eq(:),theta_ineq(:),a0,a0prec
      logical,intent(in)::alpha_fixed
      integer,intent(in)::burnin,mcmc,thin
      logical,intent(in),optional::store_theta,store_alpha
      type(mcmc_result)::res
      integer::n,jj,ii,r,c1,c2,iter,count,nsamp,j,i,ncol,col
      real(dp)::theta(size(theta_start)),alpha(size(alpha_start)),ystar(size(md,1)),mu,xx,xz,v,m,sgn,xval,zval
      logical::st,sa
      st=.true.;sa=.false.;if(present(store_theta))st=store_theta;if(present(store_alpha))sa=store_alpha
      n=size(md,1);jj=size(theta);ii=size(alpha);nsamp=mcmc/thin
      if(size(md,2)<4.or.size(theta_eq)/=jj.or.size(theta_ineq)/=jj.or. &
         nsamp<=0.or.(.not.st.and..not.sa))then
         res%status=1;return
      end if
      if(any(md(:,1)<1).or.any(md(:,1)>ii).or.any(md(:,2)<1).or. &
         any(md(:,2)>jj).or.any(md(:,3)<1).or.any(md(:,3)>jj))then
         res%status=2;return
      end if
      ncol=merge(jj,0,st)+merge(ii,0,sa);allocate(res%draws(nsamp,ncol));theta=theta_start;alpha=alpha_start;count=0
      do iter=0,burnin+mcmc-1
         do i=1,n
            r=md(i,1);c1=md(i,2);c2=md(i,3);mu=alpha(r)*(theta(c1)-theta(c2))
            if(md(i,4)==c1)then;ystar(i)=rtruncnorm(mu,1.0_dp,0.0_dp,huge(1.0_dp)/10.0_dp)
            else if(md(i,4)==c2)then;ystar(i)=rtruncnorm(mu,1.0_dp,-huge(1.0_dp)/10.0_dp,0.0_dp)
            else;ystar(i)=rnorm(mu,1.0_dp);end if
         end do
         do j=1,jj
            if(theta_eq(j)>-998.5_dp)then;theta(j)=theta_eq(j);cycle;end if
            xx=1.0_dp;xz=0.0_dp
            do i=1,n
               r=md(i,1);c1=md(i,2);c2=md(i,3)
               if(c1==j)then;sgn=1.0_dp;xval=alpha(r);zval=ystar(i)+alpha(r)*theta(c2);xx=xx+xval*xval;xz=xz+xval*zval
               else if(c2==j)then;sgn=-1.0_dp;xval=-alpha(r);zval=ystar(i)-alpha(r)*theta(c1);xx=xx+xval*xval;xz=xz+xval*zval
               end if
            end do
            v=1.0_dp/xx;m=v*xz
            if(theta_ineq(j)>0.0_dp)then;theta(j)=rtruncnorm(m,sqrt(v),0.0_dp,huge(1.0_dp)/10.0_dp)
            else if(theta_ineq(j)<0.0_dp)then;theta(j)=rtruncnorm(m,sqrt(v),-huge(1.0_dp)/10.0_dp,0.0_dp)
            else;theta(j)=rnorm(m,sqrt(v));end if
         end do
         if(.not.alpha_fixed)then
            do r=1,ii
               xx=a0prec;xz=a0prec*a0
               do i=1,n
                  if(md(i,1)/=r)cycle;c1=md(i,2);c2=md(i,3);xval=theta(c1)-theta(c2);xx=xx+xval*xval;xz=xz+xval*ystar(i)
               end do
               v=1.0_dp/xx;m=v*xz;alpha(r)=rnorm(m,sqrt(v))
            end do
         end if
         if(iter>=burnin.and.mod(iter,thin)==0)then
            count=count+1;col=0;if(st)then;res%draws(count,1:jj)=theta;col=jj;end if;if(sa)res%draws(count,col+1:col+ii)=alpha
         end if
      end do
   end function mcmc_paircompare
end module mcmcpack_latent
