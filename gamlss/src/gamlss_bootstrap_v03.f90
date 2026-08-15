! Case-resampling bootstrap and profile-likelihood confidence intervals.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_bootstrap_v03
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   use gamlss_model_selection_v02, only : profile_result_t
   use gamlss_special, only : normal_quantile
   implicit none
   private
   public :: gamlss_bootstrap_result_t, bootstrap_gamlss_cases, bootstrap_percentile_ci
   public :: profile_likelihood_ci

   type, public :: gamlss_bootstrap_result_t
      real(dp), allocatable :: beta_mu(:,:),beta_sigma(:,:),beta_nu(:,:),beta_tau(:,:)
      real(dp), allocatable :: deviance(:)
      integer, allocatable :: status(:)
      integer :: successful = 0
   end type gamlss_bootstrap_result_t
contains

   subroutine bootstrap_gamlss_cases(y,x_mu,family,n_boot,result,x_sigma,x_nu,x_tau,weights, &
      offset_mu,offset_sigma,offset_nu,offset_tau,penalty_mu,penalty_sigma,penalty_nu,penalty_tau, &
      lambda_mu,lambda_sigma,lambda_nu,lambda_tau,control)
      real(dp),intent(in)::y(:),x_mu(:,:)
      integer,intent(in)::family,n_boot
      type(gamlss_bootstrap_result_t),intent(out)::result
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      real(dp),intent(in),optional::penalty_mu(:,:),penalty_sigma(:,:),penalty_nu(:,:),penalty_tau(:,:)
      real(dp),intent(in),optional::lambda_mu,lambda_sigma,lambda_nu,lambda_tau
      type(gamlss_control_t),intent(in),optional::control
      integer,allocatable::idx(:)
      real(dp),allocatable::yb(:),xm(:,:),xs(:,:),xn(:,:),xt(:,:),wb(:),om(:),os(:),on(:),ot(:)
      type(gamlss_result_t)::fit
      integer::n,np,b,i,p1,p2,p3,p4
      real(dp)::u

      n=size(y);np=family_npar(family);p1=size(x_mu,2);p2=0;p3=0;p4=0
      if(n<1.or.size(x_mu,1)/=n.or.n_boot<1)then;allocate(result%status(0));return;end if
      if(np>=2)then;if(.not.present(x_sigma))then;allocate(result%status(0));return;end if;p2=size(x_sigma,2);end if
      if(np>=3)then;if(.not.present(x_nu))then;allocate(result%status(0));return;end if;p3=size(x_nu,2);end if
      if(np>=4)then;if(.not.present(x_tau))then;allocate(result%status(0));return;end if;p4=size(x_tau,2);end if
      allocate(result%beta_mu(p1,n_boot),result%deviance(n_boot),result%status(n_boot));result%beta_mu=0.0_dp
      if(p2>0)allocate(result%beta_sigma(p2,n_boot));if(p3>0)allocate(result%beta_nu(p3,n_boot))
      if(p4>0)allocate(result%beta_tau(p4,n_boot))
      result%deviance=huge(1.0_dp);result%status=1;result%successful=0
      allocate(idx(n),yb(n),xm(n,p1),wb(n),om(n),os(n),on(n),ot(n))
      allocate(xs(n,max(1,p2)),xn(n,max(1,p3)),xt(n,max(1,p4)))
      wb=1.0_dp;om=0.0_dp;os=0.0_dp;on=0.0_dp;ot=0.0_dp
      do b=1,n_boot
         do i=1,n
            call random_number(u);idx(i)=min(n,1+int(u*real(n,dp)))
         end do
         yb=y(idx);xm=x_mu(idx,:)
         xs=1.0_dp;xn=1.0_dp;xt=1.0_dp
         if(p2>0)xs(:,1:p2)=x_sigma(idx,:);if(p3>0)xn(:,1:p3)=x_nu(idx,:);if(p4>0)xt(:,1:p4)=x_tau(idx,:)
         wb=1.0_dp;om=0.0_dp;os=0.0_dp;on=0.0_dp;ot=0.0_dp
         if(present(weights))wb=weights(idx)
         if(present(offset_mu))om=offset_mu(idx);if(present(offset_sigma))os=offset_sigma(idx)
         if(present(offset_nu))on=offset_nu(idx);if(present(offset_tau))ot=offset_tau(idx)
         call fit_gamlss_model(yb,xm,family,fit,x_sigma=xs(:,1:max(1,p2)), &
            x_nu=xn(:,1:max(1,p3)),x_tau=xt(:,1:max(1,p4)),weights=wb, &
            offset_mu=om,offset_sigma=os,offset_nu=on,offset_tau=ot, &
            penalty_mu=penalty_mu,penalty_sigma=penalty_sigma,penalty_nu=penalty_nu,penalty_tau=penalty_tau, &
            lambda_mu=lambda_mu,lambda_sigma=lambda_sigma,lambda_nu=lambda_nu,lambda_tau=lambda_tau,control=control)
         result%status(b)=fit%status
         if(fit%status==0)then
            result%beta_mu(:,b)=fit%mu%coefficients
            if(p2>0)result%beta_sigma(:,b)=fit%sigma%coefficients
            if(p3>0)result%beta_nu(:,b)=fit%nu%coefficients
            if(p4>0)result%beta_tau(:,b)=fit%tau%coefficients
            result%deviance(b)=fit%global_deviance;result%successful=result%successful+1
         end if
      end do
   end subroutine bootstrap_gamlss_cases

   subroutine bootstrap_percentile_ci(samples,level,lower,upper,status)
      real(dp),intent(in)::samples(:,:),level
      real(dp),allocatable,intent(out)::lower(:),upper(:)
      integer,intent(out),optional::status
      real(dp),allocatable::v(:)
      integer::j,n,ilo,ihi
      if(size(samples,2)<2.or.level<=0.0_dp.or.level>=1.0_dp)then
         allocate(lower(0),upper(0));if(present(status))status=1;return
      end if
      n=size(samples,2);ilo=max(1,int(floor(0.5_dp*(1.0_dp-level)*real(n-1,dp)))+1)
      ihi=min(n,int(ceiling((1.0_dp-0.5_dp*(1.0_dp-level))*real(n-1,dp)))+1)
      allocate(lower(size(samples,1)),upper(size(samples,1)))
      do j=1,size(samples,1)
         v=samples(j,:);call sort_real(v);lower(j)=v(ilo);upper(j)=v(ihi)
      end do
      if(present(status))status=0
   end subroutine bootstrap_percentile_ci

   subroutine profile_likelihood_ci(profile,level,lower,upper,status)
      type(profile_result_t),intent(in)::profile
      real(dp),intent(in)::level
      real(dp),intent(out)::lower,upper
      integer,intent(out),optional::status
      real(dp)::cut,z,maxll
      integer::imax,i,istat
      lower=0.0_dp;upper=0.0_dp;istat=0
      if(size(profile%value)<3.or.size(profile%loglik)/=size(profile%value).or.level<=0.0_dp.or.level>=1.0_dp)then
         istat=1;if(present(status))status=istat;return
      end if
      imax=maxloc(profile%loglik,dim=1);maxll=profile%loglik(imax)
      z=normal_quantile(0.5_dp*(1.0_dp+level));cut=maxll-0.5_dp*z*z
      if(imax==1.or.imax==size(profile%value))then;istat=2;if(present(status))status=istat;return;end if
      lower=profile%value(1);upper=profile%value(size(profile%value))
      do i=imax-1,1,-1
         if(profile%loglik(i)<=cut)then
            lower=linear_cross(profile%value(i),profile%loglik(i),profile%value(i+1),profile%loglik(i+1),cut);exit
         end if
      end do
      do i=imax+1,size(profile%value)
         if(profile%loglik(i)<=cut)then
            upper=linear_cross(profile%value(i-1),profile%loglik(i-1),profile%value(i),profile%loglik(i),cut);exit
         end if
      end do
      if(present(status))status=istat
   end subroutine profile_likelihood_ci

   real(dp) function linear_cross(x1,y1,x2,y2,y) result(x)
      real(dp),intent(in)::x1,y1,x2,y2,y
      if(abs(y2-y1)<1.0e-15_dp)then;x=0.5_dp*(x1+x2);else;x=x1+(x2-x1)*(y-y1)/(y2-y1);end if
   end function linear_cross

   subroutine sort_real(a)
      real(dp),intent(inout)::a(:)
      real(dp)::key
      integer::i,j
      do i=2,size(a)
         key=a(i);j=i-1
         do while(j>=1)
            if(a(j)<=key)exit
            a(j+1)=a(j);j=j-1
         end do
         a(j+1)=key
      end do
   end subroutine sort_real
end module gamlss_bootstrap_v03
