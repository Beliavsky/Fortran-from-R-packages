! Random-intercept GAMLSS adapter using quadratic-penalty/variance-ratio updates.
! For Gaussian location models, nlme is optionally used for the starting variance ratio.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_random_effects
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar, GAMLSS_NO
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   use nlme_lme, only : fit_lme
   use nlme_types, only : lme_result
   implicit none
   private
   public :: random_intercept_result_t, fit_gamlss_random_intercept

   type,public :: random_intercept_result_t
      type(gamlss_result_t) :: model
      real(dp),allocatable :: effects(:)
      integer,allocatable :: levels(:)
      integer :: parameter=1
      real(dp) :: lambda=1.0_dp
      real(dp) :: working_variance_ratio=1.0_dp
      real(dp) :: random_edf=0.0_dp
      real(dp) :: sigma_b=0.0_dp
      integer :: status=0
   end type random_intercept_result_t
contains

   subroutine fit_gamlss_random_intercept(y,x_fixed,group,family,result,parameter, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control,use_nlme_start)
      real(dp),intent(in)::y(:),x_fixed(:,:)
      integer,intent(in)::group(:),family
      type(random_intercept_result_t),intent(out)::result
      integer,intent(in),optional::parameter
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional::control
      logical,intent(in),optional::use_nlme_start
      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:),z(:,:)
      real(dp),allocatable::pm(:,:),ps(:,:),pn(:,:),pt(:,:),zri(:,:)
      type(gamlss_control_t)::ctl
      type(lme_result)::lr
      integer::n,np,target,pfix,ng,istat
      real(dp)::lambda0
      logical::nlstart

      n=size(y);np=family_npar(family);target=1;if(present(parameter))target=parameter
      if(n<=0.or.size(x_fixed,1)/=n.or.size(group)/=n.or.target<1.or.target>np)then
         result%status=1;return
      end if
      call compressed_random_design(group,z,result%levels,istat)
      if(istat/=0)then;result%status=2;return;end if
      ng=size(z,2)
      call default_design(n,x_fixed,x_sigma,x_nu,x_tau,xm,xs,xn,xt)
      select case(target)
      case(1);pfix=size(xm,2);call augment(xm,z,pm);call zero_penalty(xs,ps);call zero_penalty(xn,pn);call zero_penalty(xt,pt)
      case(2);pfix=size(xs,2);call augment(xs,z,ps);call zero_penalty(xm,pm);call zero_penalty(xn,pn);call zero_penalty(xt,pt)
      case(3);pfix=size(xn,2);call augment(xn,z,pn);call zero_penalty(xm,pm);call zero_penalty(xs,ps);call zero_penalty(xt,pt)
      case(4);pfix=size(xt,2);call augment(xt,z,pt);call zero_penalty(xm,pm);call zero_penalty(xs,ps);call zero_penalty(xn,pn)
      end select
      lambda0=1.0_dp;nlstart=.true.;if(present(use_nlme_start))nlstart=use_nlme_start
      if(nlstart.and.family==GAMLSS_NO.and.target==1)then
         allocate(zri(n,1));zri=1.0_dp
         call fit_lme(y,x_fixed,zri,group,lr)
         if(lr%status==0.and.allocated(lr%random_covariance))then
            if(lr%random_covariance(1,1)>1.0e-12_dp.and.lr%sigma>0.0_dp) &
               lambda0=lr%sigma*lr%sigma/lr%random_covariance(1,1)
         end if
      end if
      ctl=gamlss_control_t();if(present(control))ctl=control
      select case(target)
      case(1);ctl%estimate_lambda_mu=.true.
      case(2);ctl%estimate_lambda_sigma=.true.
      case(3);ctl%estimate_lambda_nu=.true.
      case(4);ctl%estimate_lambda_tau=.true.
      end select
      call fit_gamlss_model(y,xm,family,result%model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
         offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
         penalty_mu=pm,penalty_sigma=ps,penalty_nu=pn,penalty_tau=pt, &
         lambda_mu=merge(lambda0,0.0_dp,target==1),lambda_sigma=merge(lambda0,0.0_dp,target==2), &
         lambda_nu=merge(lambda0,0.0_dp,target==3),lambda_tau=merge(lambda0,0.0_dp,target==4),control=ctl)
      result%status=result%model%status;result%parameter=target
      if(result%status/=0)return
      select case(target)
      case(1);result%lambda=result%model%mu%lambda;result%effects=result%model%mu%coefficients(pfix+1:pfix+ng)
      case(2);result%lambda=result%model%sigma%lambda;result%effects=result%model%sigma%coefficients(pfix+1:pfix+ng)
      case(3);result%lambda=result%model%nu%lambda;result%effects=result%model%nu%coefficients(pfix+1:pfix+ng)
      case(4);result%lambda=result%model%tau%lambda;result%effects=result%model%tau%coefficients(pfix+1:pfix+ng)
      end select
      result%working_variance_ratio=1.0_dp/max(1.0e-12_dp,result%lambda)
      call random_scale_summary(result%model,target,pfix,ng,result%lambda,result%effects, &
         result%random_edf,result%sigma_b)
   contains

      subroutine random_scale_summary(model,target,p0,nr,lambda,effect,edf,sigb)
         type(gamlss_result_t),intent(in)::model
         integer,intent(in)::target,p0,nr
         real(dp),intent(in)::lambda,effect(:)
         real(dp),intent(out)::edf,sigb
         real(dp),allocatable::cov(:,:)
         integer::j
         select case(target)
         case(1);cov=model%mu%covariance
         case(2);cov=model%sigma%covariance
         case(3);cov=model%nu%covariance
         case(4);cov=model%tau%covariance
         end select
         edf=real(nr,dp)
         if(size(cov,1)>=p0+nr)then
            do j=1,nr
               edf=edf-lambda*cov(p0+j,p0+j)
            end do
         end if
         edf=max(edf,1.0e-10_dp)
         sigb=sqrt(max(0.0_dp,sum(effect*effect)/edf))
      end subroutine random_scale_summary

      subroutine default_design(nn,xf,xs_in,xn_in,xt_in,a,b,c,d)
         integer,intent(in)::nn
         real(dp),intent(in)::xf(:,:)
         real(dp),intent(in),optional::xs_in(:,:),xn_in(:,:),xt_in(:,:)
         real(dp),allocatable,intent(out)::a(:,:),b(:,:),c(:,:),d(:,:)
         a=xf
         if(present(xs_in))then;b=xs_in;else;allocate(b(nn,1));b=1.0_dp;end if
         if(present(xn_in))then;c=xn_in;else;allocate(c(nn,1));c=1.0_dp;end if
         if(present(xt_in))then;d=xt_in;else;allocate(d(nn,1));d=1.0_dp;end if
      end subroutine default_design

      subroutine augment(a,zmat,pen)
         real(dp),allocatable,intent(inout)::a(:,:)
         real(dp),intent(in)::zmat(:,:)
         real(dp),allocatable,intent(out)::pen(:,:)
         real(dp),allocatable::tmp(:,:)
         integer::p,j
         p=size(a,2);allocate(tmp(size(a,1),p+size(zmat,2)));tmp(:,1:p)=a;tmp(:,p+1:)=zmat
         call move_alloc(tmp,a);allocate(pen(size(a,2),size(a,2)));pen=0.0_dp
         do j=p+1,size(a,2);pen(j,j)=1.0_dp;end do
      end subroutine augment

      subroutine zero_penalty(a,pen)
         real(dp),intent(in)::a(:,:)
         real(dp),allocatable,intent(out)::pen(:,:)
         allocate(pen(size(a,2),size(a,2)));pen=0.0_dp
      end subroutine zero_penalty
   end subroutine fit_gamlss_random_intercept

   subroutine compressed_random_design(group,z,levels,status)
      integer,intent(in)::group(:)
      real(dp),allocatable,intent(out)::z(:,:)
      integer,allocatable,intent(out)::levels(:)
      integer,intent(out)::status
      integer,allocatable::tmp(:)
      integer::i,j,nlev
      if(size(group)==0)then;allocate(z(0,0),levels(0));status=1;return;end if
      allocate(tmp(size(group)));nlev=0
      do i=1,size(group)
         if(.not.any(tmp(1:nlev)==group(i)))then;nlev=nlev+1;tmp(nlev)=group(i);end if
      end do
      allocate(levels(nlev));levels=tmp(1:nlev);allocate(z(size(group),nlev));z=0.0_dp
      do i=1,size(group)
         do j=1,nlev;if(group(i)==levels(j))then;z(i,j)=1.0_dp;exit;end if;end do
      end do
      status=0
   end subroutine compressed_random_design

end module gamlss_random_effects
