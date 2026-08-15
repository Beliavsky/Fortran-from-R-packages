! Simultaneous random intercepts on multiple GAMLSS distribution parameters.
! Each active parameter receives its own group-effect vector and independently
! estimated quadratic-penalty multiplier during the RS iterations.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_multi_random_v04
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   implicit none
   private
   public :: multi_random_intercept_result_t, fit_gamlss_multi_random_intercept

   type,public :: multi_random_intercept_result_t
      type(gamlss_result_t) :: model
      real(dp),allocatable :: effects(:,:)
      real(dp) :: lambda(4)=0.0_dp
      real(dp) :: working_variance_ratio(4)=0.0_dp
      logical :: active(4)=.false.
      integer,allocatable :: levels(:)
      integer :: status=0
   end type multi_random_intercept_result_t
contains

   subroutine fit_gamlss_multi_random_intercept(y,x_mu,group,family,result,active_parameters, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control)
      real(dp),intent(in)::y(:),x_mu(:,:)
      integer,intent(in)::group(:),family
      type(multi_random_intercept_result_t),intent(out)::result
      logical,intent(in),optional::active_parameters(4)
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional::control
      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:),z(:,:)
      real(dp),allocatable::pm(:,:),ps(:,:),pn(:,:),pt(:,:)
      type(gamlss_control_t)::ctl
      logical::act(4)
      integer::n,np,ng,istat,pfix(4),j

      n=size(y);np=family_npar(family)
      if(n<=0.or.size(x_mu,1)/=n.or.size(group)/=n.or.np<1.or.np>4)then;result%status=1;return;end if
      act=.false.;act(1:np)=.true.
      if(present(active_parameters))act=active_parameters
      if(any(act(np+1:4)))then;result%status=2;return;end if
      if(.not.any(act(1:np)))then;result%status=3;return;end if
      call compressed_design(group,z,result%levels,istat)
      if(istat/=0)then;result%status=4;return;end if
      ng=size(z,2)
      call default_design(n,x_mu,x_sigma,x_nu,x_tau,xm,xs,xn,xt)
      pfix=[size(xm,2),size(xs,2),size(xn,2),size(xt,2)]
      call prepare_block(xm,z,act(1),pm)
      call prepare_block(xs,z,act(2),ps)
      call prepare_block(xn,z,act(3),pn)
      call prepare_block(xt,z,act(4),pt)
      ctl=gamlss_control_t();if(present(control))ctl=control
      ctl%estimate_lambda_mu=act(1);ctl%estimate_lambda_sigma=act(2)
      ctl%estimate_lambda_nu=act(3);ctl%estimate_lambda_tau=act(4)
      call fit_gamlss_model(y,xm,family,result%model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
         offset_mu=offset_mu,offset_sigma=offset_sigma,offset_nu=offset_nu,offset_tau=offset_tau, &
         penalty_mu=pm,penalty_sigma=ps,penalty_nu=pn,penalty_tau=pt, &
         lambda_mu=merge(1.0_dp,0.0_dp,act(1)),lambda_sigma=merge(1.0_dp,0.0_dp,act(2)), &
         lambda_nu=merge(1.0_dp,0.0_dp,act(3)),lambda_tau=merge(1.0_dp,0.0_dp,act(4)),control=ctl)
      result%status=result%model%status;result%active=act
      if(result%status/=0)return
      allocate(result%effects(ng,4));result%effects=0.0_dp
      do j=1,np
         if(.not.act(j))cycle
         select case(j)
         case(1)
            result%lambda(j)=result%model%mu%lambda
            result%effects(:,j)=result%model%mu%coefficients(pfix(j)+1:pfix(j)+ng)
         case(2)
            result%lambda(j)=result%model%sigma%lambda
            result%effects(:,j)=result%model%sigma%coefficients(pfix(j)+1:pfix(j)+ng)
         case(3)
            result%lambda(j)=result%model%nu%lambda
            result%effects(:,j)=result%model%nu%coefficients(pfix(j)+1:pfix(j)+ng)
         case(4)
            result%lambda(j)=result%model%tau%lambda
            result%effects(:,j)=result%model%tau%coefficients(pfix(j)+1:pfix(j)+ng)
         end select
         result%working_variance_ratio(j)=1.0_dp/max(1.0e-12_dp,result%lambda(j))
      end do
   end subroutine fit_gamlss_multi_random_intercept

   subroutine default_design(n,xmu,xs_in,xn_in,xt_in,xm,xs,xn,xt)
      integer,intent(in)::n
      real(dp),intent(in)::xmu(:,:)
      real(dp),intent(in),optional::xs_in(:,:),xn_in(:,:),xt_in(:,:)
      real(dp),allocatable,intent(out)::xm(:,:),xs(:,:),xn(:,:),xt(:,:)
      xm=xmu
      if(present(xs_in))then;xs=xs_in;else;allocate(xs(n,1));xs=1.0_dp;end if
      if(present(xn_in))then;xn=xn_in;else;allocate(xn(n,1));xn=1.0_dp;end if
      if(present(xt_in))then;xt=xt_in;else;allocate(xt(n,1));xt=1.0_dp;end if
   end subroutine default_design

   subroutine prepare_block(x,z,active,penalty)
      real(dp),allocatable,intent(inout)::x(:,:)
      real(dp),intent(in)::z(:,:)
      logical,intent(in)::active
      real(dp),allocatable,intent(out)::penalty(:,:)
      real(dp),allocatable::tmp(:,:)
      integer::p,k
      p=size(x,2)
      if(active)then
         allocate(tmp(size(x,1),p+size(z,2)));tmp(:,1:p)=x;tmp(:,p+1:)=z;call move_alloc(tmp,x)
      end if
      allocate(penalty(size(x,2),size(x,2)));penalty=0.0_dp
      if(active)then
         do k=p+1,size(x,2);penalty(k,k)=1.0_dp;end do
      end if
   end subroutine prepare_block

   subroutine compressed_design(group,z,levels,status)
      integer,intent(in)::group(:)
      real(dp),allocatable,intent(out)::z(:,:)
      integer,allocatable,intent(out)::levels(:)
      integer,intent(out)::status
      integer,allocatable::tmp(:)
      integer::i,j,nlev
      if(size(group)==0)then;allocate(z(0,0),levels(0));status=1;return;end if
      allocate(tmp(size(group)));nlev=0
      do i=1,size(group)
         if(nlev==0)then
            nlev=1;tmp(1)=group(i)
         else if(.not.any(tmp(1:nlev)==group(i)))then
            nlev=nlev+1;tmp(nlev)=group(i)
         end if
      end do
      allocate(levels(nlev));levels=tmp(1:nlev);allocate(z(size(group),nlev));z=0.0_dp
      do i=1,size(group)
         do j=1,nlev
            if(group(i)==levels(j))then;z(i,j)=1.0_dp;exit;end if
         end do
      end do
      status=0
   end subroutine compressed_design
end module gamlss_multi_random_v04
