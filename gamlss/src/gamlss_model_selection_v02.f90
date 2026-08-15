! Matrix-based stepwise GAIC and coefficient profile-likelihood utilities.
! These preserve the computational part of stepGAIC/profile workflows without R formulas.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_model_selection_v02
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   implicit none
   private
   public :: stepwise_result_t, forward_gaic_mu
   public :: profile_result_t, profile_gamlss_coefficient

   type,public :: stepwise_result_t
      integer,allocatable :: selected(:)
      real(dp),allocatable :: gaic_path(:)
      type(gamlss_result_t) :: final_model
      integer :: status=0
   end type stepwise_result_t

   type,public :: profile_result_t
      real(dp),allocatable :: value(:),loglik(:),deviance(:)
      integer,allocatable :: status(:)
      integer :: parameter=1
      integer :: coefficient=1
   end type profile_result_t
contains

   subroutine forward_gaic_mu(y,x_base,x_candidates,family,result,k,x_sigma,x_nu,x_tau,weights,control)
      real(dp),intent(in)::y(:),x_base(:,:),x_candidates(:,:)
      integer,intent(in)::family
      type(stepwise_result_t),intent(out)::result
      real(dp),intent(in),optional::k,x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      type(gamlss_control_t),intent(in),optional::control
      logical,allocatable::chosen(:)
      integer,allocatable::sel(:),tmp_sel(:)
      real(dp),allocatable::path(:),tmp_path(:),xm(:,:)
      type(gamlss_result_t)::fit,bestfit
      real(dp)::kk,current,best,val
      integer::j,bestj,step,ncand,istat
      if(size(x_base,1)/=size(y).or.size(x_candidates,1)/=size(y))then;result%status=1;return;end if
      kk=2.0_dp;if(present(k))kk=k;ncand=size(x_candidates,2);allocate(chosen(ncand));chosen=.false.
      allocate(sel(0),path(0));call build_mu(x_base,x_candidates,sel,xm)
      call fit_model(xm,fit,istat);if(istat/=0)then;result%status=2;return;end if
      current=fit%global_deviance+kk*fit%df_fit;bestfit=fit
      do step=1,ncand
         best=current;bestj=0
         do j=1,ncand
            if(chosen(j))cycle
            call append_int(sel,j,tmp_sel);call build_mu(x_base,x_candidates,tmp_sel,xm)
            call fit_model(xm,fit,istat);if(istat/=0)cycle
            val=fit%global_deviance+kk*fit%df_fit
            if(val<best-1.0e-8_dp)then;best=val;bestj=j;bestfit=fit;end if
         end do
         if(bestj==0)exit
         chosen(bestj)=.true.;call append_int(sel,bestj,tmp_sel);call move_alloc(tmp_sel,sel)
         call append_real(path,best,tmp_path);call move_alloc(tmp_path,path);current=best
      end do
      result%selected=sel;result%gaic_path=path;result%final_model=bestfit;result%status=0
   contains
      subroutine fit_model(xmu,model,status)
         real(dp),intent(in)::xmu(:,:)
         type(gamlss_result_t),intent(out)::model
         integer,intent(out)::status
         real(dp),allocatable::xs(:,:),xn(:,:),xt(:,:)
         integer::np,n
         n=size(y);np=family_npar(family)
         call optional_design(n,x_sigma,xs);call optional_design(n,x_nu,xn);call optional_design(n,x_tau,xt)
         call fit_gamlss_model(y,xmu,family,model,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights,control=control)
         status=model%status
      end subroutine fit_model
   end subroutine forward_gaic_mu

   subroutine profile_gamlss_coefficient(y,family,parameter,coefficient,grid,profile,x_mu, &
      x_sigma,x_nu,x_tau,weights,offset_mu,offset_sigma,offset_nu,offset_tau,control)
      real(dp),intent(in)::y(:),grid(:),x_mu(:,:)
      integer,intent(in)::family,parameter,coefficient
      type(profile_result_t),intent(out)::profile
      real(dp),intent(in),optional::x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      real(dp),intent(in),optional::offset_mu(:),offset_sigma(:),offset_nu(:),offset_tau(:)
      type(gamlss_control_t),intent(in),optional::control
      real(dp),allocatable::xm(:,:),xs(:,:),xn(:,:),xt(:,:),om(:),os(:),on(:),ot(:)
      real(dp),allocatable::xnew(:,:),offnew(:)
      logical::fixed
      type(gamlss_result_t)::fit
      integer::i,np,p,istat,n
      n=size(y);np=family_npar(family)
      if(parameter<1.or.parameter>np.or.size(x_mu,1)/=n)then;allocate(profile%value(0));return;end if
      xm=x_mu;call optional_design(n,x_sigma,xs);call optional_design(n,x_nu,xn);call optional_design(n,x_tau,xt)
      call optional_offset(n,offset_mu,om);call optional_offset(n,offset_sigma,os)
      call optional_offset(n,offset_nu,on);call optional_offset(n,offset_tau,ot)
      select case(parameter)
      case(1);p=size(xm,2)
      case(2);p=size(xs,2)
      case(3);p=size(xn,2)
      case(4);p=size(xt,2)
      end select
      if(coefficient<1.or.coefficient>p)then;allocate(profile%value(0));return;end if
      allocate(profile%value(size(grid)),profile%loglik(size(grid)),profile%deviance(size(grid)),profile%status(size(grid)))
      profile%value=grid;profile%parameter=parameter;profile%coefficient=coefficient
      do i=1,size(grid)
         select case(parameter)
         case(1)
            call remove_fixed_column(xm,om,coefficient,grid(i),xnew,offnew,fixed);xm=xnew;om=offnew
         case(2)
            call remove_fixed_column(xs,os,coefficient,grid(i),xnew,offnew,fixed);xs=xnew;os=offnew
         case(3)
            call remove_fixed_column(xn,on,coefficient,grid(i),xnew,offnew,fixed);xn=xnew;on=offnew
         case(4)
            call remove_fixed_column(xt,ot,coefficient,grid(i),xnew,offnew,fixed);xt=xnew;ot=offnew
         end select
         call fit_gamlss_model(y,xm,family,fit,x_sigma=xs,x_nu=xn,x_tau=xt,weights=weights, &
            offset_mu=om,offset_sigma=os,offset_nu=on,offset_tau=ot, &
            fix_mu=(fixed.and.parameter==1),fix_sigma=(fixed.and.parameter==2), &
            fix_nu=(fixed.and.parameter==3),fix_tau=(fixed.and.parameter==4),control=control)
         profile%status(i)=fit%status
         if(fit%status==0)then
            profile%deviance(i)=fit%global_deviance;profile%loglik(i)=-0.5_dp*fit%global_deviance
         else
            profile%deviance(i)=huge(1.0_dp);profile%loglik(i)=-huge(1.0_dp)
         end if
         ! Restore full designs/offsets before the next grid value.
         xm=x_mu;call optional_design(n,x_sigma,xs);call optional_design(n,x_nu,xn);call optional_design(n,x_tau,xt)
         call optional_offset(n,offset_mu,om);call optional_offset(n,offset_sigma,os)
         call optional_offset(n,offset_nu,on);call optional_offset(n,offset_tau,ot)
      end do
   end subroutine profile_gamlss_coefficient

   subroutine remove_fixed_column(x,offset,j,value,xnew,offnew,all_fixed)
      real(dp),intent(in)::x(:,:),offset(:),value
      integer,intent(in)::j
      real(dp),allocatable,intent(out)::xnew(:,:),offnew(:)
      logical,intent(out)::all_fixed
      integer::p
      p=size(x,2);allocate(offnew(size(offset)));offnew=offset+value*x(:,j)
      if(p==1)then
         allocate(xnew(size(x,1),1));xnew=0.0_dp;all_fixed=.true.
      else
         allocate(xnew(size(x,1),p-1));
         if(j>1)xnew(:,1:j-1)=x(:,1:j-1)
         if(j<p)xnew(:,j:p-1)=x(:,j+1:p)
         all_fixed=.false.
      end if
   end subroutine remove_fixed_column

   subroutine optional_design(n,xin,xout)
      integer,intent(in)::n
      real(dp),intent(in),optional::xin(:,:)
      real(dp),allocatable,intent(out)::xout(:,:)
      if(present(xin))then;xout=xin;else;allocate(xout(n,1));xout=1.0_dp;end if
   end subroutine optional_design

   subroutine optional_offset(n,oin,oout)
      integer,intent(in)::n
      real(dp),intent(in),optional::oin(:)
      real(dp),allocatable,intent(out)::oout(:)
      allocate(oout(n));oout=0.0_dp;if(present(oin))oout=oin
   end subroutine optional_offset

   subroutine build_mu(base,candidates,sel,x)
      real(dp),intent(in)::base(:,:),candidates(:,:)
      integer,intent(in)::sel(:)
      real(dp),allocatable,intent(out)::x(:,:)
      integer::j,p
      p=size(base,2);allocate(x(size(base,1),p+size(sel)));x(:,1:p)=base
      do j=1,size(sel);x(:,p+j)=candidates(:,sel(j));end do
   end subroutine build_mu

   subroutine append_int(a,v,b)
      integer,intent(in)::a(:),v
      integer,allocatable,intent(out)::b(:)
      allocate(b(size(a)+1));if(size(a)>0)b(1:size(a))=a;b(size(b))=v
   end subroutine append_int

   subroutine append_real(a,v,b)
      real(dp),intent(in)::a(:),v
      real(dp),allocatable,intent(out)::b(:)
      allocate(b(size(a)+1));if(size(a)>0)b(1:size(a))=a;b(size(b))=v
   end subroutine append_real

end module gamlss_model_selection_v02
