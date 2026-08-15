! Backward and bidirectional matrix-based GAIC selection.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_model_selection_v03
   use gamlss_kinds, only : dp
   use gamlss_fit, only : family_npar
   use gamlss_types
   use gamlss_core, only : fit_gamlss_model
   use gamlss_model_selection_v02, only : stepwise_result_t
   implicit none
   private
   integer,parameter,public :: STEP_FORWARD=1,STEP_BACKWARD=2,STEP_BOTH=3
   public :: stepwise_gaic_mu
contains

   subroutine stepwise_gaic_mu(y,x_base,x_candidates,family,result,direction,start_selected,k, &
      x_sigma,x_nu,x_tau,weights,control,max_steps)
      real(dp),intent(in)::y(:),x_base(:,:),x_candidates(:,:)
      integer,intent(in)::family
      type(stepwise_result_t),intent(out)::result
      integer,intent(in),optional::direction,start_selected(:),max_steps
      real(dp),intent(in),optional::k,x_sigma(:,:),x_nu(:,:),x_tau(:,:),weights(:)
      type(gamlss_control_t),intent(in),optional::control
      logical,allocatable::selected(:),trial(:)
      integer,allocatable::sel(:)
      real(dp),allocatable::path(:),tmp_path(:),xm(:,:)
      type(gamlss_result_t)::fit,bestfit
      integer::dir,ncand,j,step,nstep,istat,bestj,bestaction
      real(dp)::kk,current,best,val

      if(size(x_base,1)/=size(y).or.size(x_candidates,1)/=size(y))then;result%status=1;return;end if
      dir=STEP_BOTH;if(present(direction))dir=direction
      if(dir<STEP_FORWARD.or.dir>STEP_BOTH)then;result%status=2;return;end if
      kk=2.0_dp;if(present(k))kk=k;ncand=size(x_candidates,2)
      allocate(selected(ncand));selected=.false.
      if(present(start_selected))then
         do j=1,size(start_selected)
            if(start_selected(j)<1.or.start_selected(j)>ncand)then;result%status=3;return;end if
            selected(start_selected(j))=.true.
         end do
      else if(dir==STEP_BACKWARD)then
         selected=.true.
      end if
      call logical_to_indices(selected,sel);call build_mu(x_base,x_candidates,sel,xm)
      call fit_model(xm,fit,istat);if(istat/=0)then;result%status=4;return;end if
      current=fit%global_deviance+kk*fit%df_fit;bestfit=fit
      allocate(path(1));path(1)=current
      nstep=max(1,2*max(1,ncand)+5);if(present(max_steps))nstep=max(1,max_steps)
      do step=1,nstep
         best=current;bestj=0;bestaction=0
         do j=1,ncand
            if((dir==STEP_FORWARD.or.dir==STEP_BOTH).and..not.selected(j))then
               trial=selected;trial(j)=.true.;call logical_to_indices(trial,sel)
               call build_mu(x_base,x_candidates,sel,xm);call fit_model(xm,fit,istat)
               if(istat==0)then
                  val=fit%global_deviance+kk*fit%df_fit
                  if(val<best-1.0e-8_dp)then;best=val;bestj=j;bestaction=1;bestfit=fit;end if
               end if
            end if
            if((dir==STEP_BACKWARD.or.dir==STEP_BOTH).and.selected(j))then
               trial=selected;trial(j)=.false.;call logical_to_indices(trial,sel)
               call build_mu(x_base,x_candidates,sel,xm);call fit_model(xm,fit,istat)
               if(istat==0)then
                  val=fit%global_deviance+kk*fit%df_fit
                  if(val<best-1.0e-8_dp)then;best=val;bestj=j;bestaction=-1;bestfit=fit;end if
               end if
            end if
         end do
         if(bestj==0)exit
         selected(bestj)=(bestaction==1);current=best
         call append_real(path,current,tmp_path);call move_alloc(tmp_path,path)
      end do
      call logical_to_indices(selected,sel);result%selected=sel;result%gaic_path=path
      result%final_model=bestfit;result%status=0
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
   end subroutine stepwise_gaic_mu

   subroutine logical_to_indices(mask,idx)
      logical,intent(in)::mask(:)
      integer,allocatable,intent(out)::idx(:)
      integer::j,k
      allocate(idx(count(mask)));k=0
      do j=1,size(mask);if(mask(j))then;k=k+1;idx(k)=j;end if;end do
   end subroutine logical_to_indices

   subroutine build_mu(base,candidates,sel,x)
      real(dp),intent(in)::base(:,:),candidates(:,:)
      integer,intent(in)::sel(:)
      real(dp),allocatable,intent(out)::x(:,:)
      integer::j,p
      p=size(base,2);allocate(x(size(base,1),p+size(sel)));x(:,1:p)=base
      do j=1,size(sel);x(:,p+j)=candidates(:,sel(j));end do
   end subroutine build_mu

   subroutine optional_design(n,xin,xout)
      integer,intent(in)::n
      real(dp),intent(in),optional::xin(:,:)
      real(dp),allocatable,intent(out)::xout(:,:)
      if(present(xin))then;xout=xin;else;allocate(xout(n,1));xout=1.0_dp;end if
   end subroutine optional_design

   subroutine append_real(a,v,b)
      real(dp),intent(in)::a(:),v
      real(dp),allocatable,intent(out)::b(:)
      allocate(b(size(a)+1));if(size(a)>0)b(1:size(a))=a;b(size(b))=v
   end subroutine append_real
end module gamlss_model_selection_v03
