! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! Uses the supplied modern Fortran translation of R's splines package.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_smoothing
   use vgam_kinds, only : dp
   use splines, only : bs_basis, natural_spline_basis
   use vgam_vglm, only : vglm_result_t, fit_vglm, family_gaussian, &
      family_poisson, family_binomial
   implicit none
   private
   integer,parameter,public::basis_bspline=1,basis_natural=2
   type,public::vgam_smooth_result_t
      type(vglm_result_t)::fit
      real(dp),allocatable::interior_knots(:)
      real(dp)::boundary_knots(2)=0.0_dp
      integer::degree=3
      integer::basis_kind=basis_bspline
      integer::n_linear=0
      real(dp)::lambda=0.0_dp
   contains
      procedure::predict=>predict_smooth
   end type vgam_smooth_result_t
   public::fit_pspline_vglm,fit_gam_gaussian,fit_gam_poisson,fit_gam_binomial
   public::difference_penalty
contains

   subroutine difference_penalty(ncoef,order,lambda,p)
      integer,intent(in)::ncoef,order
      real(dp),intent(in)::lambda
      real(dp),allocatable,intent(out)::p(:,:)
      real(dp),allocatable::d(:,:),next(:,:)
      integer::i,j,nr
      allocate(d(ncoef,ncoef));d=0.0_dp
      do i=1,ncoef;d(i,i)=1.0_dp;end do
      nr=ncoef
      do j=1,order
         if(nr<=1)exit
         allocate(next(nr-1,ncoef))
         do i=1,nr-1;next(i,:)=d(i+1,:)-d(i,:);end do
         call move_alloc(next,d);nr=nr-1
      end do
      p=lambda*matmul(transpose(d),d)
   end subroutine difference_penalty

   subroutine fit_pspline_vglm(x_smooth,y,family,result,df,lambda,x_linear, &
                               link_id,degree,basis_kind,penalty_order)
      real(dp),intent(in)::x_smooth(:),y(:)
      integer,intent(in)::family
      type(vgam_smooth_result_t),intent(out)::result
      integer,intent(in),optional::df,link_id,degree,basis_kind,penalty_order
      real(dp),intent(in),optional::lambda,x_linear(:,:)
      real(dp),allocatable::b(:,:),design(:,:),sp(:,:),pen(:,:),inner(:)
      real(dp)::lam
      integer::n,nb,nl,dff,deg,bkind,pord,stat,j0
      n=size(y)
      if(size(x_smooth)/=n)then;result%fit%status=1;return;end if
      dff=6;if(present(df))dff=df
      lam=1.0_dp;if(present(lambda))lam=lambda
      deg=3;if(present(degree))deg=degree
      bkind=basis_bspline;if(present(basis_kind))bkind=basis_kind
      pord=2;if(present(penalty_order))pord=penalty_order
      result%boundary_knots=[minval(x_smooth),maxval(x_smooth)]
      if(bkind==basis_natural)then
         call natural_spline_basis(x_smooth,b,df=dff,intercept=.false., &
            boundary_knots=result%boundary_knots,interior_knots=inner,status=stat)
      else
         call bs_basis(x_smooth,b,degree=deg,df=dff,intercept=.false., &
            boundary_knots=result%boundary_knots,interior_knots=inner,status=stat)
      end if
      if(stat/=0)then;result%fit%status=10+stat;return;end if
      nb=size(b,2);nl=0
      if(present(x_linear))then
         if(size(x_linear,1)/=n)then;result%fit%status=2;return;end if
         nl=size(x_linear,2)
      end if
      allocate(design(n,1+nl+nb));design(:,1)=1.0_dp
      if(nl>0)design(:,2:1+nl)=x_linear
      design(:,2+nl:)=b
      call difference_penalty(nb,pord,lam,sp)
      allocate(pen(size(design,2),size(design,2)));pen=0.0_dp
      j0=2+nl;pen(j0:,j0:)=sp
      if(present(link_id))then
         call fit_vglm(y,design,family,result%fit,link_id=link_id,penalty=pen)
      else
         call fit_vglm(y,design,family,result%fit,penalty=pen)
      end if
      result%interior_knots=inner
      result%degree=deg;result%basis_kind=bkind;result%n_linear=nl;result%lambda=lam
   end subroutine fit_pspline_vglm

   subroutine fit_gam_gaussian(x_smooth,y,result,df,lambda,x_linear)
      real(dp),intent(in)::x_smooth(:),y(:)
      type(vgam_smooth_result_t),intent(out)::result
      integer,intent(in),optional::df
      real(dp),intent(in),optional::lambda,x_linear(:,:)
      call fit_pspline_vglm(x_smooth,y,family_gaussian,result,df,lambda,x_linear)
   end subroutine fit_gam_gaussian

   subroutine fit_gam_poisson(x_smooth,y,result,df,lambda,x_linear)
      real(dp),intent(in)::x_smooth(:),y(:)
      type(vgam_smooth_result_t),intent(out)::result
      integer,intent(in),optional::df
      real(dp),intent(in),optional::lambda,x_linear(:,:)
      call fit_pspline_vglm(x_smooth,y,family_poisson,result,df,lambda,x_linear)
   end subroutine fit_gam_poisson

   subroutine fit_gam_binomial(x_smooth,y,result,df,lambda,x_linear)
      real(dp),intent(in)::x_smooth(:),y(:)
      type(vgam_smooth_result_t),intent(out)::result
      integer,intent(in),optional::df
      real(dp),intent(in),optional::lambda,x_linear(:,:)
      call fit_pspline_vglm(x_smooth,y,family_binomial,result,df,lambda,x_linear)
   end subroutine fit_gam_binomial

   function predict_smooth(self,x_smooth,x_linear,response) result(pred)
      class(vgam_smooth_result_t),intent(in)::self
      real(dp),intent(in)::x_smooth(:)
      real(dp),intent(in),optional::x_linear(:,:)
      logical,intent(in),optional::response
      real(dp),allocatable::pred(:),b(:,:),design(:,:)
      integer::n,nb,nl,stat
      n=size(x_smooth);nl=self%n_linear
      if(self%basis_kind==basis_natural)then
         call natural_spline_basis(x_smooth,b,knots=self%interior_knots, &
            intercept=.false.,boundary_knots=self%boundary_knots,status=stat)
      else
         call bs_basis(x_smooth,b,degree=self%degree,knots=self%interior_knots, &
            intercept=.false.,boundary_knots=self%boundary_knots,status=stat)
      end if
      if(stat/=0)then;allocate(pred(0));return;end if
      nb=size(b,2);allocate(design(n,1+nl+nb));design(:,1)=1.0_dp
      if(nl>0)then
         if(.not.present(x_linear))then;allocate(pred(0));return;end if
         design(:,2:1+nl)=x_linear
      end if
      design(:,2+nl:)=b
      if(present(response))then
         pred=self%fit%predict(design,response=response)
      else
         pred=self%fit%predict(design)
      end if
   end function predict_smooth
end module vgam_smoothing
