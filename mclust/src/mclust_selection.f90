! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_selection_mod
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit, mclust_selection, em_control
  use mclust_models, only : fit_model, default_model_names, model_supported
  use mclust_hierarchical, only : hc_responsibilities
  implicit none
  private
  public :: n_var_params, n_mclust_params, bic_value, icl_value
  public :: mclust_bic, mclust_select, pick_bic

contains

  integer function n_var_params(model,d,g) result(npar)
    character(len=*),intent(in)::model
    integer,intent(in)::d,g
    select case(trim(adjustl(model)))
    case('E','EII'); npar=1
    case('V','VII'); npar=g
    case('EEI'); npar=d
    case('VEI'); npar=g+d-1
    case('EVI'); npar=1+g*(d-1)
    case('VVI'); npar=g*d
    case('EEE'); npar=d*(d+1)/2
    case('EVE'); npar=1+g*(d-1)+d*(d-1)/2
    case('VEE'); npar=g+(d-1)+d*(d-1)/2
    case('VVE'); npar=g+g*(d-1)+d*(d-1)/2
    case('EEV'); npar=1+(d-1)+g*d*(d-1)/2
    case('VEV'); npar=g+(d-1)+g*d*(d-1)/2
    case('EVV'); npar=1-g+g*d*(d+1)/2
    case('VVV'); npar=g*d*(d+1)/2
    case default; npar=-1
    end select
  end function n_var_params

  integer function n_mclust_params(model,d,g,equal_pro,noise) result(npar)
    character(len=*),intent(in)::model
    integer,intent(in)::d,g
    logical,intent(in),optional::equal_pro,noise
    logical::eq,nz
    eq=.false.; if(present(equal_pro)) eq=equal_pro
    nz=.false.; if(present(noise)) nz=noise
    if(g==0) then
      npar=merge(1,-1,nz); return
    end if
    npar=n_var_params(model,d,g)
    if(npar<0) return
    npar=npar+g*d
    if(.not.eq) npar=npar+g-1
    if(nz) npar=npar+2
  end function n_mclust_params

  pure real(dp) function bic_value(loglik,n,nparams) result(bic)
    real(dp),intent(in)::loglik
    integer,intent(in)::n,nparams
    bic=2.0_dp*loglik-real(nparams,dp)*log(real(n,dp))
  end function bic_value

  real(dp) function icl_value(fit) result(icl)
    type(mclust_fit),intent(in)::fit
    integer::i,k
    real(dp)::s
    if(.not.allocated(fit%z)) then; icl=fit%bic; return; end if
    s=0.0_dp
    do i=1,fit%n
      k=maxloc(fit%z(i,:),dim=1)
      if(fit%z(i,k)>0.0_dp) s=s+log(fit%z(i,k))
    end do
    icl=fit%bic+2.0_dp*s
  end function icl_value

  subroutine mclust_bic(x,g_values,selection,model_names,control)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g_values(:)
    type(mclust_selection),intent(out)::selection
    character(len=*),intent(in),optional::model_names(:)
    type(em_control),intent(in),optional::control
    character(len=3),allocatable::mods(:)
    type(em_control)::ctl
    type(mclust_fit)::fit
    integer::i,j,np,d,n
    real(dp)::bestbic
    real(dp),allocatable::zinit(:,:)
    integer::hinfo

    n=size(x,1); d=size(x,2); ctl=em_control(); if(present(control)) ctl=control
    if(present(model_names)) then
      allocate(mods(size(model_names)))
      do j=1,size(model_names); mods(j)=adjustl(model_names(j)); end do
    else
      call default_model_names(d,mods)
    end if
    selection%n=n; selection%d=d; selection%n_g=size(g_values); selection%n_models=size(mods)
    allocate(selection%g_values(size(g_values)),selection%model_names(size(mods)))
    allocate(selection%bic(size(g_values),size(mods)),selection%status(size(g_values),size(mods)))
    selection%g_values=g_values; selection%model_names=mods
    selection%bic=-huge(1.0_dp); selection%status=-99
    bestbic=-huge(1.0_dp)
    do i=1,size(g_values)
      if(g_values(i)<1 .or. g_values(i)>n) cycle
      if(allocated(zinit)) deallocate(zinit)
      hinfo=-1
      if(ctl%use_hc .and. g_values(i)>1) then
        call hc_responsibilities(x,g_values(i),zinit,ctl%hc_model,hinfo,ctl%hc_use)
      end if
      do j=1,size(mods)
        if(.not.model_supported(mods(j),d)) cycle
        if(hinfo==0) then
          call fit_model(x,g_values(i),mods(j),fit,ctl,zinit)
        else
          call fit_model(x,g_values(i),mods(j),fit,ctl)
        end if
        selection%status(i,j)=fit%status
        if(fit%status<0 .or. .not.allocated(fit%z)) cycle
        np=n_mclust_params(mods(j),d,g_values(i),ctl%equal_pro)
        fit%bic=bic_value(fit%loglik,n,np)
        fit%icl=icl_value(fit)
        selection%bic(i,j)=fit%bic
        if(fit%bic>bestbic) then
          bestbic=fit%bic
          if(allocated(selection%best)) deallocate(selection%best)
          allocate(selection%best,source=fit)
        end if
      end do
    end do
  end subroutine mclust_bic

  subroutine mclust_select(x,fit,g_values,model_names,control,status)
    real(dp),intent(in)::x(:,:)
    type(mclust_fit),intent(out)::fit
    integer,intent(in),optional::g_values(:)
    character(len=*),intent(in),optional::model_names(:)
    type(em_control),intent(in),optional::control
    integer,intent(out),optional::status
    type(mclust_selection)::sel
    integer,allocatable::gs(:)
    integer::i,gmax
    gmax=min(9,max(1,size(x,1)))
    if(present(g_values)) then
      allocate(gs(size(g_values))); gs=g_values
    else
      allocate(gs(gmax)); gs=[(i,i=1,gmax)]
    end if
    if(present(model_names)) then
      call mclust_bic(x,gs,sel,model_names,control)
    else
      call mclust_bic(x,gs,sel,control=control)
    end if
    if(allocated(sel%best)) then
      fit=sel%best
      if(present(status)) status=0
    else
      fit%status=-1
      if(present(status)) status=-1
    end if
  end subroutine mclust_select

  subroutine pick_bic(selection,n_pick,g_out,model_out,bic_out)
    type(mclust_selection),intent(in)::selection
    integer,intent(in)::n_pick
    integer,allocatable,intent(out)::g_out(:)
    character(len=3),allocatable,intent(out)::model_out(:)
    real(dp),allocatable,intent(out)::bic_out(:)
    logical,allocatable::used(:,:)
    integer::q,i,j,ii,jj,np
    real(dp)::best
    np=min(n_pick,selection%n_g*selection%n_models)
    allocate(g_out(np),model_out(np),bic_out(np),used(selection%n_g,selection%n_models)); used=.false.
    do q=1,np
      best=-huge(1.0_dp); ii=0; jj=0
      do i=1,selection%n_g; do j=1,selection%n_models
        if(.not.used(i,j) .and. selection%bic(i,j)>best) then; best=selection%bic(i,j); ii=i; jj=j; end if
      end do; end do
      if(ii==0) then
        g_out(q:)=0; model_out(q:)=''; bic_out(q:)=-huge(1.0_dp); exit
      end if
      used(ii,jj)=.true.; g_out(q)=selection%g_values(ii); model_out(q)=selection%model_names(jj); bic_out(q)=best
    end do
  end subroutine pick_bic

end module mclust_selection_mod
