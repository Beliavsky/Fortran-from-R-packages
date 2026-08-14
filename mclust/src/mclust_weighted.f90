! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_weighted
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit, em_control
  use mclust_models, only : mstep_model, initialize_responsibilities, model_supported
  use mclust_math, only : mixture_posterior, map_z
  use mclust_selection_mod, only : n_mclust_params, bic_value, icl_value
  implicit none
  private
  public :: fit_model_weighted
contains
  subroutine fit_model_weighted(x,g,model,weights,fit,control,z_init)
    real(dp),intent(in)::x(:,:),weights(:)
    integer,intent(in)::g
    character(len=*),intent(in)::model
    type(mclust_fit),intent(out)::fit
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::z_init(:,:)
    type(em_control)::ctl
    real(dp),allocatable::z(:,:),zw(:,:),znew(:,:),mu(:,:),sigma(:,:,:),pro(:),ld(:),w(:)
    real(dp)::ll,prev,err,mw
    integer::n,d,it,info,np
    n=size(x,1); d=size(x,2); ctl=em_control(); if(present(control)) ctl=control
    fit%n=n; fit%d=d; fit%g=g; fit%model_name=adjustl(model)
    if(size(weights)/=n .or. any(weights<0.0_dp) .or. maxval(weights)<=0.0_dp .or. .not.model_supported(model,d)) then
      fit%status=-1; return
    end if
    allocate(w(n)); w=weights/maxval(weights); mw=sum(w)/real(n,dp)
    allocate(z(n,g),zw(n,g),znew(n,g),mu(d,g),sigma(d,d,g),pro(g),ld(n))
    if(present(z_init)) then
      z=z_init
      call normalize_rows(z,info); if(info/=0) then; fit%status=-2; return; end if
    else
      call initialize_responsibilities(x,g,z)
    end if
    prev=-huge(1.0_dp); err=huge(1.0_dp); info=0
    do it=1,ctl%max_iter
      zw=z*spread(w,2,g)
      call mstep_model(x,zw,model,mu,sigma,pro,ctl,info); if(info/=0) exit
      pro=pro/mw; pro=pro/sum(pro); if(ctl%equal_pro)pro=1.0_dp/real(g,dp)
      call mixture_posterior(x,pro,mu,sigma,znew,ld,info); if(info/=0) exit
      ll=sum(w*ld)/mw
      if(it>1) then
        err=abs(ll-prev)/(1.0_dp+abs(ll)); if(err<=ctl%tol) then; z=znew; prev=ll; exit; end if
      end if
      z=znew; prev=ll
    end do
    if(info/=0) then; fit%status=info; return; end if
    zw=z*spread(w,2,g); call mstep_model(x,zw,model,mu,sigma,pro,ctl,info)
    pro=pro/mw; pro=pro/sum(pro); if(ctl%equal_pro)pro=1.0_dp/real(g,dp)
    call mixture_posterior(x,pro,mu,sigma,znew,ld,info); ll=sum(w*ld)/mw
    fit%iterations=min(it,ctl%max_iter); fit%error=err; fit%loglik=ll; fit%status=0
    allocate(fit%pro(g),fit%mean(d,g),fit%sigma(d,d,g),fit%z(n,g),fit%classification(n),fit%uncertainty(n))
    fit%pro=pro; fit%mean=mu; fit%sigma=sigma; fit%z=znew; call map_z(znew,fit%classification,fit%uncertainty)
    np=n_mclust_params(model,d,g,ctl%equal_pro); fit%bic=bic_value(ll,n,np); fit%icl=icl_value(fit)
  end subroutine fit_model_weighted
  subroutine normalize_rows(z,status)
    real(dp),intent(inout)::z(:,:); integer,intent(out)::status; integer::i; real(dp)::s
    status=0; do i=1,size(z,1); s=sum(z(i,:)); if(s<=0.0_dp) then; status=1; return; end if; z(i,:)=z(i,:)/s; end do
  end subroutine
end module mclust_weighted
