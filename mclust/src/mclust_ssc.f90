! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_ssc
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit, em_control
  use mclust_models, only : mstep_model, initialize_responsibilities, default_model_names, model_supported
  use mclust_math, only : mixture_posterior, dmvnorm, map_z
  use mclust_selection_mod, only : n_mclust_params, bic_value, icl_value
  implicit none
  private
  public :: fit_model_ssc, mclust_ssc_select
contains

  subroutine fit_model_ssc(x,class,g,model,fit,control,status)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::class(:) ! 0 means unknown; positive integers are known classes
    integer,intent(in)::g
    character(len=*),intent(in)::model
    type(mclust_fit),intent(out)::fit
    type(em_control),intent(in),optional::control
    integer,intent(out),optional::status
    type(em_control)::ctl
    real(dp),allocatable::z(:,:),znew(:,:),mu(:,:),sigma(:,:,:),pro(:),ld(:),zu(:,:)
    integer,allocatable::known(:),unknown(:)
    real(dp)::ll,prev,err
    integer::n,d,nclass,it,info,i,k,nk,nu

    n=size(x,1); d=size(x,2); ctl=em_control(); if(present(control))ctl=control
    fit%n=n; fit%d=d; fit%g=g; fit%model_name=adjustl(model)
    if(size(class)/=n .or. any(class<0) .or. n<1 .or. g<1 .or. .not.model_supported(model,d)) then
      fit%status=-1; if(present(status))status=-1; return
    end if
    nclass=maxval(class)
    if(nclass<1 .or. g<nclass .or. .not.any(class==0)) then
      fit%status=-2; if(present(status))status=-2; return
    end if
    nk=count(class>0); nu=n-nk
    allocate(known(nk),unknown(nu)); known=pack([(i,i=1,n)],class>0); unknown=pack([(i,i=1,n)],class==0)
    allocate(z(n,g),znew(n,g),mu(d,g),sigma(d,d,g),pro(g),ld(n)); z=0.0_dp
    do i=1,nk
      z(known(i),class(known(i)))=1.0_dp
    end do
    allocate(zu(nu,g)); call initialize_responsibilities(x(unknown,:),g,zu)
    ! Upstream initializes every unlabeled row with the k-means component proportions.
    do k=1,g; z(unknown,k)=sum(zu(:,k))/real(nu,dp); end do
    do i=1,nu
      if(sum(z(unknown(i),:))<=0.0_dp) z(unknown(i),:)=1.0_dp/real(g,dp)
      z(unknown(i),:)=z(unknown(i),:)/sum(z(unknown(i),:))
    end do

    prev=-huge(1.0_dp); err=huge(1.0_dp); info=0
    do it=1,ctl%max_iter
      call mstep_model(x,z,model,mu,sigma,pro,ctl,info); if(info/=0)exit
      if(ctl%equal_pro)pro=1.0_dp/real(g,dp)
      call mixture_posterior(x,pro,mu,sigma,znew,ld,info); if(info/=0)exit
      do i=1,nk
        znew(known(i),:)=0.0_dp; znew(known(i),class(known(i)))=1.0_dp
      end do
      call ssc_loglik(x,class,pro,mu,sigma,ld,ll,info); if(info/=0)exit
      if(it>1) then
        err=abs(ll-prev)/(1.0_dp+abs(ll))
        if(err<=ctl%tol) then; z=znew; prev=ll; exit; end if
      end if
      z=znew; prev=ll
    end do
    if(info/=0) then; fit%status=info; if(present(status))status=info; return; end if
    call mstep_model(x,z,model,mu,sigma,pro,ctl,info)
    if(info/=0)then
      fit%status=info; if(present(status))status=info; return
    end if
    if(ctl%equal_pro)pro=1.0_dp/real(g,dp)
    call mixture_posterior(x,pro,mu,sigma,znew,ld,info)
    if(info/=0)then
      fit%status=info; if(present(status))status=info; return
    end if
    do i=1,nk
      znew(known(i),:)=0.0_dp; znew(known(i),class(known(i)))=1.0_dp
    end do
    call ssc_loglik(x,class,pro,mu,sigma,ld,ll,info)
    if(info/=0)then
      fit%status=info; if(present(status))status=info; return
    end if
    fit%iterations=min(it,ctl%max_iter); fit%error=err; fit%loglik=ll; fit%status=merge(1,0,it>ctl%max_iter)
    fit%bic=bic_value(ll,n,n_mclust_params(model,d,g,ctl%equal_pro));
    allocate(fit%pro(g),fit%mean(d,g),fit%sigma(d,d,g),fit%z(n,g),fit%classification(n),fit%uncertainty(n))
    fit%pro=pro; fit%mean=mu; fit%sigma=sigma; fit%z=znew
    call map_z(fit%z,fit%classification,fit%uncertainty); fit%icl=icl_value(fit)
    if(present(status))status=fit%status
  end subroutine fit_model_ssc

  subroutine ssc_loglik(x,class,pro,mu,sigma,mixture_ld,ll,status)
    real(dp),intent(in)::x(:,:),pro(:),mu(:,:),sigma(:,:,:),mixture_ld(:)
    integer,intent(in)::class(:)
    real(dp),intent(out)::ll
    integer,intent(out)::status
    real(dp),allocatable::ld(:)
    integer::k,info
    ll=sum(pack(mixture_ld,class==0)); status=0
    allocate(ld(size(x,1)))
    do k=1,size(pro)
      if(.not.any(class==k)) cycle
      call dmvnorm(x,mu(:,k),sigma(:,:,k),ld,info)
      if(info/=0) then; status=10+k; return; end if
      ll=ll+sum(pack(ld+log(max(pro(k),tiny(1.0_dp))),class==k))
    end do
  end subroutine ssc_loglik

  subroutine mclust_ssc_select(x,class,g,fit,model_names,control,status)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::class(:)
    integer,intent(in)::g
    type(mclust_fit),intent(out)::fit
    character(len=*),intent(in),optional::model_names(:)
    type(em_control),intent(in),optional::control
    integer,intent(out),optional::status
    character(len=3),allocatable::mods(:)
    type(mclust_fit)::f
    real(dp)::best
    integer::j,info
    if(present(model_names)) then
      allocate(mods(size(model_names)))
      do j=1,size(mods); mods(j)=adjustl(model_names(j)); end do
    else
      call default_model_names(size(x,2),mods)
    end if
    best=-huge(1.0_dp); info=-1
    do j=1,size(mods)
      if(.not.model_supported(mods(j),size(x,2)))cycle
      call fit_model_ssc(x,class,g,mods(j),f,control)
      if(f%status<0)cycle
      if(f%bic>best) then; best=f%bic; fit=f; info=0; end if
    end do
    if(present(status))status=info
  end subroutine mclust_ssc_select
end module mclust_ssc
