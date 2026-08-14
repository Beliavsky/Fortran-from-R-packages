! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_models
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit, em_control
  use mclust_math, only : mixture_posterior, map_z
  use mclust_legacy_interfaces
  implicit none
  private
  public :: fit_model, mstep_model, initialize_responsibilities
  public :: model_supported, default_model_names

contains

  pure logical function model_supported(model,d) result(ok)
    character(len=*),intent(in)::model
    integer,intent(in)::d
    character(len=3)::m
    m=adjustl(model)
    if(d==1) then
      ok=trim(m)=='E' .or. trim(m)=='V'
    else
      select case(trim(m))
      case('EII','VII','EEI','VEI','EVI','VVI','EEE','EVE','VEE','VVE','EEV','VEV','EVV','VVV')
        ok=.true.
      case default
        ok=.false.
      end select
    end if
  end function model_supported

  subroutine default_model_names(d,names)
    integer,intent(in)::d
    character(len=3),allocatable,intent(out)::names(:)
    if(d==1) then
      allocate(names(2)); names=['E  ','V  ']
    else
      allocate(names(14))
      names=['EII','VII','EEI','VEI','EVI','VVI','EEE','EVE','VEE','VVE','EEV','VEV','EVV','VVV']
    end if
  end subroutine default_model_names

  subroutine fit_model(x,g,model,fit,control,z_init)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g
    character(len=*),intent(in)::model
    type(mclust_fit),intent(out)::fit
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::z_init(:,:)
    type(em_control)::ctl
    real(dp),allocatable::z(:,:),znew(:,:),mu(:,:),sigma(:,:,:),pro(:),ld(:)
    real(dp)::ll,prev,err
    integer::n,d,it,info

    n=size(x,1); d=size(x,2); ctl=em_control(); if(present(control)) ctl=control
    fit%n=n; fit%d=d; fit%g=g; fit%model_name=adjustl(model)
    if(n<1 .or. d<1 .or. g<1 .or. g>n .or. .not.model_supported(model,d)) then
      fit%status=-1; return
    end if
    allocate(z(n,g),znew(n,g),mu(d,g),sigma(d,d,g),pro(g),ld(n))
    if(present(z_init)) then
      if(size(z_init,1)/=n .or. size(z_init,2)/=g) then; fit%status=-2; return; end if
      z=z_init
      call normalize_rows(z,info); if(info/=0) then; fit%status=-3; return; end if
    else
      call initialize_responsibilities(x,g,z)
    end if

    prev=-huge(1.0_dp); err=huge(1.0_dp); info=0
    do it=1,ctl%max_iter
      call mstep_model(x,z,model,mu,sigma,pro,ctl,info)
      if(info/=0) exit
      if(ctl%equal_pro) pro=1.0_dp/real(g,dp)
      call mixture_posterior(x,pro,mu,sigma,znew,ld,info)
      if(info/=0) exit
      ll=sum(ld)
      if(it>1) then
        err=abs(ll-prev)/(1.0_dp+abs(ll))
        if(err<=ctl%tol) then
          z=znew; prev=ll; exit
        end if
      end if
      z=znew; prev=ll
    end do
    if(info/=0) then
      fit%status=info; return
    end if
    ! Recompute the final M-step and posterior so parameters and z correspond.
    call mstep_model(x,z,model,mu,sigma,pro,ctl,info)
    if(info/=0) then; fit%status=info; return; end if
    if(ctl%equal_pro) pro=1.0_dp/real(g,dp)
    call mixture_posterior(x,pro,mu,sigma,znew,ld,info)
    if(info/=0) then; fit%status=info; return; end if
    ll=sum(ld)

    fit%iterations=min(it,ctl%max_iter)
    fit%error=err
    fit%loglik=ll
    fit%status=merge(1,0,it>ctl%max_iter)
    allocate(fit%pro(g),fit%mean(d,g),fit%sigma(d,d,g),fit%z(n,g))
    allocate(fit%classification(n),fit%uncertainty(n))
    fit%pro=pro; fit%mean=mu; fit%sigma=sigma; fit%z=znew
    call map_z(fit%z,fit%classification,fit%uncertainty)
  end subroutine fit_model

  subroutine mstep_model(x,z,model,mu,sigma,pro,control,status)
    real(dp),intent(in)::x(:,:),z(:,:)
    character(len=*),intent(in)::model
    real(dp),intent(out)::mu(:,:),sigma(:,:,:),pro(:)
    type(em_control),intent(in)::control
    integer,intent(out)::status
    integer::n,p,g,k,j,lwork,maxi,info,niter
    real(dp)::sigsq,scale0,tol,errin
    real(dp),allocatable::sigs(:),scale(:),shape(:,:),shape1(:),o(:,:,:),ocom(:,:),u(:,:,:),c(:,:),w(:),s(:,:),scl(:),shp(:),ww(:,:)
    character(len=3)::m

    n=size(x,1); p=size(x,2); g=size(z,2); m=adjustl(model)
    status=0; mu=0.0_dp; sigma=0.0_dp; pro=0.0_dp
    if(size(mu,1)/=p .or. size(mu,2)/=g .or. size(sigma,1)/=p .or. &
       size(sigma,2)/=p .or. size(sigma,3)/=g .or. size(pro)/=g) then
      status=-1; return
    end if

    select case(trim(m))
    case('E')
      if(p/=1) then; status=-2; return; end if
      call ms1e(x(:,1),z,n,g,mu(1,:),sigsq,pro)
      do k=1,g; sigma(1,1,k)=sigsq; end do
    case('V')
      if(p/=1) then; status=-2; return; end if
      allocate(sigs(g)); sigs=0.0_dp
      call ms1v(x(:,1),z,n,g,mu(1,:),sigs,pro)
      do k=1,g; sigma(1,1,k)=sigs(k); end do
    case('EII')
      call mseii(x,z,n,p,g,mu,sigsq,pro)
      do k=1,g; do j=1,p; sigma(j,j,k)=sigsq; end do; end do
    case('VII')
      allocate(sigs(g)); sigs=0.0_dp
      call msvii(x,z,n,p,g,mu,sigs,pro)
      do k=1,g; do j=1,p; sigma(j,j,k)=sigs(k); end do; end do
    case('EEI')
      allocate(shape1(p)); shape1=0.0_dp; scale0=0.0_dp
      call mseei(x,z,n,p,g,mu,scale0,shape1,pro)
      do k=1,g; do j=1,p; sigma(j,j,k)=scale0*shape1(j); end do; end do
    case('EVI')
      allocate(shape(p,g)); shape=0.0_dp; scale0=0.0_dp
      call msevi(x,z,n,p,g,mu,scale0,shape,pro)
      do k=1,g; do j=1,p; sigma(j,j,k)=scale0*shape(j,k); end do; end do
    case('VEI')
      allocate(scale(g),shape1(p),scl(g),shp(p),ww(p,g)); scale=0.0_dp; shape1=0.0_dp
      scl=0.0_dp; shp=0.0_dp; ww=0.0_dp
      maxi=control%inner_max_iter; tol=control%inner_tol
      call msvei(x,z,n,p,g,maxi,tol,mu,scale,shape1,pro,scl,shp,ww)
      do k=1,g; do j=1,p; sigma(j,j,k)=scale(k)*shape1(j); end do; end do
    case('VVI')
      allocate(scale(g),shape(p,g)); scale=0.0_dp; shape=0.0_dp
      call msvvi(x,z,n,p,g,mu,scale,shape,pro)
      do k=1,g; do j=1,p; sigma(j,j,k)=scale(k)*shape(j,k); end do; end do
    case('EEE')
      allocate(w(p),ocom(p,p)); w=0.0_dp; ocom=0.0_dp
      call mseee(x,z,n,p,g,w,mu,ocom,pro)
      call upper_cross(ocom,sigma(:,:,1))
      do k=2,g; sigma(:,:,k)=sigma(:,:,1); end do
    case('EEV')
      lwork=max(3*min(n,p)+max(n,p),5*min(n,p),p+g)
      allocate(w(lwork),shape1(p),o(p,p,g)); w=0.0_dp; shape1=0.0_dp; o=0.0_dp; scale0=0.0_dp
      call mseev(x,z,n,p,g,w,lwork,mu,scale0,shape1,o,pro)
      do k=1,g; call eigen_cov_transposed(o(:,:,k),shape1,scale0,sigma(:,:,k)); end do
    case('VEV')
      lwork=max(3*min(n,p)+max(n,p),5*min(n,p),p+g)
      allocate(w(lwork),scale(g),shape1(p),o(p,p,g)); w=0.0_dp; scale=0.0_dp; shape1=0.0_dp; o=0.0_dp
      maxi=control%inner_max_iter; tol=control%inner_tol
      call msvev(x,z,n,p,g,w,lwork,maxi,tol,mu,scale,shape1,o,pro)
      do k=1,g; call eigen_cov_transposed(o(:,:,k),shape1,scale(k),sigma(:,:,k)); end do
    case('EVV')
      lwork=max(3*min(n,p)+max(n,p),5*min(n,p),g)
      allocate(scale(g),shape(p,g),o(p,p,g),u(p,p,g)); scale=0.0_dp; shape=0.0_dp; o=0.0_dp; u=0.0_dp
      info=0
      call msevv(x,z,n,p,g,mu,o,u,scale,shape,pro,lwork,info,control%eps)
      if(info/=0) then; status=100+info; return; end if
      do k=1,g; call eigen_cov_transposed(o(:,:,k),shape(:,k),scale(1),sigma(:,:,k)); end do
    case('VEE')
      lwork=max(3*min(n,p)+max(n,p),5*min(n,p),p+g)
      allocate(scale(g),u(p,p,g),c(p,p)); scale=1.0_dp; u=0.0_dp; c=0.0_dp
      info=0; niter=0; errin=huge(1.0_dp)
      call msvee(x,z,n,p,g,mu,u,c,scale,pro,lwork,info,control%inner_max_iter,control%inner_tol,niter,errin)
      if(info/=0) then; status=100+info; return; end if
      do k=1,g; sigma(:,:,k)=scale(k)*c; end do
    case('EVE')
      lwork=max(3*min(n,p)+max(n,p),5*min(n,p),p+g)
      allocate(shape(p,g),u(p,p,g),ocom(p,p)); shape=1.0_dp; u=0.0_dp; ocom=0.0_dp; scale0=1.0_dp
      do j=1,p; ocom(j,j)=1.0_dp; end do
      info=0; niter=0; errin=huge(1.0_dp)
      call mseve(x,z,n,p,g,mu,u,ocom,scale0,shape,pro,lwork,info,control%inner_max_iter,control%inner_tol,niter,errin,control%eps)
      if(info/=0) then; status=100+info; return; end if
      do k=1,g; call eigen_cov_transposed(ocom,shape(:,k),scale0,sigma(:,:,k)); end do
    case('VVE')
      lwork=max(3*min(n,p)+max(n,p),5*min(n,p),p+g)
      allocate(scale(g),shape(p,g),u(p,p,g),ocom(p,p)); scale=1.0_dp; shape=1.0_dp; u=0.0_dp; ocom=0.0_dp
      do j=1,p; ocom(j,j)=1.0_dp; end do
      info=0; niter=0; errin=huge(1.0_dp)
      call msvve(x,z,n,p,g,mu,u,ocom,scale,shape,pro,lwork,info,control%inner_max_iter,control%inner_tol,niter,errin,control%eps)
      if(info/=0) then; status=100+info; return; end if
      do k=1,g; call eigen_cov_transposed(ocom,shape(:,k),scale(k),sigma(:,:,k)); end do
    case('VVV')
      allocate(w(p),u(p,p,g),s(p,p)); w=0.0_dp; u=0.0_dp; s=0.0_dp
      call msvvv(x,z,n,p,g,w,mu,u,pro,s)
      do k=1,g; call upper_cross(u(:,:,k),sigma(:,:,k)); end do
    case default
      status=-3; return
    end select

    if(any(.not.(sigma<huge(1.0_dp))) .or. any(.not.(pro<huge(1.0_dp))) .or. &
       any(pro<0.0_dp)) then
      status=200; return
    end if
    do k=1,g
      if(any([(sigma(j,j,k)<=max(control%eps,0.0_dp),j=1,p)])) then; status=201; return; end if
    end do
  end subroutine mstep_model

  subroutine initialize_responsibilities(x,g,z)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g
    real(dp),intent(out)::z(:,:)
    real(dp),allocatable::cent(:,:),d2(:),best(:)
    integer,allocatable::cl(:)
    integer::n,p,i,k,j,idx,it
    real(dp)::dist,old
    n=size(x,1); p=size(x,2)
    allocate(cent(g,p),d2(g),best(n),cl(n))
    cent(1,:)=sum(x,dim=1)/real(n,dp)
    best=huge(1.0_dp)
    do k=2,g
      do i=1,n
        dist=sum((x(i,:)-cent(k-1,:))**2); best(i)=min(best(i),dist)
      end do
      idx=maxloc(best,dim=1); cent(k,:)=x(idx,:)
    end do
    cl=1
    do it=1,30
      old=0.0_dp
      do i=1,n
        d2=[(sum((x(i,:)-cent(k,:))**2),k=1,g)]
        j=minloc(d2,dim=1); if(j/=cl(i)) old=old+1.0_dp; cl(i)=j
      end do
      do k=1,g
        if(count(cl==k)>0) then
          cent(k,:)=0.0_dp
          do i=1,n; if(cl(i)==k) cent(k,:)=cent(k,:)+x(i,:); end do
          cent(k,:)=cent(k,:)/real(count(cl==k),dp)
        end if
      end do
      if(old<0.5_dp) exit
    end do
    if(g==1) then
      z=1.0_dp
    else
      z=(0.05_dp/real(g-1,dp))
      do i=1,n; z(i,cl(i))=0.95_dp; end do
    end if
  end subroutine initialize_responsibilities

  subroutine normalize_rows(z,status)
    real(dp),intent(inout)::z(:,:)
    integer,intent(out)::status
    real(dp)::s
    integer::i
    status=0
    if(any(z<0.0_dp)) then; status=1; return; end if
    do i=1,size(z,1)
      s=sum(z(i,:)); if(s<=0.0_dp) then; status=2; return; end if; z(i,:)=z(i,:)/s
    end do
  end subroutine normalize_rows

  subroutine upper_cross(u,sigma)
    real(dp),intent(in)::u(:,:)
    real(dp),intent(out)::sigma(:,:)
    sigma=matmul(transpose(u),u)
  end subroutine upper_cross

  subroutine eigen_cov_transposed(ot,shape,scale,sigma)
    real(dp),intent(in)::ot(:,:),shape(:),scale
    real(dp),intent(out)::sigma(:,:)
    real(dp),allocatable::tmp(:,:)
    integer::j
    allocate(tmp(size(ot,1),size(ot,2))); tmp=transpose(ot)
    do j=1,size(shape); tmp(:,j)=tmp(:,j)*(scale*shape(j)); end do
    sigma=matmul(tmp,ot)
    sigma=0.5_dp*(sigma+transpose(sigma))
  end subroutine eigen_cov_transposed

end module mclust_models
