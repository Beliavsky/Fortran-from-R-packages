! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_classification
  use mclust_kinds, only : dp
  use mclust_types, only : mclust_fit, em_control
  use mclust_selection_mod, only : mclust_select
  use mclust_math, only : mixture_posterior, logsumexp
  implicit none
  private

  type, public :: mclust_da_fit
    integer :: n = 0
    integer :: d = 0
    integer :: n_classes = 0
    integer, allocatable :: labels(:)
    real(dp), allocatable :: prior(:)
    type(mclust_fit), allocatable :: model(:)
  end type mclust_da_fit

  public :: fit_mclust_da, predict_mclust_da, class_error_rate

contains

  subroutine fit_mclust_da(x,class,da,g_values,model_names,control,prior,status)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::class(:)
    type(mclust_da_fit),intent(out)::da
    integer,intent(in),optional::g_values(:)
    character(len=*),intent(in),optional::model_names(:)
    type(em_control),intent(in),optional::control
    real(dp),intent(in),optional::prior(:)
    integer,intent(out),optional::status
    integer,allocatable::labels(:),idx(:)
    real(dp),allocatable::xc(:,:)
    integer::c,i,j,nc,info

    if(size(x,1)/=size(class)) then; if(present(status)) status=-1; return; end if
    call unique_int(class,labels)
    da%n=size(x,1); da%d=size(x,2); da%n_classes=size(labels)
    allocate(da%labels(size(labels)),da%prior(size(labels)),da%model(size(labels)))
    da%labels=labels
    if(present(prior)) then
      if(size(prior)/=size(labels) .or. any(prior<0.0_dp) .or. sum(prior)<=0.0_dp) then
        if(present(status)) status=-2; return
      end if
      da%prior=prior/sum(prior)
    else
      do c=1,size(labels); da%prior(c)=real(count(class==labels(c)),dp)/real(size(class),dp); end do
    end if
    info=0
    do c=1,size(labels)
      nc=count(class==labels(c)); allocate(idx(nc),xc(nc,size(x,2)))
      j=0; do i=1,size(class); if(class(i)==labels(c)) then; j=j+1; idx(j)=i; xc(j,:)=x(i,:); end if; end do
      if(present(g_values) .and. present(model_names)) then
        call mclust_select(xc,da%model(c),g_values,model_names,control,info)
      else if(present(g_values)) then
        call mclust_select(xc,da%model(c),g_values=g_values,control=control,status=info)
      else if(present(model_names)) then
        call mclust_select(xc,da%model(c),model_names=model_names,control=control,status=info)
      else
        call mclust_select(xc,da%model(c),control=control,status=info)
      end if
      deallocate(idx,xc)
      if(info/=0) exit
    end do
    if(present(status)) status=info
  end subroutine fit_mclust_da

  subroutine predict_mclust_da(da,x,posterior,classification,log_density,status)
    type(mclust_da_fit),intent(in)::da
    real(dp),intent(in)::x(:,:)
    real(dp),allocatable,intent(out)::posterior(:,:)
    integer,allocatable,intent(out)::classification(:)
    real(dp),allocatable,intent(out),optional::log_density(:)
    integer,intent(out),optional::status
    real(dp),allocatable::z(:,:),ld(:),ldall(:,:),row(:)
    integer::c,i,info,k

    allocate(posterior(size(x,1),da%n_classes),classification(size(x,1)))
    allocate(ldall(size(x,1),da%n_classes),row(da%n_classes))
    do c=1,da%n_classes
      allocate(z(size(x,1),da%model(c)%g),ld(size(x,1)))
      call mixture_posterior(x,da%model(c)%pro,da%model(c)%mean,da%model(c)%sigma,z,ld,info)
      if(info/=0) then; if(present(status)) status=10+c; return; end if
      ldall(:,c)=ld+log(max(da%prior(c),tiny(1.0_dp)))
      deallocate(z,ld)
    end do
    if(present(log_density)) allocate(log_density(size(x,1)))
    do i=1,size(x,1)
      row=ldall(i,:); k=maxloc(row,dim=1); classification(i)=da%labels(k)
      if(present(log_density)) log_density(i)=logsumexp(row)
      posterior(i,:)=exp(row-logsumexp(row))
    end do
    if(present(status)) status=0
  end subroutine predict_mclust_da

  pure real(dp) function class_error_rate(pred,true_class) result(rate)
    integer,intent(in)::pred(:),true_class(:)
    if(size(pred)/=size(true_class) .or. size(pred)==0) then; rate=1.0_dp; return; end if
    rate=real(count(pred/=true_class),dp)/real(size(pred),dp)
  end function class_error_rate

  subroutine unique_int(x,u)
    integer,intent(in)::x(:)
    integer,allocatable,intent(out)::u(:)
    integer,allocatable::tmp(:)
    integer::i,n
    allocate(tmp(size(x))); n=0
    do i=1,size(x)
      if(n==0 .or. all(tmp(1:n)/=x(i))) then; n=n+1; tmp(n)=x(i); end if
    end do
    allocate(u(n)); u=tmp(1:n)
  end subroutine unique_int
end module mclust_classification
