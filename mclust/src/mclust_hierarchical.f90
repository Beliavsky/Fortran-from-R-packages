! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_hierarchical
  use mclust_kinds, only : dp
  use mclust_legacy_interfaces, only : hc1e, hc1v, hceii, hcvii, hceee, hcvvv
  use mclust_linalg, only : symmetric_eigen
  implicit none
  private

  type, public :: hc_result
    character(len=3) :: model_name = ''
    integer :: n = 0
    integer :: d = 0
    integer :: n_initial = 0
    integer :: minclus = 1
    integer, allocatable :: initial_partition(:)
    integer, allocatable :: merge(:,:)
  end type hc_result

  public :: hc_fit, hclass, hc_responsibilities
contains

  subroutine hc_fit(x, model, hc, minclus, partition, alpha, beta, status)
    real(dp), intent(in) :: x(:,:)
    character(len=*), intent(in) :: model
    type(hc_result), intent(out) :: hc
    integer, intent(in), optional :: minclus
    integer, intent(in), optional :: partition(:)
    real(dp), intent(in), optional :: alpha, beta
    integer, intent(out), optional :: status
    real(dp), allocatable :: xx(:,:), work(:), v(:), u(:,:), s(:,:), r(:,:)
    integer, allocatable :: ic(:), io(:), jo(:)
    integer :: n,p,ng,ns,nd,mn,i,k,info
    real(dp) :: a,b,tr
    character(len=3) :: mdl

    n=size(x,1); p=size(x,2); mdl=adjustl(model); mn=1
    if(present(minclus)) mn=minclus
    info=0
    if(n<2 .or. p<1 .or. mn<1) info=-1
    allocate(ic(n))
    if(present(partition)) then
      if(size(partition)/=n) then
        info=-2
      else
        call consecutive_partition(partition,ic,ng)
      end if
    else
      ic=[(i,i=1,n)]; ng=n
    end if
    if(info/=0) then
      if(present(status)) status=info
      return
    end if
    ns=ng-mn
    if(ns<1) then
      if(present(status)) status=-3
      return
    end if
    hc%n=n; hc%d=p; hc%n_initial=ng; hc%minclus=mn; hc%model_name=mdl
    allocate(hc%initial_partition(n),hc%merge(2,ns)); hc%initial_partition=ic; hc%merge=0
    allocate(xx(n,max(p+1,4))); xx=0.0_dp; xx(:,1:p)=x

    select case(trim(mdl))
    case('E')
      if(p/=1) then; info=-4; goto 900; end if
      nd=max(ng*(ng-1)/2,3*ns); allocate(work(max(1,nd))); work=0.0_dp
      call hc1e(xx(:,1),n,ic,ng,ns,nd,work)
      do k=1,ns
        hc%merge(1,k)=nint(xx(k,1)); hc%merge(2,k)=ic(k)
      end do
    case('V')
      if(p/=1) then; info=-4; goto 900; end if
      nd=max(ng*(ng-1)/2,3*ns); allocate(work(max(1,nd))); work=0.0_dp
      a=1.0_dp; if(present(alpha)) a=alpha
      tr=sum((x(:,1)-sum(x(:,1))/real(n,dp))**2)/real(n,dp)
      a=max(a*tr,epsilon(1.0_dp))
      call hc1v(xx(:,1),n,ic,ng,ns,a,nd,work)
      do k=1,ns
        hc%merge(1,k)=nint(xx(k,1)); hc%merge(2,k)=ic(k)
      end do
    case('EII')
      nd=max(ng*(ng-1)/2,3*ns); allocate(work(max(1,nd)),v(p)); work=0.0_dp; v=0.0_dp
      call hceii(xx,n,p,ic,ng,ns,v,nd,work)
      do k=1,ns
        hc%merge(:,k)=nint(xx(k,1:2))
      end do
    case('VII')
      nd=max(n,ng*(ng-1)/2,3*ns); allocate(work(max(1,nd)),v(p)); work=0.0_dp; v=0.0_dp
      a=1.0_dp; if(present(alpha)) a=alpha
      tr=trace_centered(x)/real(n*p,dp); a=max(a*tr,epsilon(1.0_dp))
      call hcvii(xx,n,p,ic,ng,ns,a,v,nd,work)
      do k=1,ns
        hc%merge(:,k)=nint(xx(k,1:2))
      end do
    case('VVV')
      nd=max(n,ng*(ng-1)/2+1,3*ns); allocate(work(max(1,nd)),v(p),u(p,p),s(p,p),r(p,p))
      work=0.0_dp; v=0.0_dp; u=0.0_dp; s=0.0_dp; r=0.0_dp
      a=1.0_dp; if(present(alpha)) a=alpha
      b=1.0_dp; if(present(beta)) b=beta
      tr=trace_centered(x)/real(n*p,dp); a=max(a*tr,epsilon(1.0_dp))
      call hcvvv(xx,n,p,ic,ng,ns,a,b,v,u,s,r,nd,work)
      do k=1,ns
        hc%merge(:,k)=nint(xx(k,1:2))
      end do
    case('EEE')
      allocate(io(max(1,ns)),jo(max(1,ns)),v(p),s(p,p),u(p,p),r(p,p))
      io=0; jo=0; v=0.0_dp; s=0.0_dp; u=0.0_dp; r=0.0_dp
      call hceee(xx,n,p,ic,ng,ns,io,jo,v,s,u,r)
      if(p<3) then
        hc%merge(1,:)=io(1:ns); hc%merge(2,:)=jo(1:ns)
      else if(p<4) then
        do k=1,ns; hc%merge(1,k)=nint(xx(k+1,3)); hc%merge(2,k)=jo(k); end do
      else
        do k=1,ns; hc%merge(:,k)=nint(xx(k+1,3:4)); end do
      end if
    case default
      info=-5
    end select
900 continue
    if(present(status)) status=info
  end subroutine hc_fit

  subroutine hclass(hc,g,classification,status)
    type(hc_result),intent(in)::hc
    integer,intent(in)::g
    integer,allocatable,intent(out)::classification(:)
    integer,intent(out),optional::status
    integer,allocatable::cl(:),values(:)
    integer::l,i,j,k,ng,info
    info=0; ng=hc%n_initial
    if(g<1 .or. g>ng .or. .not.allocated(hc%merge)) then
      allocate(classification(0)); if(present(status))status=-1; return
    end if
    allocate(cl(hc%n)); cl=hc%initial_partition
    do l=1,ng-g
      i=min(hc%merge(1,l),hc%merge(2,l)); j=max(hc%merge(1,l),hc%merge(2,l))
      where(cl==j) cl=i
    end do
    allocate(classification(hc%n),values(hc%n)); classification=0; values=0
    ! Relabel in first-occurrence order, as partconv(..., consec=TRUE).
    k=0
    do i=1,hc%n
      if(k==0 .or. .not.any(values(1:k)==cl(i))) then
        k=k+1; values(k)=cl(i)
      end if
    end do
    do i=1,k
      where(cl==values(i)) classification=i
    end do
    if(present(status)) status=info
  end subroutine hclass

  subroutine hc_responsibilities(x,g,z,model,status,use)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g
    real(dp),allocatable,intent(out)::z(:,:)
    character(len=*),intent(in),optional::model
    integer,intent(out),optional::status
    character(len=*),intent(in),optional::use
    type(hc_result)::hc
    integer,allocatable::cl(:)
    real(dp),allocatable::xt(:,:)
    character(len=3)::mdl
    character(len=4)::how
    integer::i,info
    if(size(x,2)==1) then
      call quantile_classes(x(:,1),g,cl)
      info=0
    else
      mdl='VVV'; if(present(model)) mdl=adjustl(model)
      how='VARS'; if(present(use)) how=adjustl(use)
      call hc_transform(x,how,xt,info)
      if(info/=0) then
        allocate(z(0,0)); if(present(status))status=info; return
      end if
      if(size(xt,1)<=size(xt,2) .and. trim(mdl)=='VVV') mdl='EII'
      call hc_fit(xt,mdl,hc,minclus=g,status=info)
      if(info==0) call hclass(hc,g,cl,info)
    end if
    if(info/=0) then
      allocate(z(0,0)); if(present(status)) status=info; return
    end if
    allocate(z(size(x,1),g)); z=0.0_dp
    do i=1,size(x,1)
      if(cl(i)>=1 .and. cl(i)<=g) z(i,cl(i))=1.0_dp
    end do
    if(present(status)) status=0
  end subroutine hc_responsibilities

  subroutine hc_transform(x,use,z,status)
    real(dp),intent(in)::x(:,:)
    character(len=*),intent(in)::use
    real(dp),allocatable,intent(out)::z(:,:)
    integer,intent(out)::status
    real(dp),allocatable::xc(:,:),xs(:,:),cross(:,:),eval(:),vec(:,:),sd(:),tmp(:,:)
    logical,allocatable::keep(:)
    real(dp)::den,scale
    integer::n,p,j,k,info,q
    character(len=4)::how
    n=size(x,1); p=size(x,2); how=adjustl(use); status=0
    if(trim(how)=='VARS') then
      allocate(z(n,p)); z=x; return
    end if
    allocate(xc(n,p),sd(p),keep(p)); xc=x-spread(sum(x,dim=1)/real(n,dp),1,n)
    den=real(max(1,n-1),dp)
    do j=1,p; sd(j)=sqrt(sum(xc(:,j)**2)/den); end do
    keep=sd>sqrt(epsilon(1.0_dp))*max(1.0_dp,maxval(sd))
    if(.not.any(keep)) then; allocate(z(0,0)); status=-10; return; end if
    q=count(keep); allocate(xs(n,q)); k=0
    do j=1,p
      if(.not.keep(j))cycle
      k=k+1
      if(trim(how)=='STD' .or. trim(how)=='PCR' .or. trim(how)=='SVD') then
        xs(:,k)=xc(:,j)/sd(j)
      else
        xs(:,k)=xc(:,j)
      end if
    end do
    if(trim(how)=='STD') then; call move_alloc(xs,z); return; end if
    allocate(cross(q,q)); cross=matmul(transpose(xs),xs)
    if(trim(how)=='SPH') cross=cross/real(n,dp)
    call symmetric_eigen(cross,eval,vec,info)
    if(info/=0) then; allocate(z(0,0)); status=info; return; end if
    select case(trim(how))
    case('PCS','PCR')
      z=matmul(xs,vec)
    case('SVD')
      allocate(tmp(n,q)); tmp=matmul(xs,vec)
      do j=1,q
        ! singular value is sqrt(eigenvalue); R divides by sqrt(singular value).
        scale=max(eval(j),tiny(1.0_dp))**0.25_dp
        tmp(:,j)=tmp(:,j)/scale
      end do
      call move_alloc(tmp,z)
    case('SPH')
      allocate(tmp(n,q)); tmp=matmul(xs,vec)
      do j=1,q
        scale=sqrt(max(eval(j),tiny(1.0_dp))); tmp(:,j)=tmp(:,j)/scale
      end do
      call move_alloc(tmp,z)
    case default
      allocate(z(0,0)); status=-11
    end select
  end subroutine hc_transform

  subroutine consecutive_partition(partition,out,ng)
    integer,intent(in)::partition(:)
    integer,intent(out)::out(:)
    integer,intent(out)::ng
    integer,allocatable::values(:)
    integer::i,j
    allocate(values(size(partition))); values=0; ng=0
    do i=1,size(partition)
      j=0
      if(ng>0) then
        do j=1,ng
          if(values(j)==partition(i)) exit
        end do
        if(j>ng) j=0
      end if
      if(j==0) then; ng=ng+1; values(ng)=partition(i); j=ng; end if
      out(i)=j
    end do
  end subroutine consecutive_partition

  pure real(dp) function trace_centered(x) result(v)
    real(dp),intent(in)::x(:,:)
    real(dp)::mu(size(x,2))
    integer::i
    mu=sum(x,dim=1)/real(size(x,1),dp); v=0.0_dp
    do i=1,size(x,1); v=v+sum((x(i,:)-mu)**2); end do
  end function trace_centered

  subroutine quantile_classes(x,g,cl)
    real(dp),intent(in)::x(:)
    integer,intent(in)::g
    integer,allocatable,intent(out)::cl(:)
    integer,allocatable::ord(:)
    integer::i,k,n
    n=size(x); allocate(ord(n),cl(n)); ord=[(i,i=1,n)]
    call sort_index(x,ord)
    do i=1,n
      k=min(g,1+int(real((i-1)*g,dp)/real(n,dp))); cl(ord(i))=k
    end do
  end subroutine quantile_classes

  subroutine sort_index(x,idx)
    real(dp),intent(in)::x(:)
    integer,intent(inout)::idx(:)
    integer::i,j,t
    do i=2,size(idx)
      t=idx(i); j=i-1
      do while(j>=1)
        if(x(idx(j))<=x(t)) exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=t
    end do
  end subroutine sort_index
end module mclust_hierarchical
