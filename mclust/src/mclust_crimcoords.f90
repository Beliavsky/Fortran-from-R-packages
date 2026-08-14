! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_crimcoords
  use mclust_kinds, only : dp
  use mclust_linalg, only : inverse_sqrt_symmetric, symmetric_eigen
  use mclust_utilities, only : unmap_classes
  implicit none
  private
  type,public :: crimcoords_fit
    integer :: n=0,d=0,g=0,numdir=0,status=0
    real(dp),allocatable :: means(:,:),between(:,:),within(:,:)
    real(dp),allocatable :: evalues(:),basis(:,:),projection(:,:)
    integer,allocatable :: classification(:),levels(:)
  end type crimcoords_fit
  public :: fit_crimcoords
contains
  subroutine fit_crimcoords(x,classification,result,numdir,unbiased,status)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::classification(:)
    type(crimcoords_fit),intent(out)::result
    integer,intent(in),optional::numdir
    logical,intent(in),optional::unbiased
    integer,intent(out),optional::status
    real(dp),allocatable::z(:,:),winvhalf(:,:),eval(:),evec(:,:),a(:,:)
    real(dp),allocatable::overall(:),group_mean(:,:),centered(:,:)
    integer,allocatable::levels(:),nk(:)
    logical::unb
    integer::n,d,g,i,k,nd,info
    n=size(x,1); d=size(x,2); unb=.false.; if(present(unbiased)) unb=unbiased
    if(size(classification)/=n .or. n<2 .or. d<1) then
      result%status=-1; if(present(status)) status=-1; return
    end if
    call unmap_classes(classification,z,levels); g=size(z,2)
    allocate(nk(g),overall(d),group_mean(g,d),centered(n,d))
    do k=1,g; nk(k)=count(classification==levels(k)); end do
    if(any(nk==0) .or. (unb .and. (n<=g .or. g<=1))) then
      result%status=-2; if(present(status)) status=-2; return
    end if
    overall=sum(x,dim=1)/real(n,dp)
    group_mean=0.0_dp
    do k=1,g
      do i=1,n
        if(classification(i)==levels(k)) group_mean(k,:)=group_mean(k,:)+x(i,:)
      end do
      group_mean(k,:)=group_mean(k,:)/real(nk(k),dp)
    end do
    allocate(result%within(d,d),result%between(d,d)); result%within=0.0_dp; result%between=0.0_dp
    do i=1,n
      k=find_level(levels,classification(i)); centered(i,:)=x(i,:)-group_mean(k,:)
      result%within=result%within+outer(centered(i,:),centered(i,:))
    end do
    do k=1,g
      result%between=result%between+real(nk(k),dp)*outer(group_mean(k,:)-overall,group_mean(k,:)-overall)
    end do
    if(unb) then
      result%within=result%within/real(n-g,dp)
      result%between=result%between/real(g-1,dp)
    else
      result%within=result%within/real(n,dp)
      result%between=result%between/real(g,dp)
    end if
    call inverse_sqrt_symmetric(result%within,winvhalf,status=info)
    if(info/=0 .or. size(winvhalf,1)==0) then
      result%status=10+info; if(present(status)) status=result%status; return
    end if
    a=matmul(winvhalf,matmul(result%between,winvhalf))
    call symmetric_eigen(a,eval,evec,info)
    if(info/=0) then; result%status=20+info; if(present(status)) status=result%status; return; end if
    eval=max(eval,0.0_dp)
    if(present(numdir)) then
      nd=max(0,min(d,numdir))
    else
      nd=count(eval>sqrt(epsilon(1.0_dp)))
      nd=min(d,nd)
    end if
    result%n=n; result%d=d; result%g=g; result%numdir=nd; result%status=0
    allocate(result%means(g,d),result%evalues(d),result%basis(d,nd),result%projection(n,nd))
    allocate(result%classification(n),result%levels(g))
    result%means=group_mean; result%evalues=eval; result%classification=classification; result%levels=levels
    if(nd>0) then
      result%basis=matmul(winvhalf,evec(:,1:nd))
      result%projection=matmul(x,result%basis)
      do k=1,nd
        if(median(result%projection(:,k))<0.0_dp) then
          result%basis(:,k)=-result%basis(:,k); result%projection(:,k)=-result%projection(:,k)
        end if
      end do
    end if
    if(present(status)) status=0
  end subroutine fit_crimcoords

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i
    do i=1,size(a); c(i,:)=a(i)*b; end do
  end function outer

  pure integer function find_level(levels,value) result(k)
    integer,intent(in)::levels(:),value
    integer::i
    k=0
    do i=1,size(levels)
      if(levels(i)==value) then; k=i; return; end if
    end do
  end function find_level

  real(dp) function median(x) result(m)
    real(dp),intent(in)::x(:)
    real(dp),allocatable::a(:)
    real(dp)::key
    integer::i,j,n
    n=size(x); allocate(a(n)); a=x
    do i=2,n
      key=a(i); j=i-1
      do while(j>=1)
        if(a(j)<=key) exit
        a(j+1)=a(j); j=j-1
      end do
      a(j+1)=key
    end do
    if(mod(n,2)==1) then
      m=a((n+1)/2)
    else
      m=0.5_dp*(a(n/2)+a(n/2+1))
    end if
  end function median
end module mclust_crimcoords
