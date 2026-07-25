! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_spatial
  use fbasics_kinds, only: dp
  use fbasics_linalg, only: matrix_inverse
  use fbasics_stats, only: sample_variance, sample_median
  implicit none
  private
  type, public :: kriging_model
    real(dp) :: range=1.0_dp,sill=1.0_dp,nugget=0.0_dp
    character(len=12) :: covariance='exponential'
  end type kriging_model
  public :: estimate_kriging_model, ordinary_kriging, triangulated_interp
contains
  subroutine estimate_kriging_model(xy,z,model,covariance)
    real(dp),intent(in)::xy(:,:),z(:)
    type(kriging_model),intent(out)::model
    character(len=*),intent(in),optional::covariance
    real(dp),allocatable::dist(:)
    integer::i,j,k,n,np
    n=size(z);np=n*(n-1)/2;allocate(dist(max(1,np)));k=0
    do i=1,n-1;do j=i+1,n;k=k+1;dist(k)=sqrt(sum((xy(i,:)-xy(j,:))**2));end do;end do
    if(k>0)then;model%range=max(sample_median(dist(1:k)),1.0e-8_dp);else;model%range=1.0_dp;end if
    model%sill=max(sample_variance(z),1.0e-10_dp);model%nugget=1.0e-8_dp*model%sill
    if(present(covariance))model%covariance=adjustl(covariance)
  end subroutine estimate_kriging_model

  subroutine ordinary_kriging(xy,z,query,model,prediction,variance,weights,ok)
    real(dp),intent(in)::xy(:,:),z(:),query(:,:)
    type(kriging_model),intent(in)::model
    real(dp),allocatable,intent(out)::prediction(:),variance(:)
    real(dp),allocatable,intent(out),optional::weights(:,:)
    logical,intent(out),optional::ok
    real(dp),allocatable::a(:,:),ainv(:,:),rhs(:),sol(:),wout(:,:)
    real(dp)::d,c0
    integer::n,m,i,j,k
    integer::info
    n=size(z);m=size(query,1);allocate(a(n+1,n+1),rhs(n+1),sol(n+1),prediction(m),variance(m),wout(n,m));a=0.0_dp
    do i=1,n
      do j=1,n
        d=sqrt(sum((xy(i,:)-xy(j,:))**2));a(i,j)=covariance_value(d,model)
        if(i==j)a(i,j)=a(i,j)+model%nugget
      end do
      a(i,n+1)=1.0_dp;a(n+1,i)=1.0_dp
    end do
    call matrix_inverse(a,ainv,info)
    if(info/=0)then;prediction=0.0_dp;variance=huge(1.0_dp);if(present(weights))weights=0.0_dp;if(present(ok))ok=.false.;return;end if
    c0=model%sill+model%nugget
    do k=1,m
      do i=1,n;d=sqrt(sum((xy(i,:)-query(k,:))**2));rhs(i)=covariance_value(d,model);end do;rhs(n+1)=1.0_dp;sol=matmul(ainv,rhs);wout(:,k)=sol(1:n);prediction(k)=dot_product(sol(1:n),z);variance(k)=max(0.0_dp,c0-dot_product(sol(1:n),rhs(1:n))+sol(n+1))
    end do
    if(present(weights))weights=wout;if(present(ok))ok=.true.
  end subroutine ordinary_kriging

  pure real(dp) function covariance_value(distance,model) result(v)
    real(dp),intent(in)::distance;type(kriging_model),intent(in)::model;real(dp)::h
    h=distance/max(model%range,1.0e-12_dp)
    select case(trim(model%covariance))
    case('gaussian');v=model%sill*exp(-h*h)
    case('spherical');if(h>=1.0_dp)then;v=0.0_dp;else;v=model%sill*(1.0_dp-1.5_dp*h+0.5_dp*h**3);end if
    case default;v=model%sill*exp(-h)
    end select
  end function covariance_value

  subroutine triangulated_interp(xy,z,query,value,inside,triangle_indices)
    real(dp),intent(in)::xy(:,:),z(:),query(:,:)
    real(dp),allocatable,intent(out)::value(:)
    logical,allocatable,intent(out),optional::inside(:)
    integer,allocatable,intent(out),optional::triangle_indices(:,:)
    integer::n,m,i,j,k,q,besti,bestj,bestk
    real(dp)::w1,w2,w3,score,bestscore,val
    logical::found
    logical,allocatable::inout(:);integer,allocatable::tout(:,:)
    n=size(z);m=size(query,1);allocate(value(m),inout(m),tout(3,m));inout=.false.;tout=0
    do q=1,m
      found=.false.;bestscore=huge(1.0_dp);val=0.0_dp;besti=0;bestj=0;bestk=0
      do i=1,n-2;do j=i+1,n-1;do k=j+1,n
        call barycentric(query(q,:),xy(i,:),xy(j,:),xy(k,:),w1,w2,w3)
        if(min(w1,min(w2,w3))>=-1.0e-10_dp.and.max(w1,max(w2,w3))<=1.0_dp+1.0e-10_dp)then
          score=sum((xy(i,:)-xy(j,:))**2)+sum((xy(i,:)-xy(k,:))**2)+sum((xy(j,:)-xy(k,:))**2)
          if(score<bestscore)then;bestscore=score;val=w1*z(i)+w2*z(j)+w3*z(k);besti=i;bestj=j;bestk=k;found=.true.;end if
        end if
      end do;end do;end do
      if(.not.found)call nearest_plane_value(xy,z,query(q,:),val,besti,bestj,bestk)
      value(q)=val;inout(q)=found;tout(:,q)=[besti,bestj,bestk]
    end do
    if(present(inside))inside=inout;if(present(triangle_indices))triangle_indices=tout
  end subroutine triangulated_interp

  subroutine barycentric(p,a,b,c,w1,w2,w3)
    real(dp),intent(in)::p(:),a(:),b(:),c(:);real(dp),intent(out)::w1,w2,w3;real(dp)::den
    den=(b(2)-c(2))*(a(1)-c(1))+(c(1)-b(1))*(a(2)-c(2))
    if(abs(den)<1.0e-14_dp)then;w1=-huge(1.0_dp);w2=w1;w3=w1;return;end if
    w1=((b(2)-c(2))*(p(1)-c(1))+(c(1)-b(1))*(p(2)-c(2)))/den
    w2=((c(2)-a(2))*(p(1)-c(1))+(a(1)-c(1))*(p(2)-c(2)))/den;w3=1.0_dp-w1-w2
  end subroutine barycentric

  subroutine nearest_plane_value(xy,z,p,val,i1,i2,i3)
    real(dp),intent(in)::xy(:,:),z(:),p(:);real(dp),intent(out)::val;integer,intent(out)::i1,i2,i3
    real(dp),allocatable::d(:);real(dp)::w1,w2,w3,ws
    integer::i,n
    n=size(z);allocate(d(n));do i=1,n;d(i)=sqrt(sum((xy(i,:)-p)**2));end do
    i1=minloc(d,dim=1);d(i1)=huge(1.0_dp);i2=minloc(d,dim=1);d(i2)=huge(1.0_dp);i3=minloc(d,dim=1)
    w1=1.0_dp/max(sqrt(sum((xy(i1,:)-p)**2)),1.0e-12_dp);w2=1.0_dp/max(sqrt(sum((xy(i2,:)-p)**2)),1.0e-12_dp);w3=1.0_dp/max(sqrt(sum((xy(i3,:)-p)**2)),1.0e-12_dp);ws=w1+w2+w3;val=(w1*z(i1)+w2*z(i2)+w3*z(i3))/ws
  end subroutine nearest_plane_value
end module fbasics_spatial
