! SPDX-License-Identifier: GPL-3.0-or-later
! Natural-spline basis compatible with splines2::naturalSpline(intercept=FALSE)
! as used by flexsurv's spline="splines2ns" option.  The implementation is a
! native Fortran implementation of the natural cubic B-spline construction.
module flexsurv_splines2ns
  use flexsurv_kinds, only : dp
  implicit none
  private
  public :: splines2ns_basis, splines2ns_dbasis
  public :: splines2ns_basis_matrix, splines2ns_dbasis_matrix

contains

  pure subroutine splines2ns_basis(knots,x,b)
    real(dp),intent(in)::knots(:),x
    real(dp),intent(out)::b(size(knots))
    real(dp),allocatable::bc(:),nullmat(:,:),nat(:)
    integer::ncol
    ncol=size(knots);b=0.0_dp
    if(ncol<2)return
    allocate(bc(ncol+2),nullmat(ncol+2,ncol),nat(ncol))
    call complete_bspline(knots,x,0,bc)
    call natural_null_matrix(knots,nullmat)
    nat=matmul(bc,nullmat)
    b(1)=1.0_dp
    if(ncol>1)b(2:)=nat(2:)
  end subroutine splines2ns_basis

  pure subroutine splines2ns_dbasis(knots,x,b)
    real(dp),intent(in)::knots(:),x
    real(dp),intent(out)::b(size(knots))
    real(dp),allocatable::bc(:),nullmat(:,:),nat(:)
    integer::ncol
    ncol=size(knots);b=0.0_dp
    if(ncol<2)return
    allocate(bc(ncol+2),nullmat(ncol+2,ncol),nat(ncol))
    call complete_bspline(knots,x,1,bc)
    call natural_null_matrix(knots,nullmat)
    nat=matmul(bc,nullmat)
    b(1)=0.0_dp
    if(ncol>1)b(2:)=nat(2:)
  end subroutine splines2ns_dbasis

  pure subroutine splines2ns_basis_matrix(knots,x,b)
    real(dp),intent(in)::knots(:),x(:)
    real(dp),intent(out)::b(size(x),size(knots))
    integer::i
    do i=1,size(x)
      call splines2ns_basis(knots,x(i),b(i,:))
    end do
  end subroutine splines2ns_basis_matrix

  pure subroutine splines2ns_dbasis_matrix(knots,x,b)
    real(dp),intent(in)::knots(:),x(:)
    real(dp),intent(out)::b(size(x),size(knots))
    integer::i
    do i=1,size(x)
      call splines2ns_dbasis(knots,x(i),b(i,:))
    end do
  end subroutine splines2ns_dbasis_matrix

  pure subroutine complete_bspline(knots,x,deriv,b)
    real(dp),intent(in)::knots(:),x
    integer,intent(in)::deriv
    real(dp),intent(out)::b(size(knots)+2)
    real(dp),allocatable::kv(:),bv(:),bd(:)
    real(dp)::xx,left,right,epsx
    integer::m,nbase
    m=size(knots)-2;nbase=m+4
    allocate(kv(m+8),bv(nbase),bd(nbase))
    call cubic_knot_vector(knots,kv)
    left=knots(1);right=knots(size(knots));xx=x
    epsx=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(left),abs(right))
    if(x<left)then
      call cubic_bspline_value_derivative(kv,left+epsx,bv,bd)
      if(deriv==0)then;b=bv+bd*(x-left);else;b=bd;end if
    else if(x>right)then
      call cubic_bspline_value_derivative(kv,right-epsx,bv,bd)
      if(deriv==0)then;b=bv+bd*(x-right);else;b=bd;end if
    else
      if(x<=left)xx=left+epsx
      if(x>=right)xx=right-epsx
      call cubic_bspline_value_derivative(kv,xx,bv,bd)
      if(deriv==0)then;b=bv;else;b=bd;end if
    end if
  end subroutine complete_bspline

  pure subroutine cubic_knot_vector(knots,kv)
    real(dp),intent(in)::knots(:)
    real(dp),intent(out)::kv(size(knots)+6)
    integer::m
    m=size(knots)-2
    kv(1:4)=knots(1)
    if(m>0)kv(5:4+m)=knots(2:size(knots)-1)
    kv(5+m:8+m)=knots(size(knots))
  end subroutine cubic_knot_vector

  pure subroutine cubic_bspline_value_derivative(kv,x,b,db)
    real(dp),intent(in)::kv(:),x
    real(dp),intent(out)::b(size(kv)-4),db(size(kv)-4)
    real(dp),allocatable::n0(:),n1(:),n2(:),n3(:)
    real(dp)::d1,d2
    integer::i,nbase
    nbase=size(kv)-4
    allocate(n0(nbase+3),n1(nbase+2),n2(nbase+1),n3(nbase))
    n0=0.0_dp
    do i=1,nbase+3
      if(kv(i)<=x.and.x<kv(i+1))n0(i)=1.0_dp
    end do
    call bspline_step(kv,x,n0,1,n1)
    call bspline_step(kv,x,n1,2,n2)
    call bspline_step(kv,x,n2,3,n3)
    b=n3;db=0.0_dp
    do i=1,nbase
      d1=kv(i+3)-kv(i)
      d2=kv(i+4)-kv(i+1)
      if(d1>0.0_dp)db(i)=db(i)+3.0_dp*n2(i)/d1
      if(d2>0.0_dp)db(i)=db(i)-3.0_dp*n2(i+1)/d2
    end do
  end subroutine cubic_bspline_value_derivative

  pure subroutine bspline_step(kv,x,prev,degree,out)
    real(dp),intent(in)::kv(:),x,prev(:)
    integer,intent(in)::degree
    real(dp),intent(out)::out(size(prev)-1)
    real(dp)::d1,d2,a,c
    integer::i
    out=0.0_dp
    do i=1,size(out)
      d1=kv(i+degree)-kv(i)
      d2=kv(i+degree+1)-kv(i+1)
      a=0.0_dp;c=0.0_dp
      if(d1>0.0_dp)a=(x-kv(i))*prev(i)/d1
      if(d2>0.0_dp)c=(kv(i+degree+1)-x)*prev(i+1)/d2
      out(i)=a+c
    end do
  end subroutine bspline_step

  pure subroutine natural_null_matrix(knots,a)
    real(dp),intent(in)::knots(:)
    real(dp),intent(out)::a(size(knots)+2,size(knots))
    integer::m,ncol,i,j
    real(dp)::w1,w2,w3,s
    m=size(knots)-2;ncol=m+2;a=0.0_dp
    if(m==0)then
      a(1,1)=3.0_dp;a(2,1)=2.0_dp;a(3,1)=1.0_dp
      a(2,2)=1.0_dp;a(3,2)=2.0_dp;a(4,2)=3.0_dp
    else if(m==1)then
      w1=knots(2)-knots(1);w2=knots(3)-knots(2);w3=w1+w2
      a(1,1)=1.0_dp+w1/w3;a(2,1)=1.0_dp
      a(2,2)=w1/(w1+w3);a(3,2)=1.0_dp;a(4,2)=w2/(w2+w3)
      a(4,3)=1.0_dp;a(5,3)=1.0_dp+w2/w3
    else
      do i=0,2
        a(i+1,1)=1.0_dp
        a(ncol-i+2,ncol)=1.0_dp
      end do
      a(2,2)=1.0_dp
      a(3,2)=1.0_dp+(knots(3)-knots(1))/(knots(2)-knots(1))
      a(ncol,ncol-1)=1.0_dp+(knots(size(knots))-knots(size(knots)-2))/ &
        (knots(size(knots))-knots(size(knots)-1))
      a(ncol+1,ncol-1)=1.0_dp
      if(ncol>4)then
        do j=0,ncol-5
          a(j+4,j+3)=1.0_dp
        end do
      end if
    end if
    do j=1,ncol
      s=sum(a(:,j))
      if(s>0.0_dp)a(:,j)=a(:,j)/s
    end do
  end subroutine natural_null_matrix

end module flexsurv_splines2ns
