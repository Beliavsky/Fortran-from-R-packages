! SPDX-License-Identifier: BSD-3-Clause
! Minimal BLAS level-1 compatibility routines used by the vendored solvers.
! Straightforward reference implementations; no external BLAS is required.
subroutine dcopy(n, dx, incx, dy, incy)
  implicit none
  integer, intent(in) :: n, incx, incy
  double precision, intent(in) :: dx(*)
  double precision, intent(out) :: dy(*)
  integer :: i,ix,iy
  ix=1;iy=1
  do i=1,n;dy(iy)=dx(ix);ix=ix+incx;iy=iy+incy;end do
end subroutine dcopy

subroutine daxpy(n, da, dx, incx, dy, incy)
  implicit none
  integer,intent(in)::n,incx,incy
  double precision,intent(in)::da,dx(*)
  double precision,intent(inout)::dy(*)
  integer::i,ix,iy
  ix=1;iy=1
  do i=1,n;dy(iy)=dy(iy)+da*dx(ix);ix=ix+incx;iy=iy+incy;end do
end subroutine daxpy

double precision function ddot(n,dx,incx,dy,incy)
  implicit none
  integer,intent(in)::n,incx,incy
  double precision,intent(in)::dx(*),dy(*)
  integer::i,ix,iy
  ddot=0.0d0;ix=1;iy=1
  do i=1,n;ddot=ddot+dx(ix)*dy(iy);ix=ix+incx;iy=iy+incy;end do
end function ddot

subroutine dscal(n,da,dx,incx)
  implicit none
  integer,intent(in)::n,incx
  double precision,intent(in)::da
  double precision,intent(inout)::dx(*)
  integer::i,ix
  ix=1
  do i=1,n;dx(ix)=da*dx(ix);ix=ix+incx;end do
end subroutine dscal

integer function idamax(n,dx,incx)
  implicit none
  integer,intent(in)::n,incx
  double precision,intent(in)::dx(*)
  integer::i,ix
  double precision::best
  if(n<=0)then;idamax=0;return;end if
  idamax=1;best=abs(dx(1));ix=1+incx
  do i=2,n
    if(abs(dx(ix))>best)then;best=abs(dx(ix));idamax=i;end if
    ix=ix+incx
  end do
end function idamax

double precision function dnrm2(n,x,incx)
  implicit none
  integer,intent(in)::n,incx
  double precision,intent(in)::x(*)
  integer::i,ix
  double precision::scale,ssq,ax
  if(n<=0)then;dnrm2=0.0d0;return;end if
  scale=0.0d0;ssq=1.0d0;ix=1
  do i=1,n
    if(x(ix)/=0.0d0)then
      ax=abs(x(ix))
      if(scale<ax)then;ssq=1.0d0+ssq*(scale/ax)**2;scale=ax
      else;ssq=ssq+(ax/scale)**2;end if
    end if
    ix=ix+incx
  end do
  dnrm2=scale*sqrt(ssq)
end function dnrm2

subroutine zcopy(n,zx,incx,zy,incy)
  implicit none
  integer,intent(in)::n,incx,incy
  complex(kind=kind(1.0d0)),intent(in)::zx(*)
  complex(kind=kind(1.0d0)),intent(out)::zy(*)
  integer::i,ix,iy
  ix=1;iy=1
  do i=1,n;zy(iy)=zx(ix);ix=ix+incx;iy=iy+incy;end do
end subroutine zcopy

subroutine zaxpy(n,za,zx,incx,zy,incy)
  implicit none
  integer,intent(in)::n,incx,incy
  complex(kind=kind(1.0d0)),intent(in)::za,zx(*)
  complex(kind=kind(1.0d0)),intent(inout)::zy(*)
  integer::i,ix,iy
  ix=1;iy=1
  do i=1,n;zy(iy)=zy(iy)+za*zx(ix);ix=ix+incx;iy=iy+incy;end do
end subroutine zaxpy

complex(kind=kind(1.0d0)) function zdotc(n,zx,incx,zy,incy)
  implicit none
  integer,intent(in)::n,incx,incy
  complex(kind=kind(1.0d0)),intent(in)::zx(*),zy(*)
  integer::i,ix,iy
  zdotc=(0.0d0,0.0d0);ix=1;iy=1
  do i=1,n;zdotc=zdotc+conjg(zx(ix))*zy(iy);ix=ix+incx;iy=iy+incy;end do
end function zdotc

subroutine zscal(n,za,zx,incx)
  implicit none
  integer,intent(in)::n,incx
  complex(kind=kind(1.0d0)),intent(in)::za
  complex(kind=kind(1.0d0)),intent(inout)::zx(*)
  integer::i,ix
  ix=1
  do i=1,n;zx(ix)=za*zx(ix);ix=ix+incx;end do
end subroutine zscal

integer function izamax(n,zx,incx)
  implicit none
  integer,intent(in)::n,incx
  complex(kind=kind(1.0d0)),intent(in)::zx(*)
  integer::i,ix
  double precision::best,cur
  if(n<=0)then;izamax=0;return;end if
  izamax=1;best=abs(dble(zx(1)))+abs(aimag(zx(1)));ix=1+incx
  do i=2,n
    cur=abs(dble(zx(ix)))+abs(aimag(zx(ix)))
    if(cur>best)then;best=cur;izamax=i;end if
    ix=ix+incx
  end do
end function izamax
