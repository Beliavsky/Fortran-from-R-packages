! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_interp
  use fbasics_kinds, only: dp
  implicit none
  private
  public :: linear_interp, bilinear_interp, idw_interp, local_plane_interp
contains
  real(dp) function linear_interp(x,y,xout) result(v)
    real(dp),intent(in)::x(:),y(:),xout
    integer::i,n
    n=size(x)
    if(n/=size(y).or.n==0)then;v=0.0_dp;return;end if
    if(xout<=x(1))then;v=y(1);return;end if
    if(xout>=x(n))then;v=y(n);return;end if
    do i=1,n-1
      if(xout>=x(i).and.xout<=x(i+1))then
        v=y(i)+(y(i+1)-y(i))*(xout-x(i))/(x(i+1)-x(i));return
      end if
    end do
    v=y(n)
  end function

  real(dp) function bilinear_interp(x,y,z,xout,yout) result(v)
    real(dp),intent(in)::x(:),y(:),z(:,:),xout,yout
    integer::i,j,nx,ny
    real(dp)::tx,ty
    nx=size(x);ny=size(y)
    i=1;do while(i<nx-1.and.xout>x(i+1));i=i+1;end do
    j=1;do while(j<ny-1.and.yout>y(j+1));j=j+1;end do
    tx=(min(max(xout,x(i)),x(i+1))-x(i))/max(x(i+1)-x(i),epsilon(1.0_dp))
    ty=(min(max(yout,y(j)),y(j+1))-y(j))/max(y(j+1)-y(j),epsilon(1.0_dp))
    v=(1.0_dp-tx)*(1.0_dp-ty)*z(i,j)+tx*(1.0_dp-ty)*z(i+1,j)+ &
      (1.0_dp-tx)*ty*z(i,j+1)+tx*ty*z(i+1,j+1)
  end function

  real(dp) function idw_interp(x,y,z,xout,yout,power,k_nearest) result(v)
    real(dp),intent(in)::x(:),y(:),z(:),xout,yout
    real(dp),intent(in),optional::power
    integer,intent(in),optional::k_nearest
    real(dp)::pw,d,w,sumw
    real(dp),allocatable::dist(:)
    integer,allocatable::idx(:)
    integer::i,j,k,n,kk,it
    n=size(x);pw=2.0_dp;if(present(power))pw=power;kk=n;if(present(k_nearest))kk=min(k_nearest,n)
    allocate(dist(n),idx(n));do i=1,n;dist(i)=hypot(x(i)-xout,y(i)-yout);idx(i)=i;end do
    do i=2,n
      d=dist(i);it=idx(i);j=i-1
      do while(j>=1)
        if(dist(j)<=d)exit
        dist(j+1)=dist(j);idx(j+1)=idx(j);j=j-1
      end do
      dist(j+1)=d;idx(j+1)=it
    end do
    if(dist(1)<=epsilon(1.0_dp))then;v=z(idx(1));return;end if
    v=0.0_dp;sumw=0.0_dp
    do k=1,kk;w=1.0_dp/max(dist(k),epsilon(1.0_dp))**pw;v=v+w*z(idx(k));sumw=sumw+w;end do
    v=v/sumw
  end function

  real(dp) function local_plane_interp(x,y,z,xout,yout,k_nearest) result(v)
    real(dp),intent(in)::x(:),y(:),z(:),xout,yout
    integer,intent(in),optional::k_nearest
    integer::n,k,i,j,it,info
    integer,allocatable::idx(:)
    real(dp),allocatable::dist(:),a(:,:),ata(:,:),rhs(:),ainv(:,:)
    real(dp)::d
    n=size(x);k=min(12,n);if(present(k_nearest))k=min(k_nearest,n)
    allocate(dist(n),idx(n));do i=1,n;dist(i)=hypot(x(i)-xout,y(i)-yout);idx(i)=i;end do
    do i=2,n;d=dist(i);it=idx(i);j=i-1;do while(j>=1);if(dist(j)<=d)exit;dist(j+1)=dist(j);idx(j+1)=idx(j);j=j-1;end do;dist(j+1)=d;idx(j+1)=it;end do
    if(dist(1)<=epsilon(1.0_dp))then;v=z(idx(1));return;end if
    allocate(a(k,3),ata(3,3),rhs(3));do i=1,k;a(i,:)=[1.0_dp,x(idx(i))-xout,y(idx(i))-yout];end do
    ata=matmul(transpose(a),a);rhs=matmul(transpose(a),z(idx(1:k)))
    call inverse3(ata,ainv,info)
    if(info==0)then;v=dot_product(ainv(1,:),rhs);else;v=idw_interp(x,y,z,xout,yout,2.0_dp,k);end if
  contains
    subroutine inverse3(m,mi,istat)
      real(dp),intent(in)::m(3,3);real(dp),allocatable,intent(out)::mi(:,:);integer,intent(out)::istat
      real(dp)::det
      allocate(mi(3,3));det=m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2))-m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1))+m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
      if(abs(det)<1.0e-14_dp)then;istat=1;mi=0.0_dp;return;end if
      mi(1,1)=(m(2,2)*m(3,3)-m(2,3)*m(3,2))/det;mi(1,2)=(m(1,3)*m(3,2)-m(1,2)*m(3,3))/det;mi(1,3)=(m(1,2)*m(2,3)-m(1,3)*m(2,2))/det
      mi(2,1)=(m(2,3)*m(3,1)-m(2,1)*m(3,3))/det;mi(2,2)=(m(1,1)*m(3,3)-m(1,3)*m(3,1))/det;mi(2,3)=(m(1,3)*m(2,1)-m(1,1)*m(2,3))/det
      mi(3,1)=(m(2,1)*m(3,2)-m(2,2)*m(3,1))/det;mi(3,2)=(m(1,2)*m(3,1)-m(1,1)*m(3,2))/det;mi(3,3)=(m(1,1)*m(2,2)-m(1,2)*m(2,1))/det;istat=0
    end subroutine
  end function
end module fbasics_interp
