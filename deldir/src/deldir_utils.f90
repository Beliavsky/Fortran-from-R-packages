module deldir_utils
use deldir_kinds, only: dp
use deldir_kernel, only: binsrt
implicit none
private
public :: deldir_bin_sort, duplicated_xy, corner_indices, midpoint_inside, find_new_in_old
contains
subroutine deldir_bin_sort(x,y,rw,xsorted,ysorted,ind,rind)
    real(dp),intent(in)::x(:),y(:),rw(4)
    real(dp),allocatable,intent(out)::xsorted(:),ysorted(:)
    integer,allocatable,intent(out)::ind(:),rind(:)
    real(dp),allocatable::tx(:),ty(:)
    integer,allocatable::ilst(:)
    integer::n
    if(size(x)/=size(y))error stop 'deldir_bin_sort: length mismatch'
    n=size(x);allocate(xsorted(n),ysorted(n),tx(n),ty(n),ind(n),rind(n),ilst(n))
    xsorted=x;ysorted=y
    call binsrt(xsorted,ysorted,rw,n,ind,rind,tx,ty,ilst)
end subroutine
subroutine duplicated_xy(x,y,mask)
    real(dp),intent(in)::x(:),y(:)
    logical,allocatable,intent(out)::mask(:)
    integer::i,j
    if(size(x)/=size(y))error stop 'duplicated_xy: length mismatch'
    allocate(mask(size(x)));mask=.false.
    do i=2,size(x)
        do j=1,i-1
            if(x(i)<=x(j).and.x(i)>=x(j).and.y(i)<=y(j).and.y(i)>=y(j))then;mask(i)=.true.;exit;end if
        end do
    end do
end subroutine
subroutine corner_indices(x,y,rw,index)
    real(dp),intent(in)::x(:),y(:),rw(4)
    integer,intent(out)::index(4)
    real(dp)::cx(4),cy(4),d,best
    integer::i,j
    if(size(x)/=size(y).or.size(x)==0)error stop 'corner_indices: invalid points'
    cx=[rw(1),rw(2),rw(2),rw(1)];cy=[rw(3),rw(3),rw(4),rw(4)]
    do j=1,4
        best=huge(1.0_dp);index(j)=1
        do i=1,size(x)
            d=(x(i)-cx(j))**2+(y(i)-cy(j))**2
            if(d<best)then;best=d;index(j)=i;end if
        end do
    end do
end subroutine
pure logical function midpoint_inside(x,y,rx,ry)
    real(dp),intent(in)::x(2),y(2),rx(2),ry(2)
    real(dp)::xm,ym
    xm=0.5_dp*(x(1)+x(2));ym=0.5_dp*(y(1)+y(2))
    midpoint_inside=rx(1)<xm.and.xm<rx(2).and.ry(1)<ym.and.ym<ry(2)
end function
subroutine find_new_in_old(xnew,ynew,xold,yold,index,tolerance)
    real(dp),intent(in)::xnew(:),ynew(:),xold(:),yold(:)
    integer,allocatable,intent(out)::index(:)
    real(dp),intent(in),optional::tolerance
    real(dp)::tol
    integer::i,j
    if(size(xnew)/=size(ynew).or.size(xold)/=size(yold))error stop 'find_new_in_old: length mismatch'
    tol=sqrt(epsilon(1.0_dp));if(present(tolerance))tol=tolerance
    allocate(index(size(xnew)));index=0
    do i=1,size(xnew)
        do j=1,size(xold)
            if(abs(xnew(i)-xold(j))<=tol.and.abs(ynew(i)-yold(j))<=tol)then;index(i)=j;exit;end if
        end do
    end do
end subroutine
end module deldir_utils
