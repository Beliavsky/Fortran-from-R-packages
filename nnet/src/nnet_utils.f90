! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module nnet_utils
use r_compat, only: dp
implicit none
private
public :: class_ind, which_is_max, which_is_max_rows, summarize_rows
contains

pure function class_ind(labels, n_classes) result(x)
integer, intent(in) :: labels(:)
integer, intent(in), optional :: n_classes
real(dp), allocatable :: x(:,:)
integer :: ncl, i
ncl=maxval(labels)
if(present(n_classes)) ncl=n_classes
if(any(labels<1) .or. any(labels>ncl)) error stop "class_ind: labels out of range"
allocate(x(size(labels),ncl))
x=0.0_dp
do i=1,size(labels)
x(i,labels(i))=1.0_dp
end do
end function class_ind

integer function which_is_max(x) result(idx)
real(dp), intent(in) :: x(:)
integer, allocatable :: ties(:)
real(dp) :: u
integer :: k
ties=pack([(k,k=1,size(x))],x==maxval(x))
if(size(ties)==1) then
   idx=ties(1)
else
   call random_number(u)
   idx=ties(min(size(ties),1+int(u*size(ties))))
end if
end function which_is_max

function which_is_max_rows(x) result(idx)
real(dp), intent(in) :: x(:,:)
integer, allocatable :: idx(:)
integer :: i
allocate(idx(size(x,1)))
do i=1,size(x,1)
idx(i)=which_is_max(x(i,:))
end do
end function which_is_max_rows

subroutine summarize_rows(x,y,x_out,y_out)
! Equivalent computational kernel to VR_summ2: lexicographically sort by X,
! then sum Y for rows with identical X.
real(dp), intent(in) :: x(:,:), y(:,:)
real(dp), allocatable, intent(out) :: x_out(:,:), y_out(:,:)
integer, allocatable :: ord(:)
real(dp), allocatable :: xs(:,:), ys(:,:), xt(:,:), yt(:,:)
integer :: i,nout
if(size(x,1)/=size(y,1)) error stop "summarize_rows: row counts differ"
allocate(ord(size(x,1)))
ord=[(i,i=1,size(x,1))]
call sort_indices_lex(x,ord,1,size(ord))
xs=x(ord,:)
ys=y(ord,:)
allocate(xt(size(x,1),size(x,2)),yt(size(y,1),size(y,2)))
nout=0
do i=1,size(x,1)
   if(i==1) then
      nout=nout+1
      xt(nout,:)=xs(i,:)
      yt(nout,:)=ys(i,:)
   else if(any(xs(i,:) /= xs(i-1,:))) then
      nout=nout+1
      xt(nout,:)=xs(i,:)
      yt(nout,:)=ys(i,:)
   else
      yt(nout,:)=yt(nout,:)+ys(i,:)
   end if
end do
allocate(x_out(nout,size(x,2)),y_out(nout,size(y,2)))
x_out=xt(1:nout,:)
y_out=yt(1:nout,:)
end subroutine summarize_rows

recursive subroutine sort_indices_lex(x,ord,lo,hi)
real(dp), intent(in) :: x(:,:)
integer, intent(inout) :: ord(:)
integer, intent(in) :: lo,hi
integer :: i,j,piv,tmp
if(lo>=hi) return
piv=ord((lo+hi)/2)
i=lo
j=hi
do
   do while(lex_less(x(ord(i),:),x(piv,:)))
   i=i+1
   end do
   do while(lex_less(x(piv,:),x(ord(j),:)))
   j=j-1
   end do
   if(i<=j) then
      tmp=ord(i)
      ord(i)=ord(j)
      ord(j)=tmp
      i=i+1
      j=j-1
   end if
   if(i>j) exit
end do
if(lo<j) call sort_indices_lex(x,ord,lo,j)
if(i<hi) call sort_indices_lex(x,ord,i,hi)
end subroutine sort_indices_lex

pure logical function lex_less(a,b) result(ans)
real(dp), intent(in) :: a(:),b(:)
integer :: k
ans=.false.
do k=1,size(a)
   if(a(k)<b(k)) then
   ans=.true.
   return
   else if(a(k)>b(k)) then
   return
   end if
end do
end function lex_less
end module nnet_utils
