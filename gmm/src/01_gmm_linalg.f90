! SPDX-License-Identifier: GPL-2.0-or-later
module gmm_linalg
use r_compat, only: dp
use r_linalg, only: shared_inverse_matrix => inverse_matrix
use r_linalg, only: shared_least_squares => least_squares
use r_linalg, only: shared_solve_system => solve_system
implicit none
private
public :: solve_linear, invert_matrix, least_squares, quadratic_inverse
public :: colmeans_mat, center_columns, symmetrize, kron, block_diag

contains

subroutine solve_linear(a,b,x,info)
real(dp),intent(in)::a(:,:),b(:)
real(dp),allocatable,intent(out)::x(:)
integer,intent(out)::info
integer::n
n=size(a,1)
allocate(x(n))
call shared_solve_system(a,b,x,info)
if(info/=0) x=0.0_dp
end subroutine solve_linear

subroutine invert_matrix(a,ainv,info)
real(dp),intent(in)::a(:,:)
real(dp),allocatable,intent(out)::ainv(:,:)
integer,intent(out)::info
call shared_inverse_matrix(a,ainv,info)
if(info/=0) ainv=0.0_dp
end subroutine invert_matrix

subroutine least_squares(a,b,x,info)
real(dp),intent(in)::a(:,:),b(:)
real(dp),allocatable,intent(out)::x(:)
integer,intent(out)::info
integer::n
n=size(a,2)
allocate(x(n))
call shared_least_squares(a,b,x,info)
if(info/=0) x=0.0_dp
end subroutine least_squares

function quadratic_inverse(a,x,info) result(v)
real(dp),intent(in)::a(:,:),x(:)
integer,intent(out)::info
real(dp)::v
real(dp),allocatable::z(:)
call solve_linear(a,x,z,info)
if(info==0) then
v=dot_product(x,z)
else
v=huge(1.0_dp)
end if
end function quadratic_inverse

pure function colmeans_mat(x) result(m)
real(dp),intent(in)::x(:,:)
real(dp)::m(size(x,2))
if(size(x,1)==0) then
m=0.0_dp
else
m=sum(x,dim=1)/real(size(x,1),dp)
end if
end function colmeans_mat

pure function center_columns(x) result(y)
real(dp),intent(in)::x(:,:)
real(dp)::y(size(x,1),size(x,2)),m(size(x,2))
m=colmeans_mat(x)
y=x-spread(m,1,size(x,1))
end function center_columns

pure function symmetrize(a) result(b)
real(dp),intent(in)::a(:,:)
real(dp)::b(size(a,1),size(a,2))
b=0.5_dp*(a+transpose(a))
end function symmetrize

pure function kron(a,b) result(c)
real(dp),intent(in)::a(:,:),b(:,:)
real(dp)::c(size(a,1)*size(b,1),size(a,2)*size(b,2))
integer::i,j,rb,cb
rb=size(b,1)
cb=size(b,2)
do j=1,size(a,2)
   do i=1,size(a,1)
      c((i-1)*rb+1:i*rb,(j-1)*cb+1:j*cb)=a(i,j)*b
   end do
end do
end function kron

subroutine block_diag(blocks,c)
real(dp),intent(in)::blocks(:,:,:)
real(dp),allocatable,intent(out)::c(:,:)
integer::nb,r,cc,i,r0,c0
! blocks(:,:,i) must have common rectangular dimensions.
r=size(blocks,1)
cc=size(blocks,2)
nb=size(blocks,3)
allocate(c(r*nb,cc*nb))
c=0.0_dp
r0=0
c0=0
do i=1,nb
   c(r0+1:r0+r,c0+1:c0+cc)=blocks(:,:,i)
   r0=r0+r
   c0=c0+cc
end do
end subroutine block_diag

end module gmm_linalg
