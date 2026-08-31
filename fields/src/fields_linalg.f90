! GPL-2.0-or-later. See upstream/LICENSE.note and NOTICE.md.
module fields_linalg
use fields_kinds, only: dp
use la_lapack, only: dgesv => gesv, dpotrf => potrf, dpotri => potri, dpotrs => potrs, dsyev => syev
implicit none
private
public :: chol_upper, solve_spd, solve_spd_mat, inverse_spd, logdet_spd, &
          solve_general, symmetric_eigen, mvn_sample, identity_matrix

contains

function identity_matrix(n) result(a)
integer, intent(in) :: n
real(dp), allocatable :: a(:,:)
integer :: i
allocate(a(n,n)); a=0.0_dp
do i=1,n
   a(i,i)=1.0_dp
end do
end function identity_matrix

subroutine chol_upper(a,r,info)
real(dp), intent(in) :: a(:,:)
real(dp), allocatable, intent(out) :: r(:,:)
integer, intent(out) :: info
integer :: n,i
if(size(a,1)/=size(a,2)) error stop 'chol_upper: square matrix required'
n=size(a,1); allocate(r(n,n)); r=a
call dpotrf('U',n,r,n,info)
if(info==0) then
   do i=1,n
      if(i<n) r(i+1:n,i)=0.0_dp
   end do
end if
end subroutine chol_upper

function solve_spd(a,b,info) result(x)
real(dp), intent(in) :: a(:,:),b(:)
integer, intent(out), optional :: info
real(dp), allocatable :: x(:)
real(dp), allocatable :: r(:,:),bb(:,:)
integer :: ierr,n
n=size(a,1)
if(size(a,2)/=n .or. size(b)/=n) error stop 'solve_spd: dimension mismatch'
allocate(r(n,n)); r=a
call dpotrf('U',n,r,n,ierr)
allocate(bb(n,1)); bb(:,1)=b
if(ierr==0) call dpotrs('U',n,1,r,n,bb,n,ierr)
allocate(x(n))
if(ierr==0) then; x=bb(:,1); else; x=huge(1.0_dp); end if
if(present(info)) info=ierr
end function solve_spd

function solve_spd_mat(a,b,info) result(x)
real(dp), intent(in) :: a(:,:),b(:,:)
integer, intent(out), optional :: info
real(dp), allocatable :: x(:,:)
real(dp), allocatable :: r(:,:)
integer :: ierr,n,nrhs
n=size(a,1); nrhs=size(b,2)
if(size(a,2)/=n .or. size(b,1)/=n) error stop 'solve_spd_mat: dimension mismatch'
allocate(r(n,n),x(n,nrhs)); r=a; x=b
call dpotrf('U',n,r,n,ierr)
if(ierr==0) call dpotrs('U',n,nrhs,r,n,x,n,ierr)
if(ierr/=0) x=huge(1.0_dp)
if(present(info)) info=ierr
end function solve_spd_mat

function inverse_spd(a,info) result(ai)
real(dp), intent(in) :: a(:,:)
integer, intent(out), optional :: info
real(dp), allocatable :: ai(:,:)
integer :: ierr,n,i,j
n=size(a,1)
if(size(a,2)/=n) error stop 'inverse_spd: square matrix required'
allocate(ai(n,n)); ai=a
call dpotrf('U',n,ai,n,ierr)
if(ierr==0) call dpotri('U',n,ai,n,ierr)
if(ierr==0) then
   do j=1,n
      do i=j+1,n
         ai(i,j)=ai(j,i)
      end do
   end do
else
   ai=huge(1.0_dp)
end if
if(present(info)) info=ierr
end function inverse_spd

real(dp) function logdet_spd(a,info) result(v)
real(dp), intent(in) :: a(:,:)
integer, intent(out), optional :: info
real(dp), allocatable :: r(:,:)
integer :: ierr,n,i
n=size(a,1); allocate(r(n,n)); r=a
call dpotrf('U',n,r,n,ierr)
if(ierr/=0) then
   v=-huge(1.0_dp)
else
   v=0.0_dp
   do i=1,n
      v=v+2.0_dp*log(r(i,i))
   end do
end if
if(present(info)) info=ierr
end function logdet_spd

function solve_general(a,b,info) result(x)
real(dp), intent(in) :: a(:,:),b(:,:)
integer, intent(out), optional :: info
real(dp), allocatable :: x(:,:),aa(:,:)
integer, allocatable :: ipiv(:)
integer :: ierr,n,nrhs
n=size(a,1); nrhs=size(b,2)
if(size(a,2)/=n .or. size(b,1)/=n) error stop 'solve_general: dimension mismatch'
allocate(aa(n,n),x(n,nrhs),ipiv(n)); aa=a; x=b
call dgesv(n,nrhs,aa,n,ipiv,x,n,ierr)
if(ierr/=0) x=huge(1.0_dp)
if(present(info)) info=ierr
end function solve_general

subroutine symmetric_eigen(a,values,vectors,info)
real(dp), intent(in) :: a(:,:)
real(dp), allocatable, intent(out) :: values(:),vectors(:,:)
integer, intent(out) :: info
real(dp), allocatable :: work(:),tmp(:,:)
real(dp) :: qwork(1)
integer :: n,lwork
n=size(a,1)
if(size(a,2)/=n) error stop 'symmetric_eigen: square matrix required'
allocate(tmp(n,n),values(n)); tmp=a
lwork=-1
call dsyev('V','U',n,tmp,n,values,qwork,lwork,info)
if(info/=0) then; allocate(vectors(0,0)); return; end if
lwork=max(1,int(qwork(1))); allocate(work(lwork))
call dsyev('V','U',n,tmp,n,values,work,lwork,info)
allocate(vectors(n,n)); vectors=tmp
end subroutine symmetric_eigen

function mvn_sample(mean,cov,nsim,info) result(x)
use r_mod, only: rnorm1
real(dp), intent(in) :: mean(:),cov(:,:)
integer, intent(in) :: nsim
integer, intent(out), optional :: info
real(dp), allocatable :: x(:,:)
real(dp), allocatable :: r(:,:),z(:)
integer :: ierr,n,i,j
n=size(mean)
if(size(cov,1)/=n .or. size(cov,2)/=n) error stop 'mvn_sample: dimension mismatch'
call chol_upper(cov,r,ierr)
allocate(x(n,nsim),z(n))
if(ierr/=0) then
   x=huge(1.0_dp)
else
   do j=1,nsim
      do i=1,n; z(i)=rnorm1(); end do
      x(:,j)=mean+matmul(transpose(r),z)
   end do
end if
if(present(info)) info=ierr
end function mvn_sample

end module fields_linalg
