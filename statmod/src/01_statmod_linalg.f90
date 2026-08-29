module statmod_linalg
use, intrinsic :: iso_fortran_env, only: real64
use r_linalg, only: shared_cholesky_factor => cholesky_factor
use r_linalg, only: shared_full_svd => full_svd
use r_linalg, only: shared_least_squares_svd => least_squares_svd
use r_linalg, only: shared_solve_spd => solve_spd
use r_linalg, only: shared_spd_inverse_logdet => spd_inverse_logdet
implicit none
private
integer, parameter, public :: dp = real64
public :: least_squares, weighted_least_squares, symmetric_inverse, solve_spd
public :: svd_full_u, column_space_basis, weighted_hat_basis, logdet_xtwx
public :: mean_vec, variance_vec

contains

subroutine least_squares(x,y,beta,fitted,resid,rank,info)
real(dp), intent(in) :: x(:,:), y(:)
real(dp), allocatable, intent(out) :: beta(:), fitted(:), resid(:)
integer, intent(out) :: rank, info
integer :: m,n
m=size(x,1)
n=size(x,2)
allocate(beta(n),fitted(m),resid(m))
call shared_least_squares_svd(x,y,beta,rank,info,rcond=sqrt(epsilon(1.0_dp)))
if(info/=0) then
beta=0
fitted=0
resid=y
return
end if
fitted=matmul(x,beta)
resid=y-fitted
end subroutine least_squares

subroutine weighted_least_squares(x,y,w,beta,fitted,resid,rank,info)
real(dp), intent(in) :: x(:,:),y(:),w(:)
real(dp), allocatable, intent(out) :: beta(:),fitted(:),resid(:)
integer, intent(out) :: rank,info
real(dp), allocatable :: xs(:,:),ys(:),fitw(:),resw(:)
integer :: i
allocate(xs(size(x,1),size(x,2)),ys(size(y)))
do i=1,size(y)
   xs(i,:)=sqrt(max(w(i),0.0_dp))*x(i,:)
   ys(i)=sqrt(max(w(i),0.0_dp))*y(i)
end do
call least_squares(xs,ys,beta,fitw,resw,rank,info)
allocate(fitted(size(y)),resid(size(y)))
fitted=matmul(x,beta)
resid=y-fitted
end subroutine weighted_least_squares

subroutine solve_spd(a,b,x,info)
real(dp), intent(in) :: a(:,:),b(:)
real(dp), allocatable, intent(out) :: x(:)
integer, intent(out) :: info
integer :: n
n=size(b)
allocate(x(n))
call shared_solve_spd(a,b,x,info,upper=.true.)
if(info/=0) x=0.0_dp
end subroutine solve_spd

subroutine symmetric_inverse(a,ainv,info)
real(dp), intent(in) :: a(:,:)
real(dp), allocatable, intent(out) :: ainv(:,:)
integer, intent(out) :: info
real(dp) :: unused_logdet
call shared_spd_inverse_logdet(a,ainv,unused_logdet,info)
end subroutine symmetric_inverse

subroutine svd_full_u(a,u,s,rank,info)
real(dp), intent(in) :: a(:,:)
real(dp), allocatable, intent(out) :: u(:,:),s(:)
integer, intent(out) :: rank,info
integer :: m,n,k
real(dp), allocatable :: vt(:,:)
real(dp) :: tol
m=size(a,1)
n=size(a,2)
k=min(m,n)
call shared_full_svd(a,u,s,vt,info)
if(info/=0) then
rank=0
return
end if
if(k==0) then
   rank=0
else
   tol=max(m,n)*epsilon(1.0_dp)*maxval(s)
   rank=count(s>tol)
end if
end subroutine svd_full_u

subroutine column_space_basis(x,q,rank,info)
real(dp), intent(in) :: x(:,:)
real(dp), allocatable, intent(out) :: q(:,:)
integer, intent(out) :: rank,info
real(dp), allocatable :: u(:,:),s(:)
call svd_full_u(x,u,s,rank,info)
if(info/=0) then
allocate(q(size(x,1),0))
return
end if
allocate(q(size(x,1),rank))
if(rank>0) q=u(:,1:rank)
end subroutine column_space_basis

subroutine weighted_hat_basis(x,w,q,rank,info)
real(dp), intent(in) :: x(:,:),w(:)
real(dp), allocatable, intent(out) :: q(:,:)
integer, intent(out) :: rank,info
real(dp), allocatable :: xs(:,:)
integer :: i
allocate(xs(size(x,1),size(x,2)))
do i=1,size(x,1)
   xs(i,:)=sqrt(max(w(i),0.0_dp))*x(i,:)
end do
call column_space_basis(xs,q,rank,info)
end subroutine weighted_hat_basis

function logdet_xtwx(x,w,info) result(val)
real(dp), intent(in) :: x(:,:),w(:)
integer, intent(out) :: info
real(dp) :: val
real(dp), allocatable :: a(:,:),factor(:,:)
integer :: i,j,n
n=size(x,2)
allocate(a(n,n))
a=0.0_dp
do i=1,size(x,1)
   do j=1,n
      a(j,:)=a(j,:)+w(i)*x(i,j)*x(i,:)
   end do
end do
call shared_cholesky_factor(a,factor,info,upper=.true.)
if(info/=0) then
   val=huge(1.0_dp)
else
   val=0.0_dp
   do i=1,n
      val=val+2.0_dp*log(abs(factor(i,i)))
   end do
end if
end function logdet_xtwx

pure function mean_vec(x) result(v)
real(dp),intent(in)::x(:)
real(dp)::v
if(size(x)==0) then
v=0.0_dp
else
v=sum(x)/real(size(x),dp)
end if
end function mean_vec

pure function variance_vec(x) result(v)
real(dp),intent(in)::x(:)
real(dp)::v,m
if(size(x)<=1) then
v=0.0_dp
return
end if
m=sum(x)/real(size(x),dp)
v=sum((x-m)**2)/real(size(x)-1,dp)
end function variance_vec

end module statmod_linalg
