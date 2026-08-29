module spam_precision
use spam_kinds, only: dp
use spam_types, only: csr_matrix
use spam_csr, only: csr_identity, csr_diag, csr_transpose, csr_matmul, csr_add, csr_subtract, &
                    csr_kronecker, csr_diff, csr_toeplitz, csr_from_dense
implicit none
private
public :: precmat_rw1, precmat_rw2, precmat_rwn, precmat_season
public :: precmat_igmrf_reglat, precmat_igmrf_irreglat, precmat_gmrf_reglat
contains

function precmat_rw1(n) result(q)
integer,intent(in)::n
type(csr_matrix)::q
real(dp),allocatable::x(:,:)
integer::i
if(n<2) error stop 'precmat_rw1: n must be >= 2'
allocate(x(n,n))
x=0.0_dp
x(1,1)=1.0_dp
x(n,n)=1.0_dp
if(n>2) then
   do i=2,n-1
   x(i,i)=2.0_dp
   end do
end if
do i=1,n-1
x(i,i+1)=-1.0_dp
x(i+1,i)=-1.0_dp
end do
q=csr_from_dense(x)
end function precmat_rw1

function precmat_rw2(n) result(q)
integer,intent(in)::n
type(csr_matrix)::q
real(dp),allocatable::x(:,:)
integer::i
if(n<4) error stop 'precmat_rw2: n must be >= 4'
allocate(x(n,n))
x=0.0_dp
! D2' D2, D2 rows contain (1,-2,1).
do i=1,n-2
   x(i,i)=x(i,i)+1.0_dp
   x(i,i+1)=x(i,i+1)-2.0_dp
   x(i,i+2)=x(i,i+2)+1.0_dp
   x(i+1,i)=x(i+1,i)-2.0_dp
   x(i+1,i+1)=x(i+1,i+1)+4.0_dp
   x(i+1,i+2)=x(i+1,i+2)-2.0_dp
   x(i+2,i)=x(i+2,i)+1.0_dp
   x(i+2,i+1)=x(i+2,i+1)-2.0_dp
   x(i+2,i+2)=x(i+2,i+2)+1.0_dp
end do
q=csr_from_dense(x)
end function precmat_rw2

function precmat_rwn(n,order) result(q)
integer,intent(in)::n,order
type(csr_matrix)::q,d,dt
if(order<1 .or. order>=n) error stop 'precmat_rwn: order out of range'
d=csr_diff(csr_identity(n),lag=1,differences=order)
dt=csr_transpose(d)
q=csr_matmul(dt,d)
end function precmat_rwn

function precmat_season(n,season) result(q)
integer,intent(in)::n,season
type(csr_matrix)::q
real(dp),allocatable::x(:,:)
integer::r,i,j
if(season<1 .or. n<2*season) error stop 'precmat_season: require n >= 2*season'
! Upstream construction is B'B where every row of B contains a block
! of `season` consecutive ones. This reproduces the exact edge weights.
allocate(x(n,n))
x=0.0_dp
do r=1,n-season+1
   do i=r,r+season-1
      do j=r,r+season-1
         x(i,j)=x(i,j)+1.0_dp
      end do
   end do
end do
q=csr_from_dense(x)
end function precmat_season

function precmat_igmrf_reglat(n,m,order,anisotropy) result(q)
integer,intent(in)::n,m
integer,intent(in),optional::order
real(dp),intent(in),optional::anisotropy
type(csr_matrix)::q,a,b,c,d
integer::ord
real(dp)::an
if(n<2 .or. m<2) error stop 'precmat_igmrf_reglat: n,m must exceed 1'
ord=1
if(present(order))ord=order
an=1.0_dp
if(present(anisotropy))an=anisotropy
if(ord==1) then
   if(an<0.0_dp .or. an>2.0_dp) error stop 'precmat_igmrf_reglat: anisotropy not in [0,2]'
   a=precmat_rw1(m)
   b=csr_identity(n,2.0_dp-an)
   c=csr_identity(m,an)
   d=precmat_rw1(n)
else if(ord==2) then
   a=precmat_rw2(m)
   b=csr_identity(n)
   c=csr_identity(m)
   d=precmat_rw2(n)
else
   if(ord<1 .or. ord>=min(n,m)) error stop 'precmat_igmrf_reglat: order out of range'
   a=precmat_rwn(m,ord)
   b=csr_identity(n)
   c=csr_identity(m)
   d=precmat_rwn(n,ord)
end if
q=csr_add(csr_kronecker(a,b),csr_kronecker(c,d))
end function precmat_igmrf_reglat

function precmat_igmrf_irreglat(adj,eps) result(q)
type(csr_matrix),intent(in)::adj
real(dp),intent(in),optional::eps
type(csr_matrix)::q,a,d
real(dp),allocatable::deg(:)
integer::i,k
real(dp)::tol
if(adj%nrow/=adj%ncol) error stop 'precmat_igmrf_irreglat: adjacency must be square'
tol=epsilon(1.0_dp)*100.0_dp
if(present(eps))tol=eps*100.0_dp
a=adj
if(allocated(a%entries)) then
   where(abs(a%entries)>0.0_dp) a%entries=1.0_dp
end if
allocate(deg(a%nrow))
deg=0.0_dp
do i=1,a%nrow
   do k=a%rowpointers(i),a%rowpointers(i+1)-1
      deg(i)=deg(i)+1.0_dp
   end do
end do
d=csr_diag(deg)
q=csr_subtract(d,a,tol)
end function precmat_igmrf_irreglat

function precmat_gmrf_reglat(n,m,par,model,eps) result(q)
integer,intent(in)::n,m
real(dp),intent(in)::par(:)
character(len=*),intent(in),optional::model
real(dp),intent(in),optional::eps
type(csr_matrix)::q,p1,p2,p3,im,tx,ty,tz,tw
real(dp),allocatable::x(:),y(:),z(:),w(:),zero(:)
character(len=8)::modl
real(dp)::tol
if(n<2 .or. m<2) error stop 'precmat_gmrf_reglat: n,m must exceed 1'
modl='m1p1'
if(present(model))modl=trim(adjustl(model))
tol=0.0_dp
if(present(eps))tol=eps
allocate(x(n),y(max(n,n*m)),z(m),w(n),zero(m))
x=0
y=0
z=0
w=0
zero=0
im=csr_identity(m)
select case(trim(modl))
case('m1p1')
   if(size(par)<1) error stop 'precmat_gmrf_reglat: par too short'
   x(1:2)=[1.0_dp,-par(1)]
   y(n+1)=-par(1)
   p1=csr_kronecker(im,csr_toeplitz(x,eps=tol))
   p2=csr_toeplitz(y(:n*m),eps=tol)
   q=csr_add(p1,p2,tol)
case('m1p2')
   if(size(par)<2) error stop 'precmat_gmrf_reglat: par too short'
   x(1:2)=[1.0_dp,-par(1)]
   y(n+1)=-par(2)
   q=csr_add(csr_kronecker(im,csr_toeplitz(x,eps=tol)),csr_toeplitz(y(:n*m),eps=tol),tol)
case('m2p3')
   if(size(par)<3) error stop 'precmat_gmrf_reglat: par too short'
   x(1:2)=[1.0_dp,-par(1)]
   y(1:2)=[-par(2),-par(3)]
   z(2)=1.0_dp
   p1=csr_kronecker(im,csr_toeplitz(x,eps=tol))
   p2=csr_kronecker(csr_toeplitz(z,eps=tol),csr_toeplitz(y(:n),eps=tol))
   q=csr_add(p1,p2,tol)
case('m2p4')
   if(size(par)<4) error stop 'precmat_gmrf_reglat: par too short'
   x(1:2)=[1.0_dp,-par(1)]
   y(1:2)=[-par(2),-par(3)]
   w(1:2)=[-par(2),-par(4)]
   z(2)=1.0_dp
   p1=csr_kronecker(im,csr_toeplitz(x,eps=tol))
   p2=csr_kronecker(csr_toeplitz(z,zero,eps=tol),csr_toeplitz(y(:n),w,eps=tol))
   p3=csr_kronecker(csr_toeplitz(zero,z,eps=tol),csr_toeplitz(w,y(:n),eps=tol))
   q=csr_add(csr_add(p1,p2,tol),p3,tol)
case default
   error stop 'precmat_gmrf_reglat: unknown model'
end select
end function precmat_gmrf_reglat

end module spam_precision
