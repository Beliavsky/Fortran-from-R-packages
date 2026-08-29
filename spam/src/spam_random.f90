module spam_random
use spam_kinds, only: dp
use spam_types, only: csr_matrix, spam_chol
use spam_csr, only: csr_to_dense
use spam_cholesky, only: spam_chol_factor, spam_solve, spam_backsolve
use r_compat, only: rnorm_vec, rchisq
implicit none
private
public :: rmvnorm_cov, rmvnorm_prec, rmvnorm_canonical, rmvt_cov
public :: rmvnorm_cov_const, rmvnorm_prec_const, rmvnorm_canonical_const
public :: rmvnorm_conditional
interface
 subroutine dpotrf(uplo,n,a,lda,info)
  character(len=1) uplo
  integer n,lda,info
  double precision a(lda,*)
 end subroutine
 subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  character(len=1) side,uplo,transa,diag
  integer m,n,lda,ldb
  double precision alpha,a(lda,*),b(ldb,*)
 end subroutine
 subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  integer n,nrhs,lda,ldb,ipiv(*),info
  double precision a(lda,*),b(ldb,*)
 end subroutine
end interface
contains

function rmvnorm_cov(n,mu,sigma) result(x)
integer,intent(in)::n
real(dp),intent(in)::mu(:)
type(csr_matrix),intent(in)::sigma
real(dp),allocatable::x(:,:)
real(dp),allocatable::l(:,:),z(:)
integer::d,j,info
d=sigma%nrow
if(sigma%ncol/=d.or.size(mu)/=d)error stop 'rmvnorm_cov: dimension mismatch'
l=csr_to_dense(sigma)
call dpotrf('L',d,l,d,info)
if(info/=0)error stop 'rmvnorm_cov: covariance not SPD'
do j=1,d-1
l(j,j+1:d)=0.0_dp
end do
allocate(x(n,d),z(d))
do j=1,n
 z=rnorm_vec(d)
 x(j,:)=mu+matmul(l,z)
end do
end function rmvnorm_cov

function rmvnorm_prec(n,mu,q) result(x)
integer,intent(in)::n
real(dp),intent(in)::mu(:)
type(csr_matrix),intent(in)::q
real(dp),allocatable::x(:,:),z(:,:),s(:,:)
type(spam_chol)::f
integer::d,j
d=q%nrow
if(q%ncol/=d.or.size(mu)/=d)error stop 'rmvnorm_prec: dimension mismatch'
f=spam_chol_factor(q)
if(f%info/=0)error stop 'rmvnorm_prec: precision factorization failed'
allocate(z(d,n))
do j=1,n
z(:,j)=rnorm_vec(d)
end do
s=spam_backsolve(f,z)
allocate(x(n,d))
do j=1,n
x(j,:)=s(:,j)+mu
end do
end function rmvnorm_prec

function rmvnorm_canonical(n,b,q) result(x)
integer,intent(in)::n
real(dp),intent(in)::b(:)
type(csr_matrix),intent(in)::q
real(dp),allocatable::x(:,:),mu(:)
type(spam_chol)::f
f=spam_chol_factor(q)
if(f%info/=0)error stop 'rmvnorm_canonical: factorization failed'
mu=spam_solve(f,b)
x=rmvnorm_prec(n,mu,q)
end function rmvnorm_canonical

function rmvt_cov(n,sigma,df,delta,kshirsagar) result(x)
integer,intent(in)::n
type(csr_matrix),intent(in)::sigma
real(dp),intent(in)::df,delta(:)
logical,intent(in),optional::kshirsagar
real(dp),allocatable::x(:,:),z(:,:),cs(:)
logical::ks
integer::i,d
d=sigma%nrow
if(size(delta)/=d)error stop 'rmvt_cov: dimension mismatch'
if(df==0.0_dp .or. (df>0.0_dp .and. df>=huge(1.0_dp)/2))then
x=rmvnorm_cov(n,delta,sigma)
return
end if
if(df<=0.0_dp)error stop 'rmvt_cov: df must be positive'
ks=.false.
if(present(kshirsagar))ks=kshirsagar
z=rmvnorm_cov(n,merge(delta,spread(0.0_dp,1,d),ks),sigma)
cs=rchisq(n,df)
allocate(x(n,d))
do i=1,n
 if(ks)then
 x(i,:)=z(i,:)/sqrt(cs(i)/df)
 else
 x(i,:)=delta+z(i,:)/sqrt(cs(i)/df)
 end if
end do
end function rmvt_cov

function rmvnorm_cov_const(n,mu,sigma,a,aval) result(x)
integer,intent(in)::n
real(dp),intent(in)::mu(:),a(:,:),aval(:)
type(csr_matrix),intent(in)::sigma
real(dp),allocatable::x(:,:),s(:,:),u(:,:),w(:,:),rhs(:,:),corr(:,:),sd(:,:)
integer::d,k,i,info
integer,allocatable::ipiv(:)
d=sigma%nrow
k=size(a,1)
if(size(a,2)/=d.or.size(aval)/=k)error stop 'rmvnorm_cov_const: dimension mismatch'
x=rmvnorm_cov(n,mu,sigma)
sd=csr_to_dense(sigma)
u=matmul(sd,transpose(a))
w=matmul(a,u)
allocate(ipiv(k))
rhs=transpose(u)
call dgesv(k,d,w,k,ipiv,rhs,k,info)
if(info/=0) error stop 'rmvnorm_cov_const: singular constraint'
! rhs = W^-1 A Sigma, hence U = rhs^T
allocate(corr(k,n))
corr=matmul(a,transpose(x))-spread(aval,2,n)
x=x-transpose(matmul(transpose(rhs),corr))
end function rmvnorm_cov_const

function rmvnorm_prec_const(n,mu,q,a,aval) result(x)
integer,intent(in)::n
real(dp),intent(in)::mu(:),a(:,:),aval(:)
type(csr_matrix),intent(in)::q
real(dp),allocatable::x(:,:),v(:,:),w(:,:),rhs(:,:),corr(:,:)
type(spam_chol)::f
integer::d,k,info
integer,allocatable::ipiv(:)
d=q%nrow
k=size(a,1)
if(size(a,2)/=d.or.size(aval)/=k)error stop 'rmvnorm_prec_const: dimension mismatch'
f=spam_chol_factor(q)
if(f%info/=0)error stop 'rmvnorm_prec_const: factorization failed'
x=rmvnorm_prec(n,mu,q)
v=spam_solve(f,transpose(a))
w=matmul(a,v)
allocate(ipiv(k))
rhs=transpose(v)
call dgesv(k,d,w,k,ipiv,rhs,k,info)
if(info/=0) error stop 'rmvnorm_prec_const: singular constraint'
allocate(corr(k,n))
corr=matmul(a,transpose(x))-spread(aval,2,n)
x=x-transpose(matmul(transpose(rhs),corr))
end function rmvnorm_prec_const

function rmvnorm_canonical_const(n,b,q,a,aval) result(x)
integer,intent(in)::n
real(dp),intent(in)::b(:),a(:,:),aval(:)
type(csr_matrix),intent(in)::q
real(dp),allocatable::x(:,:),mu(:)
type(spam_chol)::f
f=spam_chol_factor(q)
if(f%info/=0)error stop 'rmvnorm_canonical_const: factorization failed'
mu=spam_solve(f,b)
x=rmvnorm_prec_const(n,mu,q,a,aval)
end function rmvnorm_canonical_const

function rmvnorm_conditional(n,y,mu,sxx,syy,sxy) result(x)
integer,intent(in)::n
real(dp),intent(in)::y(:),mu(:)
type(csr_matrix),intent(in)::sxx,syy,sxy
real(dp),allocatable::x(:,:),joint(:,:),sd(:,:),xy(:,:),yy(:,:),tmp(:,:),rhs(:,:)
type(csr_matrix)::sj
integer::nx,ny,i,info
integer,allocatable::ipiv(:)
nx=sxx%nrow
ny=syy%nrow
if(size(y)/=ny.or.size(mu)/=nx+ny)error stop 'rmvnorm_conditional: dimension mismatch'
! Use the exact conditional Gaussian identity: mu_x + Sigma_xy Sigma_yy^-1(y-mu_y),
! covariance Sigma_xx-Sigma_xy Sigma_yy^-1 Sigma_yx.
xy=csr_to_dense(sxy)
yy=csr_to_dense(syy)
allocate(ipiv(ny))
allocate(rhs(ny,nx))
rhs=transpose(xy)
call dgesv(ny,nx,yy,ny,ipiv,rhs,ny,info)
if(info/=0)error stop 'rmvnorm_conditional: SigmaYY singular'
tmp=csr_to_dense(sxx)-matmul(xy,rhs)
sj=dense_to_csr_local(tmp)
x=rmvnorm_cov(n,mu(:nx)+matmul(xy,solve_dense_vec(csr_to_dense(syy),y-mu(nx+1:))),sj)
end function rmvnorm_conditional

function dense_to_csr_local(d) result(a)
use spam_csr, only: csr_from_dense
real(dp),intent(in)::d(:,:)
type(csr_matrix)::a
a=csr_from_dense(d)
end function
function solve_dense_vec(a,b) result(x)
real(dp),intent(in)::a(:,:),b(:)
real(dp),allocatable::x(:)
real(dp),allocatable::aa(:,:),bb(:,:)
integer,allocatable::ipiv(:)
integer::n,info
n=size(b)
aa=a
allocate(bb(n,1))
bb(:,1)=b
allocate(ipiv(n))
call dgesv(n,1,aa,n,ipiv,bb,n,info)
if(info/=0)error stop 'solve_dense_vec: singular matrix'
allocate(x(n))
x=bb(:,1)
end function
end module spam_random
