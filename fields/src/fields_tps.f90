! GPL-2.0-or-later. Thin-plate spline numerical kernels translated from fields.
module fields_tps
use fields_kinds, only: dp
use fields_covariance, only: radial_covariance
use fields_polynomial, only: polynomial_basis, choose_int
use fields_linalg, only: solve_general, identity_matrix
implicit none
private

public :: tps_fit_type, tps_covariance, tps_covariance_marginal
public :: tps_fit, tps_fit_gcv, tps_predict, tps_default_m

type :: tps_fit_type
   integer :: n=0
   integer :: dimension=0
   integer :: m=2
   real(dp) :: lambda=0.0_dp
   real(dp) :: trace_a=0.0_dp
   real(dp) :: gcv=0.0_dp
   integer :: info=0
   real(dp), allocatable :: x(:,:), y(:), weights(:), k(:,:), t(:,:)
   real(dp), allocatable :: c(:), beta(:), fitted(:), residuals(:), diag_a(:)
end type tps_fit_type

contains

pure integer function tps_default_m(d) result(m)
integer, intent(in) :: d
m=max(2,int(ceiling(0.5_dp*real(d,dp)+0.1_dp)))
end function tps_default_m

function tps_covariance(x1,x2,cardinal_x,m) result(cov)
real(dp), intent(in) :: x1(:,:),x2(:,:),cardinal_x(:,:)
integer, intent(in), optional :: m
real(dp), allocatable :: cov(:,:)
real(dp), allocatable :: pcard(:,:),pcoef(:,:),p1(:,:),p2(:,:),e12(:,:),e1c(:,:),e2c(:,:),ecc(:,:)
integer :: mm,info
mm=tps_default_m(size(x1,2)); if(present(m)) mm=m
if(size(x2,2)/=size(x1,2) .or. size(cardinal_x,2)/=size(x1,2)) error stop 'tps_covariance: dimension mismatch'
call polynomial_basis(cardinal_x,mm,pcard,info=info)
if(info/=0 .or. size(pcard,1)/=size(pcard,2)) error stop 'tps_covariance: cardinal_x must give a square nonsingular polynomial basis'
pcoef=solve_general(pcard,identity_matrix(size(pcard,1)),info)
if(info/=0) error stop 'tps_covariance: singular cardinal polynomial matrix'
call polynomial_basis(x1,mm,p1,info=info); p1=matmul(p1,pcoef)
call polynomial_basis(x2,mm,p2,info=info); p2=matmul(p2,pcoef)
e12=radial_covariance(x1,x2,mm)
e1c=radial_covariance(x1,cardinal_x,mm)
e2c=radial_covariance(x2,cardinal_x,mm)
ecc=radial_covariance(cardinal_x,cardinal_x,mm)
cov=e12-matmul(e1c,transpose(p2))-matmul(p1,transpose(e2c))+ &
    matmul(matmul(p1,ecc),transpose(p2))+matmul(p1,transpose(p2))
end function tps_covariance

function tps_covariance_marginal(x,cardinal_x,m) result(v)
real(dp), intent(in) :: x(:,:),cardinal_x(:,:)
integer, intent(in), optional :: m
real(dp), allocatable :: v(:),p(:,:),pcard(:,:),pcoef(:,:),e(:,:),ecc(:,:)
integer :: mm,info,i
mm=tps_default_m(size(x,2)); if(present(m)) mm=m
call polynomial_basis(cardinal_x,mm,pcard,info=info)
if(size(pcard,1)/=size(pcard,2)) error stop 'tps_covariance_marginal: invalid cardinal set'
pcoef=solve_general(pcard,identity_matrix(size(pcard,1)),info)
call polynomial_basis(x,mm,p,info=info); p=matmul(p,pcoef)
e=radial_covariance(x,cardinal_x,mm); ecc=radial_covariance(cardinal_x,cardinal_x,mm)
allocate(v(size(x,1)))
do i=1,size(x,1)
   v(i)=-2.0_dp*dot_product(e(i,:),p(i,:))+dot_product(p(i,:),matmul(ecc,p(i,:)))+dot_product(p(i,:),p(i,:))
end do
end function tps_covariance_marginal

function tps_fit(x,y,lambda,m,weights) result(fit)
real(dp), intent(in) :: x(:,:),y(:),lambda
integer, intent(in), optional :: m
real(dp), intent(in), optional :: weights(:)
type(tps_fit_type) :: fit
real(dp), allocatable :: k(:,:),t(:,:),w(:),aug(:,:),rhs(:,:),sol(:,:),g(:,:),a(:,:),design(:,:)
integer :: n,p,mm,i,info
real(dp) :: den
n=size(x,1); if(size(y)/=n .or. lambda<0.0_dp) error stop 'tps_fit: invalid input'
mm=tps_default_m(size(x,2)); if(present(m)) mm=m
if(2*mm-size(x,2)<=0) error stop 'tps_fit: m too small for thin-plate spline'
call polynomial_basis(x,mm,t,info=info); if(info/=0) error stop 'tps_fit: polynomial basis failed'
p=size(t,2); if(n<=p) error stop 'tps_fit: need more observations than null-space terms'
k=radial_covariance(x,x,mm)
allocate(w(n)); w=1.0_dp
if(present(weights)) then
   if(size(weights)/=n .or. any(weights<=0.0_dp)) error stop 'tps_fit: invalid weights'
   w=weights
end if
allocate(aug(n+p,n+p)); aug=0.0_dp
aug(:n,:n)=k
do i=1,n; aug(i,i)=aug(i,i)+lambda/w(i); end do
aug(:n,n+1:n+p)=t; aug(n+1:n+p,:n)=transpose(t)
allocate(rhs(n+p,1)); rhs=0.0_dp; rhs(:n,1)=y
sol=solve_general(aug,rhs,info)
if(info/=0) then; fit%info=info; return; end if
fit%c=sol(:n,1); fit%beta=sol(n+1:n+p,1)
fit%fitted=matmul(k,fit%c)+matmul(t,fit%beta); fit%residuals=y-fit%fitted
! The first n columns of inv(aug) map y to coefficients. This yields the exact smoother trace.
g=solve_general(aug,identity_matrix(n+p),info)
allocate(design(n,n+p)); design(:,:n)=k; design(:,n+1:n+p)=t
a=matmul(design,g(:,1:n))
fit%trace_a=0.0_dp
allocate(fit%diag_a(n))
do i=1,n
   fit%diag_a(i)=a(i,i)
   fit%trace_a=fit%trace_a+a(i,i)
end do
den=max(1.0e-14_dp,1.0_dp-fit%trace_a/real(n,dp))
fit%gcv=sum(w*fit%residuals**2)/real(n,dp)/den**2
fit%n=n; fit%dimension=size(x,2); fit%m=mm; fit%lambda=lambda; fit%info=0
fit%x=x; fit%y=y; fit%weights=w; fit%k=k; fit%t=t
end function tps_fit

function tps_fit_gcv(x,y,m,weights,lambda_min,lambda_max,tol,maxiter) result(fit)
real(dp), intent(in) :: x(:,:),y(:)
integer, intent(in), optional :: m,maxiter
real(dp), intent(in), optional :: weights(:),lambda_min,lambda_max,tol
type(tps_fit_type) :: fit
real(dp) :: lo,hi,a,b,c,d,fc,fd,gr,eps
integer :: it,nit,mm
lo=1.0e-10_dp; if(present(lambda_min)) lo=lambda_min
hi=1.0e10_dp; if(present(lambda_max)) hi=lambda_max
eps=1.0e-6_dp; if(present(tol)) eps=tol
nit=100; if(present(maxiter)) nit=maxiter
mm=tps_default_m(size(x,2)); if(present(m)) mm=m
a=log(lo); b=log(hi); gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
c=b-gr*(b-a); d=a+gr*(b-a); fc=tps_gcv_at(c,x,y,mm,weights); fd=tps_gcv_at(d,x,y,mm,weights)
do it=1,nit
   if(abs(b-a)<=eps*(1.0_dp+abs(a)+abs(b))) exit
   if(fc<=fd) then
      b=d; d=c; fd=fc; c=b-gr*(b-a); fc=tps_gcv_at(c,x,y,mm,weights)
   else
      a=c; c=d; fc=fd; d=a+gr*(b-a); fd=tps_gcv_at(d,x,y,mm,weights)
   end if
end do
if(fc<=fd) then; fit=tps_fit(x,y,exp(c),mm,weights); else; fit=tps_fit(x,y,exp(d),mm,weights); end if
end function tps_fit_gcv

real(dp) function tps_gcv_at(h,x,y,m,weights) result(v)
real(dp), intent(in) :: h,x(:,:),y(:)
integer, intent(in) :: m
real(dp), intent(in), optional :: weights(:)
type(tps_fit_type) :: f
f=tps_fit(x,y,exp(h),m,weights)
if(f%info/=0) then; v=huge(1.0_dp); else; v=f%gcv; end if
end function tps_gcv_at

function tps_predict(fit,xnew) result(pred)
type(tps_fit_type), intent(in) :: fit
real(dp), intent(in) :: xnew(:,:)
real(dp), allocatable :: pred(:),knew(:,:),tnew(:,:)
integer :: info
if(size(xnew,2)/=fit%dimension) error stop 'tps_predict: coordinate dimension mismatch'
knew=radial_covariance(xnew,fit%x,fit%m)
call polynomial_basis(xnew,fit%m,tnew,info=info)
if(info/=0) error stop 'tps_predict: polynomial basis failed'
pred=matmul(knew,fit%c)+matmul(tnew,fit%beta)
end function tps_predict

end module fields_tps
