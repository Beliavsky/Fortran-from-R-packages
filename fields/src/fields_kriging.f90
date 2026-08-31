! GPL-2.0-or-later. Universal kriging and spatial-process numerical core.
module fields_kriging
use fields_kinds, only: dp
use fields_linalg, only: solve_spd, solve_spd_mat, inverse_spd, logdet_spd, mvn_sample
use fields_covariance, only: stationary_covariance
use fields_polynomial, only: polynomial_basis, polynomial_power_table, polynomial_gradient, choose_int
implicit none
private

public :: krig_fit, krig_multi_fit
public :: krig_fit_covariance, krig_fit_stationary, krig_fit_stationary_gcv, krig_fit_stationary_reml
public :: krig_predict, krig_predict_se, krig_predict_covariance, krig_simulate, krig_predict_gradient_stationary
public :: krig_gcv, krig_profile_reml, krig_profile_mle
public :: krig_multi_fit_covariance, krig_multi_predict

type :: krig_fit
   integer :: n=0
   integer :: p=0
   real(dp) :: lambda=0.0_dp
   real(dp) :: sigma2=1.0_dp
   real(dp) :: tau2=0.0_dp
   real(dp) :: trace_a=0.0_dp
   real(dp) :: gcv=0.0_dp
   real(dp) :: profile_reml=0.0_dp
   real(dp) :: profile_mle=0.0_dp
   integer :: info=0
   real(dp), allocatable :: y(:), weights(:), k(:,:), t(:,:)
   real(dp), allocatable :: q(:,:), h_inv(:,:), beta(:), c(:), fitted(:), residuals(:)
end type krig_fit

type :: krig_multi_fit
   integer :: n=0
   integer :: p=0
   integer :: nresponse=0
   real(dp) :: lambda=0.0_dp
   integer :: info=0
   real(dp), allocatable :: weights(:), k(:,:), t(:,:), q(:,:), h_inv(:,:)
   real(dp), allocatable :: beta(:,:), c(:,:), fitted(:,:), residuals(:,:), sigma2(:), tau2(:)
end type krig_multi_fit

contains

function krig_fit_covariance(y,k,lambda,t,weights) result(fit)
real(dp), intent(in) :: y(:),k(:,:),lambda
real(dp), intent(in), optional :: t(:,:),weights(:)
type(krig_fit) :: fit
real(dp), allocatable :: tt(:,:),w(:),m(:,:),qt(:,:),rhs(:),a(:,:),tmp(:,:),r(:)
real(dp) :: rss,den
integer :: n,p,i,info1,info2
n=size(y)
if(size(k,1)/=n .or. size(k,2)/=n) error stop 'krig_fit_covariance: K dimension mismatch'
if(lambda<0.0_dp) error stop 'krig_fit_covariance: lambda must be nonnegative'
allocate(w(n)); w=1.0_dp
if(present(weights)) then
   if(size(weights)/=n .or. any(weights<=0.0_dp)) error stop 'krig_fit_covariance: invalid weights'
   w=weights
end if
if(present(t)) then
   if(size(t,1)/=n) error stop 'krig_fit_covariance: T dimension mismatch'
   tt=t
else
   allocate(tt(n,1)); tt=1.0_dp
end if
p=size(tt,2)
allocate(m(n,n)); m=k
do i=1,n; m(i,i)=m(i,i)+lambda/w(i); end do
fit%q=inverse_spd(m,info1)
if(info1/=0) then; fit%info=info1; return; end if
qt=matmul(fit%q,tt)
fit%h_inv=inverse_spd(matmul(transpose(tt),qt),info2)
if(info2/=0) then; fit%info=1000+info2; return; end if
rhs=matmul(transpose(tt),matmul(fit%q,y))
fit%beta=matmul(fit%h_inv,rhs)
r=y-matmul(tt,fit%beta)
fit%c=matmul(fit%q,r)
fit%fitted=matmul(tt,fit%beta)+matmul(k,fit%c)
fit%residuals=y-fit%fitted
rss=dot_product(r,fit%c)
fit%sigma2=rss/max(1.0_dp,real(n-p,dp))
fit%tau2=lambda*fit%sigma2
! Smoother matrix: KQ + (T-KQT)(T'QT)^-1 T'Q.
tmp=matmul(k,fit%q)
a=tmp+matmul(matmul(tt-matmul(k,qt),fit%h_inv),matmul(transpose(tt),fit%q))
fit%trace_a=0.0_dp
do i=1,n; fit%trace_a=fit%trace_a+a(i,i); end do
den=max(1.0e-14_dp,1.0_dp-fit%trace_a/real(n,dp))
fit%gcv=sum(w*fit%residuals**2)/real(n,dp)/den**2
fit%profile_reml=profile_reml_from_parts(m,tt,fit%q,r,p)
fit%profile_mle=profile_mle_from_parts(m,fit%q,r)
fit%n=n; fit%p=p; fit%lambda=lambda; fit%info=0
fit%y=y; fit%weights=w; fit%k=k; fit%t=tt
end function krig_fit_covariance

function krig_fit_stationary(x,y,lambda,model,a_range,smoothness,power,phi,t,weights,vmat) result(fit)
real(dp), intent(in) :: x(:,:),y(:),lambda
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: a_range,smoothness,power,phi,t(:,:),weights(:),vmat(:,:)
type(krig_fit) :: fit
real(dp), allocatable :: k(:,:)
k=stationary_covariance(x,x,model,a_range,smoothness,power,phi,vmat)
fit=krig_fit_covariance(y,k,lambda,t,weights)
end function krig_fit_stationary

function krig_fit_stationary_gcv(x,y,model,a_range,smoothness,power,phi,t,weights,vmat, &
                                 lambda_min,lambda_max,tol,maxiter) result(fit)
real(dp), intent(in) :: x(:,:),y(:)
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: a_range,smoothness,power,phi,t(:,:),weights(:),vmat(:,:),lambda_min,lambda_max,tol
integer, intent(in), optional :: maxiter
type(krig_fit) :: fit
real(dp), allocatable :: k(:,:)
real(dp) :: hbest
k=stationary_covariance(x,x,model,a_range,smoothness,power,phi,vmat)
hbest=optimize_log_lambda(y,k,t,weights,'gcv',lambda_min,lambda_max,tol,maxiter)
fit=krig_fit_covariance(y,k,exp(hbest),t,weights)
end function krig_fit_stationary_gcv

function krig_fit_stationary_reml(x,y,model,a_range,smoothness,power,phi,t,weights,vmat, &
                                  lambda_min,lambda_max,tol,maxiter) result(fit)
real(dp), intent(in) :: x(:,:),y(:)
character(len=*), intent(in), optional :: model
real(dp), intent(in), optional :: a_range,smoothness,power,phi,t(:,:),weights(:),vmat(:,:),lambda_min,lambda_max,tol
integer, intent(in), optional :: maxiter
type(krig_fit) :: fit
real(dp), allocatable :: k(:,:)
real(dp) :: hbest
k=stationary_covariance(x,x,model,a_range,smoothness,power,phi,vmat)
hbest=optimize_log_lambda(y,k,t,weights,'reml',lambda_min,lambda_max,tol,maxiter)
fit=krig_fit_covariance(y,k,exp(hbest),t,weights)
end function krig_fit_stationary_reml

real(dp) function optimize_log_lambda(y,k,t,weights,criterion,lambda_min,lambda_max,tol,maxiter) result(hbest)
real(dp), intent(in) :: y(:),k(:,:)
real(dp), intent(in), optional :: t(:,:),weights(:),lambda_min,lambda_max,tol
character(len=*), intent(in) :: criterion
integer, intent(in), optional :: maxiter
real(dp) :: a,b,c,d,fc,fd,gr,eps,lo,hi
integer :: it,nit
lo=1.0e-10_dp; if(present(lambda_min)) lo=lambda_min
hi=1.0e10_dp; if(present(lambda_max)) hi=lambda_max
if(lo<=0.0_dp .or. hi<=lo) error stop 'optimize_log_lambda: invalid bounds'
eps=1.0e-6_dp; if(present(tol)) eps=tol
nit=120; if(present(maxiter)) nit=maxiter
a=log(lo); b=log(hi); gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
c=b-gr*(b-a); d=a+gr*(b-a)
fc=krig_objective(c,y,k,t,weights,criterion); fd=krig_objective(d,y,k,t,weights,criterion)
do it=1,nit
   if(abs(b-a)<=eps*(1.0_dp+abs(a)+abs(b))) exit
   if(fc<=fd) then
      b=d; d=c; fd=fc; c=b-gr*(b-a); fc=krig_objective(c,y,k,t,weights,criterion)
   else
      a=c; c=d; fc=fd; d=a+gr*(b-a); fd=krig_objective(d,y,k,t,weights,criterion)
   end if
end do
if(fc<=fd) then; hbest=c; else; hbest=d; end if
end function optimize_log_lambda

real(dp) function krig_objective(h,y,k,t,weights,criterion) result(v)
real(dp), intent(in) :: h,y(:),k(:,:)
real(dp), intent(in), optional :: t(:,:),weights(:)
character(len=*), intent(in) :: criterion
type(krig_fit) :: fit
fit=krig_fit_covariance(y,k,exp(h),t,weights)
if(fit%info/=0) then; v=huge(1.0_dp); return; end if
select case(lower(trim(criterion)))
case('gcv'); v=fit%gcv
case('mle','ml'); v=fit%profile_mle
case default; v=fit%profile_reml
end select
end function krig_objective

real(dp) function krig_gcv(y,k,lambda,t,weights) result(v)
real(dp), intent(in) :: y(:),k(:,:),lambda
real(dp), intent(in), optional :: t(:,:),weights(:)
type(krig_fit) :: fit
fit=krig_fit_covariance(y,k,lambda,t,weights); v=fit%gcv
end function krig_gcv

real(dp) function krig_profile_reml(y,k,lambda,t,weights) result(v)
real(dp), intent(in) :: y(:),k(:,:),lambda
real(dp), intent(in), optional :: t(:,:),weights(:)
type(krig_fit) :: fit
fit=krig_fit_covariance(y,k,lambda,t,weights); v=fit%profile_reml
end function krig_profile_reml

real(dp) function krig_profile_mle(y,k,lambda,t,weights) result(v)
real(dp), intent(in) :: y(:),k(:,:),lambda
real(dp), intent(in), optional :: t(:,:),weights(:)
type(krig_fit) :: fit
fit=krig_fit_covariance(y,k,lambda,t,weights); v=fit%profile_mle
end function krig_profile_mle

function krig_predict(fit,k_new,t_new) result(pred)
type(krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:)
real(dp), intent(in), optional :: t_new(:,:)
real(dp), allocatable :: pred(:)
real(dp), allocatable :: tt(:,:)
integer :: m
if(size(k_new,2)/=fit%n) error stop 'krig_predict: K_new must be n_new by n_data'
m=size(k_new,1)
if(present(t_new)) then
   if(size(t_new,1)/=m .or. size(t_new,2)/=fit%p) error stop 'krig_predict: T_new dimension mismatch'
   tt=t_new
else
   allocate(tt(m,fit%p)); tt=0.0_dp
   if(fit%p==1) tt(:,1)=1.0_dp
end if
pred=matmul(tt,fit%beta)+matmul(k_new,fit%c)
end function krig_predict

function krig_predict_covariance(fit,k_new,k_newnew,t_new,include_nugget,new_weights) result(v)
type(krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:),k_newnew(:,:)
real(dp), intent(in), optional :: t_new(:,:),new_weights(:)
logical, intent(in), optional :: include_nugget
real(dp), allocatable :: v(:,:)
real(dp), allocatable :: tt(:,:),qk(:,:),u(:,:),tmp(:,:)
logical :: inc
integer :: m,i
m=size(k_new,1)
if(size(k_new,2)/=fit%n .or. size(k_newnew,1)/=m .or. size(k_newnew,2)/=m) &
   error stop 'krig_predict_covariance: dimension mismatch'
if(present(t_new)) then
   if(size(t_new,1)/=m .or. size(t_new,2)/=fit%p) error stop 'krig_predict_covariance: T_new dimension mismatch'
   tt=t_new
else
   allocate(tt(m,fit%p)); tt=0.0_dp; if(fit%p==1) tt(:,1)=1.0_dp
end if
qk=matmul(fit%q,transpose(k_new))
u=tt-matmul(k_new,matmul(fit%q,fit%t))
v=fit%sigma2*(k_newnew-matmul(k_new,qk)+matmul(matmul(u,fit%h_inv),transpose(u)))
! Symmetrize roundoff.
v=0.5_dp*(v+transpose(v))
inc=.false.; if(present(include_nugget)) inc=include_nugget
if(inc) then
   if(present(new_weights)) then
      if(size(new_weights)/=m .or. any(new_weights<=0.0_dp)) error stop 'krig_predict_covariance: invalid new weights'
      do i=1,m; v(i,i)=v(i,i)+fit%tau2/new_weights(i); end do
   else
      do i=1,m; v(i,i)=v(i,i)+fit%tau2; end do
   end if
end if
end function krig_predict_covariance

function krig_predict_se(fit,k_new,k_newnew,t_new,include_nugget,new_weights) result(se)
type(krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:),k_newnew(:,:)
real(dp), intent(in), optional :: t_new(:,:),new_weights(:)
logical, intent(in), optional :: include_nugget
real(dp), allocatable :: se(:),v(:,:)
integer :: i
v=krig_predict_covariance(fit,k_new,k_newnew,t_new,include_nugget,new_weights)
allocate(se(size(v,1)))
do i=1,size(v,1); se(i)=sqrt(max(0.0_dp,v(i,i))); end do
end function krig_predict_se

function krig_simulate(fit,k_new,k_newnew,nsim,t_new,include_nugget,new_weights,info) result(draws)
type(krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:),k_newnew(:,:)
integer, intent(in) :: nsim
real(dp), intent(in), optional :: t_new(:,:),new_weights(:)
logical, intent(in), optional :: include_nugget
integer, intent(out), optional :: info
real(dp), allocatable :: draws(:,:),mean(:),cov(:,:)
integer :: ierr
mean=krig_predict(fit,k_new,t_new)
cov=krig_predict_covariance(fit,k_new,k_newnew,t_new,include_nugget,new_weights)
draws=mvn_sample(mean,cov,nsim,ierr)
if(present(info)) info=ierr
end function krig_simulate

function krig_multi_fit_covariance(y,k,lambda,t,weights) result(fit)
real(dp), intent(in) :: y(:,:),k(:,:),lambda
real(dp), intent(in), optional :: t(:,:),weights(:)
type(krig_multi_fit) :: fit
type(krig_fit) :: one
integer :: j,n,r,p
n=size(y,1); r=size(y,2)
one=krig_fit_covariance(y(:,1),k,lambda,t,weights)
if(one%info/=0) then; fit%info=one%info; return; end if
p=one%p
fit%n=n; fit%nresponse=r; fit%p=p; fit%lambda=lambda; fit%info=0
fit%weights=one%weights; fit%k=k; fit%t=one%t; fit%q=one%q; fit%h_inv=one%h_inv
allocate(fit%beta(p,r),fit%c(n,r),fit%fitted(n,r),fit%residuals(n,r),fit%sigma2(r),fit%tau2(r))
fit%beta(:,1)=one%beta; fit%c(:,1)=one%c; fit%fitted(:,1)=one%fitted; fit%residuals(:,1)=one%residuals
fit%sigma2(1)=one%sigma2; fit%tau2(1)=one%tau2
do j=2,r
   call fit_multi_column(fit,y(:,j),j)
end do
end function krig_multi_fit_covariance

subroutine fit_multi_column(fit,y,j)
type(krig_multi_fit), intent(inout) :: fit
real(dp), intent(in) :: y(:)
integer, intent(in) :: j
real(dp), allocatable :: r(:),rhs(:)
rhs=matmul(transpose(fit%t),matmul(fit%q,y))
fit%beta(:,j)=matmul(fit%h_inv,rhs)
r=y-matmul(fit%t,fit%beta(:,j))
fit%c(:,j)=matmul(fit%q,r)
fit%fitted(:,j)=matmul(fit%t,fit%beta(:,j))+matmul(fit%k,fit%c(:,j))
fit%residuals(:,j)=y-fit%fitted(:,j)
fit%sigma2(j)=dot_product(r,fit%c(:,j))/max(1.0_dp,real(fit%n-fit%p,dp))
fit%tau2(j)=fit%lambda*fit%sigma2(j)
end subroutine fit_multi_column

function krig_multi_predict(fit,k_new,t_new) result(pred)
type(krig_multi_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:)
real(dp), intent(in), optional :: t_new(:,:)
real(dp), allocatable :: pred(:,:),tt(:,:)
integer :: m
m=size(k_new,1); if(size(k_new,2)/=fit%n) error stop 'krig_multi_predict: dimension mismatch'
if(present(t_new)) then
   tt=t_new
else
   allocate(tt(m,fit%p)); tt=0.0_dp; if(fit%p==1) tt(:,1)=1.0_dp
end if
pred=matmul(tt,fit%beta)+matmul(k_new,fit%c)
end function krig_multi_predict

real(dp) function profile_reml_from_parts(m,t,q,r,p) result(v)
real(dp), intent(in) :: m(:,:),t(:,:),q(:,:),r(:)
integer, intent(in) :: p
real(dp), allocatable :: h(:,:)
real(dp) :: rss,ldm,ldh,s2
integer :: n,info1,info2
n=size(r); rss=dot_product(r,matmul(q,r)); ldm=logdet_spd(m,info1)
h=matmul(transpose(t),matmul(q,t)); ldh=logdet_spd(h,info2)
if(info1/=0 .or. info2/=0 .or. rss<=0.0_dp .or. n<=p) then; v=huge(1.0_dp); return; end if
s2=rss/real(n-p,dp)
v=0.5_dp*(real(n-p,dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp))+log(s2))+ldm+ldh)
end function profile_reml_from_parts

real(dp) function profile_mle_from_parts(m,q,r) result(v)
real(dp), intent(in) :: m(:,:),q(:,:),r(:)
real(dp) :: rss,ldm,s2
integer :: n,info
n=size(r); rss=dot_product(r,matmul(q,r)); ldm=logdet_spd(m,info)
if(info/=0 .or. rss<=0.0_dp) then; v=huge(1.0_dp); return; end if
s2=rss/real(n,dp)
v=0.5_dp*(real(n,dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp))+log(s2))+ldm)
end function profile_mle_from_parts

pure function lower(s) result(t)
character(len=*),intent(in)::s
character(len=len(s))::t
integer::i,c
t=s
do i=1,len(s); c=iachar(t(i:i)); if(c>=65 .and. c<=90)t(i:i)=achar(c+32); end do
end function lower

function krig_predict_gradient_stationary(fit,xdata,xnew,model,a_range,smoothness,power,phi,m,step) result(g)
type(krig_fit),intent(in)::fit
real(dp),intent(in)::xdata(:,:),xnew(:,:)
character(len=*),intent(in),optional::model
real(dp),intent(in),optional::a_range,smoothness,power,phi,step
integer,intent(in),optional::m
real(dp),allocatable::g(:,:),xp(:,:),xm(:,:),kp(:,:),km(:,:),gp(:,:)
integer,allocatable::ptab(:,:)
integer::d,j,mm
real(dp)::h,scale
d=size(xnew,2);if(size(xdata,2)/=d .or. size(xdata,1)/=fit%n)error stop 'krig_predict_gradient_stationary: dimensions'
h=1.0e-5_dp;if(present(step))h=step
if(h<=0.0_dp)error stop 'krig_predict_gradient_stationary: step must be positive'
allocate(g(size(xnew,1),d));g=0.0_dp
! Polynomial/null-space derivative when fit%t has the standard fields polynomial basis.
mm=0;if(present(m))then
 mm=m
else
 do j=1,20
  if(choose_int(j+d-1,d)==fit%p)then;mm=j;exit;end if
 end do
end if
if(mm>0 .and. fit%p>0)then
 ptab=polynomial_power_table(d,mm)
 if(size(ptab,1)==fit%p)then
  gp=polynomial_gradient(xnew,fit%beta,ptab);g=g+gp
 end if
end if
! Numerical covariance derivative mirrors predictDerivative.Krig while supporting all stationary models.
do j=1,d
 xp=xnew;xm=xnew
 scale=max(1.0_dp,maxval(abs(xnew(:,j))))
 xp(:,j)=xp(:,j)+h*scale;xm(:,j)=xm(:,j)-h*scale
 kp=stationary_covariance(xp,xdata,model,a_range,smoothness,power,phi)
 km=stationary_covariance(xm,xdata,model,a_range,smoothness,power,phi)
 g(:,j)=g(:,j)+matmul((kp-km)/(2.0_dp*h*scale),fit%c)
end do
end function krig_predict_gradient_stationary

end module fields_kriging
