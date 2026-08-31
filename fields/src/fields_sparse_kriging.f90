! GPL-2.0-or-later. Sparse Kriging/mKrig numerical path using the translated spam package.
module fields_sparse_kriging
use fields_kinds, only: dp
use fields_linalg, only: inverse_spd, logdet_spd
use fields_covariance, only: wendland_covariance_sparse, stationary_taper_sparse
use spam_types, only: csr_matrix, spam_chol
use spam_csr, only: csr_add, csr_diag, csr_matvec, csr_matmat
use spam_cholesky, only: spam_chol_factor, spam_solve, spam_logdet
implicit none
private

public :: sparse_krig_fit, sparse_krig_multi_fit
public :: sparse_krig_fit_covariance, sparse_krig_fit_wendland, sparse_krig_fit_taper
public :: sparse_krig_predict, sparse_krig_predict_covariance, sparse_krig_predict_se
public :: sparse_krig_multi_fit_covariance, sparse_krig_multi_predict

type :: sparse_krig_fit
   integer :: n=0,p=0,info=0
   real(dp) :: lambda=0.0_dp,sigma2=1.0_dp,tau2=0.0_dp,profile_reml=0.0_dp
   type(csr_matrix) :: k
   type(spam_chol) :: factor
   real(dp), allocatable :: y(:),weights(:),t(:,:),qt(:,:),h_inv(:,:),beta(:),c(:),fitted(:),residuals(:)
end type sparse_krig_fit

type :: sparse_krig_multi_fit
   integer :: n=0,p=0,nresponse=0,info=0
   real(dp) :: lambda=0.0_dp
   type(csr_matrix) :: k
   type(spam_chol) :: factor
   real(dp), allocatable :: weights(:),t(:,:),qt(:,:),h_inv(:,:)
   real(dp), allocatable :: beta(:,:),c(:,:),fitted(:,:),residuals(:,:),sigma2(:),tau2(:)
end type sparse_krig_multi_fit

contains

function sparse_krig_fit_covariance(y,k,lambda,t,weights,pivot) result(fit)
real(dp), intent(in) :: y(:),lambda
type(csr_matrix), intent(in) :: k
real(dp), intent(in), optional :: t(:,:),weights(:)
character(len=*), intent(in), optional :: pivot
type(sparse_krig_fit) :: fit
real(dp), allocatable :: w(:),tt(:,:),d(:),qy(:),h(:,:),rhs(:),r(:)
type(csr_matrix) :: md
integer :: n,p,info
real(dp) :: rss,ldh
n=size(y)
if(k%nrow/=n .or. k%ncol/=n .or. lambda<0.0_dp) error stop 'sparse_krig_fit_covariance: invalid dimensions/lambda'
allocate(w(n)); w=1.0_dp
if(present(weights)) then
   if(size(weights)/=n .or. any(weights<=0.0_dp)) error stop 'sparse_krig_fit_covariance: invalid weights'
   w=weights
end if
if(present(t)) then
   if(size(t,1)/=n) error stop 'sparse_krig_fit_covariance: T dimension mismatch'
   tt=t
else
   allocate(tt(n,1)); tt=1.0_dp
end if
p=size(tt,2); allocate(d(n)); d=lambda/w
md=csr_add(k,csr_diag(d))
if(present(pivot)) then; fit%factor=spam_chol_factor(md,pivot=pivot); else; fit%factor=spam_chol_factor(md); end if
if(fit%factor%info/=0) then; fit%info=fit%factor%info; return; end if
qy=spam_solve(fit%factor,y); fit%qt=spam_solve(fit%factor,tt)
h=matmul(transpose(tt),fit%qt); fit%h_inv=inverse_spd(h,info)
if(info/=0) then; fit%info=1000+info; return; end if
rhs=matmul(transpose(tt),qy); fit%beta=matmul(fit%h_inv,rhs)
r=y-matmul(tt,fit%beta); fit%c=spam_solve(fit%factor,r)
fit%fitted=matmul(tt,fit%beta)+csr_matvec(k,fit%c); fit%residuals=y-fit%fitted
rss=dot_product(r,fit%c); fit%sigma2=rss/max(1.0_dp,real(n-p,dp)); fit%tau2=lambda*fit%sigma2
ldh=logdet_spd(h,info)
if(info==0 .and. rss>0.0_dp .and. n>p) then
   fit%profile_reml=0.5_dp*(real(n-p,dp)*(1.0_dp+log(2.0_dp*acos(-1.0_dp))+ &
      log(rss/real(n-p,dp)))+spam_logdet(fit%factor)+ldh)
else
   fit%profile_reml=huge(1.0_dp)
end if
fit%n=n; fit%p=p; fit%lambda=lambda; fit%info=0; fit%k=k; fit%y=y; fit%weights=w; fit%t=tt
end function sparse_krig_fit_covariance

function sparse_krig_fit_wendland(x,y,lambda,a_range,korder,t,weights,pivot) result(fit)
real(dp), intent(in) :: x(:,:),y(:),lambda,a_range
integer, intent(in), optional :: korder
real(dp), intent(in), optional :: t(:,:),weights(:)
character(len=*), intent(in), optional :: pivot
type(sparse_krig_fit) :: fit
type(csr_matrix) :: k
k=wendland_covariance_sparse(x,x,a_range,korder)
fit=sparse_krig_fit_covariance(y,k,lambda,t,weights,pivot)
end function sparse_krig_fit_wendland

function sparse_krig_fit_taper(x,y,lambda,a_range,taper_range,model,smoothness,power,korder,t,weights,pivot) result(fit)
real(dp), intent(in) :: x(:,:),y(:),lambda,a_range,taper_range
character(len=*), intent(in), optional :: model,pivot
real(dp), intent(in), optional :: smoothness,power,t(:,:),weights(:)
integer, intent(in), optional :: korder
type(sparse_krig_fit) :: fit
type(csr_matrix) :: k
k=stationary_taper_sparse(x,a_range,taper_range,model,smoothness,power,korder)
fit=sparse_krig_fit_covariance(y,k,lambda,t,weights,pivot)
end function sparse_krig_fit_taper

function sparse_krig_predict(fit,k_new,t_new) result(pred)
type(sparse_krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:)
real(dp), intent(in), optional :: t_new(:,:)
real(dp), allocatable :: pred(:),tt(:,:)
integer :: m
m=size(k_new,1); if(size(k_new,2)/=fit%n) error stop 'sparse_krig_predict: K_new mismatch'
if(present(t_new)) then
   if(size(t_new,1)/=m .or. size(t_new,2)/=fit%p) error stop 'sparse_krig_predict: T_new mismatch'
   tt=t_new
else
   allocate(tt(m,fit%p)); tt=0.0_dp; if(fit%p==1) tt(:,1)=1.0_dp
end if
pred=matmul(tt,fit%beta)+matmul(k_new,fit%c)
end function sparse_krig_predict

function sparse_krig_predict_covariance(fit,k_new,k_newnew,t_new,include_nugget,new_weights) result(v)
type(sparse_krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:),k_newnew(:,:)
real(dp), intent(in), optional :: t_new(:,:),new_weights(:)
logical, intent(in), optional :: include_nugget
real(dp), allocatable :: v(:,:),tt(:,:),qk(:,:),u(:,:)
integer :: m,i
logical :: inc
m=size(k_new,1)
if(size(k_new,2)/=fit%n .or. size(k_newnew,1)/=m .or. size(k_newnew,2)/=m) error stop 'sparse_krig_predict_covariance: mismatch'
if(present(t_new)) then; tt=t_new; else; allocate(tt(m,fit%p)); tt=0.0_dp; if(fit%p==1) tt(:,1)=1.0_dp; end if
qk=spam_solve(fit%factor,transpose(k_new))
u=tt-matmul(k_new,fit%qt)
v=fit%sigma2*(k_newnew-matmul(k_new,qk)+matmul(matmul(u,fit%h_inv),transpose(u)))
v=0.5_dp*(v+transpose(v))
inc=.false.; if(present(include_nugget)) inc=include_nugget
if(inc) then
   if(present(new_weights)) then
      if(size(new_weights)/=m .or. any(new_weights<=0.0_dp)) error stop 'sparse_krig_predict_covariance: invalid weights'
      do i=1,m; v(i,i)=v(i,i)+fit%tau2/new_weights(i); end do
   else
      do i=1,m; v(i,i)=v(i,i)+fit%tau2; end do
   end if
end if
end function sparse_krig_predict_covariance

function sparse_krig_predict_se(fit,k_new,k_newnew,t_new,include_nugget,new_weights) result(se)
type(sparse_krig_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:),k_newnew(:,:)
real(dp), intent(in), optional :: t_new(:,:),new_weights(:)
logical, intent(in), optional :: include_nugget
real(dp), allocatable :: se(:),v(:,:)
integer :: i
v=sparse_krig_predict_covariance(fit,k_new,k_newnew,t_new,include_nugget,new_weights)
allocate(se(size(v,1))); do i=1,size(v,1); se(i)=sqrt(max(0.0_dp,v(i,i))); end do
end function sparse_krig_predict_se

function sparse_krig_multi_fit_covariance(y,k,lambda,t,weights,pivot) result(fit)
real(dp), intent(in) :: y(:,:),lambda
type(csr_matrix), intent(in) :: k
real(dp), intent(in), optional :: t(:,:),weights(:)
character(len=*), intent(in), optional :: pivot
type(sparse_krig_multi_fit) :: fit
type(sparse_krig_fit) :: one
real(dp), allocatable :: r(:),rhs(:)
integer :: j
one=sparse_krig_fit_covariance(y(:,1),k,lambda,t,weights,pivot)
if(one%info/=0) then; fit%info=one%info; return; end if
fit%n=one%n; fit%p=one%p; fit%nresponse=size(y,2); fit%lambda=lambda; fit%info=0
fit%k=k; fit%factor=one%factor; fit%weights=one%weights; fit%t=one%t; fit%qt=one%qt; fit%h_inv=one%h_inv
allocate(fit%beta(fit%p,fit%nresponse),fit%c(fit%n,fit%nresponse),fit%fitted(fit%n,fit%nresponse), &
         fit%residuals(fit%n,fit%nresponse),fit%sigma2(fit%nresponse),fit%tau2(fit%nresponse))
fit%beta(:,1)=one%beta; fit%c(:,1)=one%c; fit%fitted(:,1)=one%fitted; fit%residuals(:,1)=one%residuals
fit%sigma2(1)=one%sigma2; fit%tau2(1)=one%tau2
do j=2,fit%nresponse
   rhs=matmul(transpose(fit%t),spam_solve(fit%factor,y(:,j)))
   fit%beta(:,j)=matmul(fit%h_inv,rhs)
   r=y(:,j)-matmul(fit%t,fit%beta(:,j)); fit%c(:,j)=spam_solve(fit%factor,r)
   fit%fitted(:,j)=matmul(fit%t,fit%beta(:,j))+csr_matvec(k,fit%c(:,j)); fit%residuals(:,j)=y(:,j)-fit%fitted(:,j)
   fit%sigma2(j)=dot_product(r,fit%c(:,j))/max(1.0_dp,real(fit%n-fit%p,dp)); fit%tau2(j)=lambda*fit%sigma2(j)
end do
end function sparse_krig_multi_fit_covariance

function sparse_krig_multi_predict(fit,k_new,t_new) result(pred)
type(sparse_krig_multi_fit), intent(in) :: fit
real(dp), intent(in) :: k_new(:,:)
real(dp), intent(in), optional :: t_new(:,:)
real(dp), allocatable :: pred(:,:),tt(:,:)
integer :: m
m=size(k_new,1); if(size(k_new,2)/=fit%n) error stop 'sparse_krig_multi_predict: mismatch'
if(present(t_new)) then; tt=t_new; else; allocate(tt(m,fit%p)); tt=0.0_dp; if(fit%p==1) tt(:,1)=1.0_dp; end if
pred=matmul(tt,fit%beta)+matmul(k_new,fit%c)
end function sparse_krig_multi_predict

end module fields_sparse_kriging
