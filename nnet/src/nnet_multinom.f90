! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Matrix-oriented translation of the computational parts of multinom.R and vcovmultinom.R.
module nnet_multinom
use r_compat, only: dp, qr_fit_t, qr, qr_rank
use nnet_types, only: nnet_model_t, multinom_model_t
use nnet_fit_mod, only: nnet_fit, nnet_predict
use nnet_utils, only: class_ind
use r_linalg, only: shared_full_svd => full_svd
implicit none
private
public :: multinom_fit_labels, multinom_fit_counts, multinom_predict_proba, multinom_predict_class
public :: multinom_information, multinom_covariance, multinom_loglik
contains

subroutine multinom_fit_labels(fit,x,labels,case_weights,offset,decay,maxit,reltol,hessian)
type(multinom_model_t), intent(out) :: fit
real(dp), intent(in) :: x(:,:)
integer, intent(in) :: labels(:)
real(dp), intent(in), optional :: case_weights(:), offset(:,:), decay, reltol
integer, intent(in), optional :: maxit
logical, intent(in), optional :: hessian
real(dp), allocatable :: y(:,:), cw(:)
integer :: k
if(size(labels)/=size(x,1)) error stop "multinom_fit_labels: row counts differ"
k=maxval(labels)
if(k<2 .or. any(labels<1)) error stop "multinom_fit_labels: need labels 1..K, K>=2"
allocate(cw(size(x,1)))
cw=1.0_dp
if(present(case_weights)) cw=case_weights
if(k==2) then
   allocate(y(size(x,1),1))
   y(:,1)=real(labels-1,dp)
   call multinom_binary_core(fit,x,y,cw,offset,decay,maxit,reltol,hessian)
else
   y=class_ind(labels,k)
   call multinom_softmax_core(fit,x,y,cw,.false.,offset,decay,maxit,reltol,hessian)
end if
end subroutine multinom_fit_labels

subroutine multinom_fit_counts(fit,x,counts,case_weights,censored,offset,decay,maxit,reltol,hessian)
type(multinom_model_t), intent(out) :: fit
real(dp), intent(in) :: x(:,:), counts(:,:)
real(dp), intent(in), optional :: case_weights(:), offset(:,:), decay, reltol
logical, intent(in), optional :: censored, hessian
integer, intent(in), optional :: maxit
real(dp), allocatable :: y(:,:), cw(:), totals(:)
logical :: cens
integer :: i
if(size(counts,1)/=size(x,1) .or. size(counts,2)<2) error stop "multinom_fit_counts: invalid counts shape"
cens=.false.
if(present(censored)) cens=censored
allocate(cw(size(x,1)))
cw=1.0_dp
if(present(case_weights)) cw=case_weights
y=counts
if(.not.cens) then
   allocate(totals(size(x,1)))
   totals=sum(counts,dim=2)
   if(any(totals<=0.0_dp)) error stop "multinom_fit_counts: each row must have positive total"
   do i=1,size(x,1)
   y(i,:)=counts(i,:)/totals(i)
   end do
   cw=cw*totals
end if
call multinom_softmax_core(fit,x,y,cw,cens,offset,decay,maxit,reltol,hessian)
end subroutine multinom_fit_counts

subroutine multinom_binary_core(fit,x,y,cw,offset,decay,maxit,reltol,hessian)
type(multinom_model_t), intent(out) :: fit
real(dp), intent(in) :: x(:,:), y(:,:), cw(:)
real(dp), intent(in), optional :: offset(:,:), decay, reltol
integer, intent(in), optional :: maxit
logical, intent(in), optional :: hessian
real(dp), allocatable :: xx(:,:), init(:), decv(:)
logical, allocatable :: mask(:)
integer :: p,nw
p=size(x,2)
if(present(offset)) then
   if(size(offset,1)/=size(x,1) .or. size(offset,2)/=1) error stop "binary multinom offset must be n x 1"
   allocate(xx(size(x,1),p+1))
   xx(:,1:p)=x
   xx(:,p+1)=offset(:,1)
else
   xx=x
end if
nw=size(xx,2)+1
allocate(init(nw),mask(nw))
init=0.0_dp
mask=.true.
mask(1)=.false.
if(present(offset)) then
init(nw)=1.0_dp
mask(nw)=.false.
end if
allocate(decv(1))
decv=0.0_dp
if(present(decay)) decv=decay
call nnet_fit(fit%net,xx,y,0,case_weights=cw,initial_wts=init,mask=mask,entropy=.true.,skip=.true., &
   rang=0.0_dp,decay=decv,maxit=maxit_value(maxit),reltol=reltol_value(reltol),hessian=hess_value(hessian))
fit%n_classes=2
call finish_multinom(fit,x,cw)
end subroutine multinom_binary_core

subroutine multinom_softmax_core(fit,x,y,cw,censored,offset,decay,maxit,reltol,hessian)
type(multinom_model_t), intent(out) :: fit
real(dp), intent(in) :: x(:,:), y(:,:), cw(:)
logical, intent(in) :: censored
real(dp), intent(in), optional :: offset(:,:), decay, reltol
integer, intent(in), optional :: maxit
logical, intent(in), optional :: hessian
real(dp), allocatable :: xx(:,:), init(:), decv(:)
logical, allocatable :: mask(:)
integer :: p,k,nin,block,c,pos
p=size(x,2)
k=size(y,2)
if(present(offset)) then
   if(size(offset,1)/=size(x,1) .or. size(offset,2)/=k) error stop "softmax offset must be n x K"
   allocate(xx(size(x,1),p+k))
   xx(:,1:p)=x
   xx(:,p+1:p+k)=offset
else
   xx=x
end if
nin=size(xx,2)
block=nin+1
allocate(init(k*block),mask(k*block))
init=0.0_dp
mask=.false.
! Baseline class 1 is fixed at zero. Other classes have free original-X coefficients,
! fixed bias zero, and (when supplied) fixed own-class offset coefficient one.
do c=2,k
   pos=(c-1)*block
   mask(pos+2:pos+1+p)=.true.
end do
if(present(offset)) then
   do c=1,k
      pos=(c-1)*block
      init(pos+1+p+c)=1.0_dp
   end do
end if
allocate(decv(1))
decv=0.0_dp
if(present(decay)) decv=decay
call nnet_fit(fit%net,xx,y,0,case_weights=cw,initial_wts=init,mask=mask,softmax=.true.,censored=censored,skip=.true., &
   rang=0.0_dp,decay=decv,maxit=maxit_value(maxit),reltol=reltol_value(reltol),hessian=hess_value(hessian))
fit%n_classes=k
call finish_multinom(fit,x,cw)
end subroutine multinom_softmax_core

subroutine finish_multinom(fit,x,row_weights)
type(multinom_model_t), intent(inout) :: fit
real(dp), intent(in) :: x(:,:), row_weights(:)
type(qr_fit_t) :: qrf
real(dp), allocatable :: probs(:,:)
integer :: p,k,c,block,pos
p=size(x,2)
k=fit%n_classes
qrf=qr(x)
fit%rank=qr_rank(qrf)
fit%edf=(k-1)*fit%rank
fit%deviance=2.0_dp*fit%net%value
fit%aic=fit%deviance+2.0_dp*real(fit%edf,dp)
allocate(fit%coefficients(k-1,p))
fit%coefficients=0.0_dp
if(k==2 .and. fit%net%n_outputs==1) then
   fit%coefficients(1,:)=fit%net%wts(2:1+p)
   allocate(probs(size(x,1),2))
   probs(:,2)=fit%net%fitted(:,1)
   probs(:,1)=1.0_dp-probs(:,2)
else
   block=fit%net%n_inputs+1
   do c=2,k
      pos=(c-1)*block
      fit%coefficients(c-1,:)=fit%net%wts(pos+2:pos+1+p)
   end do
   probs=fit%net%fitted
end if
fit%information=multinom_information(x,probs,row_weights)
fit%covariance=multinom_covariance(fit%information)
end subroutine finish_multinom

pure function multinom_information(x,probs,row_weights) result(info)
real(dp), intent(in) :: x(:,:), probs(:,:), row_weights(:)
real(dp), allocatable :: info(:,:)
integer :: p,k,a,b,i,ia,ib
real(dp) :: fac
p=size(x,2)
k=size(probs,2)
allocate(info((k-1)*p,(k-1)*p))
info=0.0_dp
do i=1,size(x,1)
   do a=2,k
      ia=(a-2)*p
      do b=2,k
         ib=(b-2)*p
         fac=-probs(i,a)*probs(i,b)
         if(a==b) fac=fac+probs(i,a)
         info(ia+1:ia+p,ib+1:ib+p)=info(ia+1:ia+p,ib+1:ib+p)+ &
            row_weights(i)*fac*outer_product(x(i,:),x(i,:))
      end do
   end do
end do
end function multinom_information

function multinom_covariance(info_mat) result(vcov)
real(dp), intent(in) :: info_mat(:,:)
real(dp), allocatable :: vcov(:,:)
real(dp), allocatable :: sval(:), u(:,:), vt(:,:), v(:,:), invs(:)
real(dp) :: tol
integer :: n,lapack_info,i
n=size(info_mat,1)
if(size(info_mat,2)/=n) error stop "multinom_covariance: information matrix must be square"
allocate(vcov(n,n))
vcov=0.0_dp
if(n==0) return
call shared_full_svd(info_mat,u,sval,vt,lapack_info)
if(lapack_info/=0) error stop "multinom_covariance: SVD failed"
tol=sqrt(epsilon(1.0_dp))*max(maxval(sval),0.0_dp)
allocate(invs(n))
invs=0.0_dp
do i=1,n
   if(sval(i)>tol) invs(i)=1.0_dp/sval(i)
end do
v=transpose(vt)
vcov=matmul(v*spread(invs,1,n),transpose(u))
end function multinom_covariance

pure function multinom_predict_proba(fit,x,offset) result(probs)
type(multinom_model_t), intent(in) :: fit
real(dp), intent(in) :: x(:,:)
real(dp), intent(in), optional :: offset(:,:)
real(dp), allocatable :: probs(:,:),xx(:,:),raw(:,:)
integer :: p,k
p=size(x,2)
k=fit%n_classes
if(present(offset)) then
   if(k==2 .and. fit%net%n_inputs==p+1) then
      allocate(xx(size(x,1),p+1))
      xx(:,1:p)=x
      xx(:,p+1)=offset(:,1)
   else if(fit%net%n_inputs==p+k) then
      allocate(xx(size(x,1),p+k))
      xx(:,1:p)=x
      xx(:,p+1:p+k)=offset
   else
      error stop "multinom_predict_proba: incompatible offset"
   end if
else
   if(size(x,2)/=fit%net%n_inputs) error stop "multinom_predict_proba: model requires offsets"
   xx=x
end if
raw=nnet_predict(fit%net,xx)
if(k==2 .and. size(raw,2)==1) then
   allocate(probs(size(x,1),2))
   probs(:,2)=raw(:,1)
   probs(:,1)=1.0_dp-raw(:,1)
else
   probs=raw
end if
end function multinom_predict_proba

function multinom_predict_class(fit,x,offset) result(cls)
type(multinom_model_t), intent(in) :: fit
real(dp), intent(in) :: x(:,:)
real(dp), intent(in), optional :: offset(:,:)
integer, allocatable :: cls(:)
real(dp), allocatable :: p(:,:)
integer :: i
p=multinom_predict_proba(fit,x,offset)
allocate(cls(size(x,1)))
do i=1,size(x,1)
cls(i)=maxloc(p(i,:),dim=1)
end do
end function multinom_predict_class

pure real(dp) function multinom_loglik(fit) result(ll)
type(multinom_model_t), intent(in) :: fit
ll=-0.5_dp*fit%deviance
end function multinom_loglik

pure function outer_product(a,b) result(c)
real(dp), intent(in) :: a(:),b(:)
real(dp) :: c(size(a),size(b))
c=spread(a,2,size(b))*spread(b,1,size(a))
end function outer_product

pure integer function maxit_value(x) result(v)
integer,intent(in),optional::x
v=100
if(present(x))v=x
end function
pure real(dp) function reltol_value(x) result(v)
real(dp),intent(in),optional::x
v=1.0e-8_dp
if(present(x))v=x
end function
pure logical function hess_value(x) result(v)
logical,intent(in),optional::x
v=.false.
if(present(x))v=x
end function

end module nnet_multinom
