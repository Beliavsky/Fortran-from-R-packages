! GPL-2.0-or-later. fastTps computational wrapper using Wendland + spam.
module fields_fast_tps
use fields_kinds, only: dp
use fields_polynomial, only: polynomial_basis
use fields_sparse_kriging, only: sparse_krig_fit, sparse_krig_fit_wendland, sparse_krig_predict
use fields_covariance, only: wendland_covariance
implicit none
private
public :: fast_tps_fit, fast_tps_predict

contains

function fast_tps_fit(x,y,a_range,lambda,m,korder,weights,pivot) result(fit)
real(dp),intent(in)::x(:,:),y(:),a_range,lambda
integer,intent(in),optional::m,korder
real(dp),intent(in),optional::weights(:)
character(len=*),intent(in),optional::pivot
type(sparse_krig_fit)::fit
real(dp),allocatable::t(:,:)
integer::mm,kk,info,d
if(a_range<=0.0_dp)error stop 'fast_tps_fit: a_range must be positive'
d=size(x,2);mm=max(2,int(ceiling(0.5_dp*real(d,dp)+0.1_dp)));if(present(m))mm=m
kk=max(2,2*mm-d);if(present(korder))kk=korder
call polynomial_basis(x,mm,t,info=info);if(info/=0)error stop 'fast_tps_fit: polynomial basis failed'
fit=sparse_krig_fit_wendland(x,y,lambda,a_range,kk,t,weights,pivot)
end function fast_tps_fit

function fast_tps_predict(fit,xtrain,xnew,a_range,m,korder) result(pred)
type(sparse_krig_fit),intent(in)::fit
real(dp),intent(in)::xtrain(:,:),xnew(:,:),a_range
integer,intent(in),optional::m,korder
real(dp),allocatable::pred(:),tnew(:,:),knew(:,:)
integer::mm,kk,d,info
d=size(xtrain,2);mm=max(2,int(ceiling(0.5_dp*real(d,dp)+0.1_dp)));if(present(m))mm=m
kk=max(2,2*mm-d);if(present(korder))kk=korder
call polynomial_basis(xnew,mm,tnew,info=info);if(info/=0)error stop 'fast_tps_predict: polynomial basis failed'
knew=wendland_covariance(xnew,xtrain,a_range,kk)
pred=sparse_krig_predict(fit,knew,tnew)
end function fast_tps_predict

end module fields_fast_tps
