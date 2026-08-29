program test_multinom
use nnet, only: dp,multinom_model_t,multinom_fit_labels,multinom_predict_proba,multinom_predict_class
implicit none
type(multinom_model_t)::fit
real(dp)::x(12,3)
integer::lab(12)
real(dp),allocatable::p(:,:),id(:,:)
integer,allocatable::cl(:)
x(:,1)=1.0_dp
x(:,2)=[-2.0_dp,-1.5_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp,-1.2_dp,0.2_dp,1.2_dp]
x(:,3)=[0.0_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,-1.5_dp,1.2_dp,-1.2_dp]
lab=[1,1,2,2,3,2,2,3,3,1,3,2]
call multinom_fit_labels(fit,x,lab,maxit=500,reltol=1e-9_dp,hessian=.true.)
p=multinom_predict_proba(fit,x)
cl=multinom_predict_class(fit,x)
if(maxval(abs(sum(p,dim=2)-1.0_dp))>1e-12_dp) error stop 'multinom probability sum'
if(count(cl==lab)<9) then
print *,cl
error stop 'multinom poor fit'
end if
if(fit%edf/=6) error stop 'multinom edf'
if(maxval(abs(fit%information-transpose(fit%information)))>1e-10_dp) error stop 'information symmetry'
allocate(id(size(fit%information,1),size(fit%information,2)))
id=matmul(fit%information,fit%covariance)
if(maxval(abs(id-diag_identity(size(id,1))))>2e-5_dp) then
 print *,'inverse err',maxval(abs(id-diag_identity(size(id,1))))
 error stop 'multinom covariance'
end if
print *,'test_multinom passed; objective=',fit%net%value,' correct=',count(cl==lab)
contains
pure function diag_identity(n) result(a)
integer,intent(in)::n
real(dp)::a(n,n)
integer::j
a=0.0_dp
do j=1,n
a(j,j)=1.0_dp
end do
end function
end program
