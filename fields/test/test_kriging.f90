program test_kriging
use fields, only: dp,stationary_covariance,krig_fit,krig_fit_covariance,krig_predict,krig_predict_se
implicit none
integer,parameter::n=8
real(dp)::x(n,1),y(n),t(n,2),tn(3,2),xn(3,1)
real(dp),allocatable::k(:,:),kn(:,:),knn(:,:),pred(:),se(:)
type(krig_fit)::fit
integer::i
do i=1,n;x(i,1)=real(i-1,dp)/4.0_dp;y(i)=1.0_dp+0.5_dp*x(i,1)+sin(2*x(i,1));end do
t(:,1)=1.0_dp;t(:,2)=x(:,1)
k=stationary_covariance(x,x,model='exponential',a_range=0.7_dp)
fit=krig_fit_covariance(y,k,0.08_dp,t)
call check(fit%info==0,'fit info')
call check(fit%sigma2>0 .and. fit%tau2>0,'variance estimates')
call check(fit%trace_a>2 .and. fit%trace_a<real(n,dp),'trace range')
xn(:,1)=[0.15_dp,0.75_dp,1.4_dp];tn(:,1)=1;tn(:,2)=xn(:,1)
kn=stationary_covariance(xn,x,model='exponential',a_range=0.7_dp)
knn=stationary_covariance(xn,xn,model='exponential',a_range=0.7_dp)
pred=krig_predict(fit,kn,tn);se=krig_predict_se(fit,kn,knn,tn)
call check(all(pred==pred),'prediction finite')
call check(all(se>=0.0_dp) .and. maxval(se)<10.0_dp,'prediction se')
call check(abs(dot_product(fit%t(:,1),fit%c))<1e-9_dp,'universal kriging constraint intercept')
call check(abs(dot_product(fit%t(:,2),fit%c))<1e-9_dp,'universal kriging constraint slope')
print *,'test_kriging: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok)then;print *,'FAIL: ',msg;error stop 1;end if
end subroutine
end program
