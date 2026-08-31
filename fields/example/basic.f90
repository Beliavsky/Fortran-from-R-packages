program fields_basic
use fields, only: dp,krig_fit,krig_fit_stationary_reml,stationary_covariance,krig_predict
implicit none
integer,parameter::n=20
real(dp)::x(n,1),y(n),xnew(5,1)
real(dp),allocatable::knew(:,:),pred(:)
type(krig_fit)::fit
integer::i
do i=1,n;x(i,1)=real(i-1,dp)/real(n-1,dp);y(i)=sin(6*x(i,1))+0.1_dp*cos(17*x(i,1));end do
fit=krig_fit_stationary_reml(x,y,model='matern',a_range=0.25_dp,smoothness=1.5_dp)
do i=1,5;xnew(i,1)=real(i-1,dp)/4.0_dp;end do
knew=stationary_covariance(xnew,x,model='matern',a_range=0.25_dp,smoothness=1.5_dp)
pred=krig_predict(fit,knew)
write(*,'(a,es12.4)') 'lambda = ',fit%lambda
write(*,'(a,5f10.5)') 'predictions = ',pred
end program
