program test_additional
use fields
implicit none
real(dp)::ar,x(8,1),y(8),t(8,2),gtrue(2),xn(2,1),yy(6,2)
real(dp),allocatable::grad(:,:)
integer::i,lev(6)
type(krig_fit)::fit
type(fast_oneway_result)::ow
ar=matern_cor_to_range(1.0_dp,0.5_dp,0.5_dp)
call check(abs(ar-1.0_dp/log(2.0_dp))<1e-6_dp,'matern correlation range')
do i=1,8
 x(i,1)=real(i-1,dp)/7.0_dp;y(i)=2.0_dp+3.0_dp*x(i,1)
end do
t(:,1)=1.0_dp;t(:,2)=x(:,1)
fit=krig_fit_stationary(x,y,1.0_dp,model='gaussian',a_range=0.4_dp,t=t)
xn(:,1)=[0.25_dp,0.75_dp]
grad=krig_predict_gradient_stationary(fit,x,xn,model='gaussian',a_range=0.4_dp,m=2)
gtrue=3.0_dp
call check(maxval(abs(grad(:,1)-gtrue))<2e-6_dp,'krig derivative')
lev=[2,2,1,1,2,1]
yy(:,1)=[1d0,3d0,2d0,4d0,5d0,6d0];yy(:,2)=2d0*yy(:,1)
ow=fast_oneway(lev,yy)
call check(ow%ngroups==2 .and. all(ow%tags==[2,1]),'fast oneway tags')
call check(abs(ow%means(1,1)-3d0)<1e-14_dp,'fast oneway mean')
print *,'test_additional: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(*),intent(in)::msg
if(.not.ok)then;print *,'FAIL ',trim(msg);error stop;end if
end subroutine
end program
