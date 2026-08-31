program test_tps
use fields, only: dp,tps_fit_type,tps_fit,tps_predict,tps_covariance
implicit none
integer,parameter::n=9
real(dp)::x(n,2),y(n),xn(2,2),card(3,2)
real(dp),allocatable::pred(:),c(:,:)
type(tps_fit_type)::fit
integer::i,j,k
k=0
do j=0,2;do i=0,2;k=k+1;x(k,:)=[real(i,dp),real(j,dp)];y(k)=1+2*x(k,1)-0.5_dp*x(k,2);end do;end do
fit=tps_fit(x,y,0.3_dp,m=2)
call check(fit%info==0,'tps solve')
call check(maxval(abs(fit%fitted-y))<1e-9_dp,'plane null space')
xn=reshape([0.25_dp,1.4_dp,0.75_dp,0.2_dp],[2,2])
pred=tps_predict(fit,xn)
call check(maxval(abs(pred-(1+2*xn(:,1)-0.5_dp*xn(:,2))))<1e-8_dp,'tps predict plane')
card=reshape([0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp],[3,2])
c=tps_covariance(xn,xn,card,m=2)
call check(maxval(abs(c-transpose(c)))<1e-10_dp,'Tps.cov symmetry')
print *,'test_tps: PASS'
contains
subroutine check(ok,msg)
logical,intent(in)::ok;character(len=*),intent(in)::msg
if(.not.ok)then;print *,'FAIL: ',msg;error stop 1;end if
end subroutine
end program
