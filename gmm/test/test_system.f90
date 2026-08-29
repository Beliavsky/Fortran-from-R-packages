program test_system
use gmm, only: dp,system_gmm_result_t,system_gmm_fit,SYS_CONDHOM
implicit none
real(dp)::y(8,2),x(8,2,2),z(8,3,2),t(8)
type(system_gmm_result_t)::r
integer::i
t=[(real(i-1,dp),i=1,8)]
x(:,1,1)=1.0_dp
x(:,2,1)=t
x(:,1,2)=1.0_dp
x(:,2,2)=[0.0_dp,1.0_dp,1.0_dp,2.0_dp,3.0_dp,5.0_dp,5.0_dp,7.0_dp]
z(:,1,1)=1.0_dp
z(:,2,1)=x(:,2,1)
z(:,3,1)=[1.0_dp,-1.0_dp,1.0_dp,-1.0_dp,1.0_dp,-1.0_dp,1.0_dp,-1.0_dp]
z(:,1,2)=1.0_dp
z(:,2,2)=x(:,2,2)
z(:,3,2)=[0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp]
y(:,1)=1.0_dp+0.8_dp*t+[0.2_dp,-0.1_dp,0.1_dp,0.0_dp,0.3_dp,-0.2_dp,0.1_dp,-0.1_dp]
y(:,2)=-0.5_dp+1.2_dp*x(:,2,2)+[-0.1_dp,0.2_dp,-0.2_dp,0.1_dp,0.0_dp,0.1_dp,-0.1_dp,0.0_dp]
call system_gmm_fit(y,x,z,r,covariance=SYS_CONDHOM)
call ac(r%coefficients(1),1.11659311_dp,4e-5_dp)
call ac(r%coefficients(2),0.77740197_dp,4e-5_dp)
call ac(r%coefficients(3),-0.50915008_dp,4e-5_dp)
call ac(r%coefficients(4),1.20305003_dp,4e-5_dp)
print '(a)','test_system: ok'
contains
subroutine ac(a,b,t)
real(dp),intent(in)::a,b,t
if(abs(a-b)>t)then
print *,a,b
error stop 1
end if
end subroutine
end program
