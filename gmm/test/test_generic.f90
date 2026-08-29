program test_generic
use gmm, only: dp,gmm_fit_result_t,gmm_fit,GMM_TWO_STEP,COV_MDS
implicit none
real(dp)::data(8,2),theta0(1)
type(gmm_fit_result_t)::r
integer::i
data(:,1)=[(real(i,dp),i=1,8)]
data(:,2)=[0.5_dp,1.2_dp,0.7_dp,1.5_dp,0.9_dp,1.8_dp,1.1_dp,2.0_dp]
theta0=4.0_dp
call gmm_fit(mom,data,theta0,r,method=GMM_TWO_STEP,covariance=COV_MDS,gradient=grad)
call assert_close(r%theta(1),4.67039613_dp,2e-5_dp,'generic GMM')
if(r%df/=1)error stop 'wrong df'
print '(a)','test_generic: ok'
contains
pure function mom(theta,d) result(gt)
real(dp),intent(in)::theta(:),d(:,:)
real(dp),allocatable::gt(:,:)
allocate(gt(size(d,1),2))
gt(:,1)=d(:,1)-theta(1)
gt(:,2)=(d(:,1)-theta(1))*d(:,2)
end function
pure function grad(theta,d) result(g)
real(dp),intent(in)::theta(:),d(:,:)
real(dp),allocatable::g(:,:)
allocate(g(2,1))
g(1,1)=-1.0_dp
g(2,1)=-sum(d(:,2))/real(size(d,1),dp)
end function
subroutine assert_close(a,b,tol,msg)
real(dp),intent(in)::a,b,tol
character(len=*),intent(in)::msg
if(abs(a-b)>tol)then
print *,msg,a,b
error stop 1
end if
end subroutine
end program
