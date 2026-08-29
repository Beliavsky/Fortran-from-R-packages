program test_linear
use gmm, only: dp,linear_gmm_result_t,linear_gmm_fit,tsls_fit,LINEAR_TWO_STEP,LINEAR_ITERATIVE,LINEAR_CUE
implicit none
real(dp)::y(8),x(8,2),z(8,3)
type(linear_gmm_result_t)::r,r2,r3
integer::i
y=[1.2_dp,1.9_dp,2.7_dp,3.8_dp,4.1_dp,5.2_dp,5.9_dp,7.1_dp]
x(:,1)=1.0_dp
x(:,2)=[(real(i-1,dp),i=1,8)]
z(:,1)=1.0_dp
z(:,2)=x(:,2)
z(:,3)=[1.0_dp,-1.0_dp,1.0_dp,-1.0_dp,1.0_dp,-1.0_dp,1.0_dp,-1.0_dp]
call linear_gmm_fit(y,x,z,r,method=LINEAR_TWO_STEP,covariance=1)
call assert_close(r%coefficients(1),1.18830752_dp,2e-6_dp,'linear two-step intercept')
call assert_close(r%coefficients(2),0.81164113_dp,2e-6_dp,'linear two-step slope')
call tsls_fit(y,x,z,r2)
call assert_close(r2%coefficients(1),1.10833333_dp,2e-6_dp,'tsls intercept')
call assert_close(r2%coefficients(2),0.82261905_dp,2e-6_dp,'tsls slope')
if(r%j_stat<0.0_dp .or. r%j_pvalue<0.0_dp .or. r%j_pvalue>1.0_dp) error stop 'bad J test'

call linear_gmm_fit(y,x,z,r3,method=LINEAR_ITERATIVE,covariance=1,tol=1e-11_dp)
call assert_close(r3%coefficients(1),1.12195599_dp,3e-5_dp,'linear iterative intercept')
call assert_close(r3%coefficients(2),0.83462069_dp,3e-5_dp,'linear iterative slope')
call linear_gmm_fit(y,x,z,r3,method=LINEAR_CUE,covariance=1)
call assert_close(r3%coefficients(1),1.12216693_dp,4e-4_dp,'linear CUE intercept')
call assert_close(r3%coefficients(2),0.84219069_dp,4e-4_dp,'linear CUE slope')

print '(a)','test_linear: ok'
contains
subroutine assert_close(a,b,tol,msg)
real(dp),intent(in)::a,b,tol
character(len=*),intent(in)::msg
if(abs(a-b)>tol)then
print *,msg,a,b
error stop 1
end if
end subroutine
end program
