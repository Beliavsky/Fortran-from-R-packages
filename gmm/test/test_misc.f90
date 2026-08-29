program test_misc
use gmm, only: dp,char_stable,gel_result_t,gel_lambda,GEL_EL,GEL_CUE,kernel_weight,moment_covariance
implicit none
real(dp)::theta(4),tau(2),gt(6,2),s(2,2)
complex(dp)::cf(2)
type(gel_result_t)::lr
theta=[1.5_dp,0.3_dp,2.0_dp,0.5_dp]
tau=[0.2_dp,-0.7_dp]
cf=char_stable(theta,tau,1)
call ac(real(cf(1)),0.77625611_dp,2e-8_dp)
call ac(aimag(cf(1)),0.01871554_dp,2e-8_dp)
call ac(real(cf(2)),0.18874872_dp,2e-8_dp)
call ac(aimag(cf(2)),0.02793815_dp,2e-8_dp)
gt=reshape([ &
 -1.0_dp,-0.5_dp, 0.2_dp,0.8_dp,1.0_dp,1.4_dp, &
  0.5_dp,-0.2_dp,-0.4_dp,0.3_dp,0.8_dp,-0.6_dp],[6,2])
call gel_lambda(gt,GEL_EL,lr,tol_lambda=1e-11_dp,maxit=200)
if(lr%convergence/=0)error stop 'EL lambda failed'
if(maxval(abs(sum(gt*spread(lr%prob,2,size(gt,2)),dim=1)))>1e-7_dp)error stop 'EL moments not balanced'
call gel_lambda(gt,GEL_CUE,lr)
if(lr%convergence/=0)error stop 'CUE lambda failed'
call ac(kernel_weight(0.5_dp,'Bartlett'),0.5_dp,1e-14_dp)
s=moment_covariance(gt,.true.)
if(s(1,1)<=0.0_dp .or. s(2,2)<=0.0_dp)error stop 'bad covariance'
print '(a)','test_misc: ok'
contains
subroutine ac(a,b,t)
real(dp),intent(in)::a,b,t
if(abs(a-b)>t)then
print *,a,b
error stop 1
end if
end subroutine
end program
