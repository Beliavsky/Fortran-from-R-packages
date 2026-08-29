program test_multinom_counts
use nnet, only: dp,multinom_model_t,multinom_fit_counts,multinom_predict_proba
implicit none
type(multinom_model_t)::a,b,c
real(dp)::x(8,2),cnt(8,3),allowed(8,3),off(8,3)
real(dp),allocatable::pa(:,:),pb(:,:),pc(:,:)
integer::i
x(:,1)=1.0_dp
x(:,2)=[-2._dp,-1.5_dp,-1._dp,-0.3_dp,0.2_dp,0.8_dp,1.4_dp,2._dp]
cnt=transpose(reshape([5._dp,2._dp,1._dp,4._dp,3._dp,1._dp,2._dp,5._dp,1._dp,2._dp,3._dp,3._dp, &
 1._dp,4._dp,3._dp,1._dp,3._dp,4._dp,1._dp,2._dp,5._dp,1._dp,1._dp,6._dp],[3,8]))
call multinom_fit_counts(a,x,cnt,maxit=400,reltol=1e-9_dp)
pa=multinom_predict_proba(a,x)
if(maxval(abs(sum(pa,dim=2)-1._dp))>1e-12_dp) error stop 'counts probs'
! Independent SciPy BFGS reference for this baseline-category multinomial likelihood.
if(abs(a%net%value-61.44892130753712_dp)>5e-7_dp) then
   print *,'reference objective',a%net%value
   error stop 'counts reference objective'
end if
if(maxval(abs(a%coefficients-reshape([0.53347279_dp,0.41497049_dp, &
   0.40661046_dp,1.08252913_dp],[2,2])))>6e-5_dp) error stop 'counts reference coefficients' 
off=0._dp
call multinom_fit_counts(b,x,cnt,offset=off,maxit=400,reltol=1e-9_dp)
pb=multinom_predict_proba(b,x,off)
if(maxval(abs(pa-pb))>2e-5_dp) then
print *,maxval(abs(pa-pb))
error stop 'zero offset parity'
end if
allowed=0._dp
do i=1,8
 allowed(i,mod(i-1,3)+1)=1._dp
 allowed(i,mod(i,3)+1)=1._dp
end do
call multinom_fit_counts(c,x,allowed,censored=.true.,maxit=200,reltol=1e-8_dp)
pc=multinom_predict_proba(c,x)
if(any(pc<0._dp).or.maxval(abs(sum(pc,dim=2)-1._dp))>1e-12_dp) error stop 'censored fit'
print *,'test_multinom_counts passed'
end program
