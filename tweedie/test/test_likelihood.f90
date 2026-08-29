program test_likelihood
use tweedie
implicit none
real(dp)::y(5),mu(5),phi,p,ll1,ll2,fd,analytic,phat,aic
integer::failures
failures=0
y=[0.2_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp]
mu=1.0_dp
phi=0.8_dp
p=1.5_dp
ll1=tweedie_loglik(y,mu,phi,p,tweedie_method_series)
ll2=tweedie_loglik(y,mu,phi+1.0e-6_dp,p,tweedie_method_series)
fd=(ll2-ll1)/1.0e-6_dp
analytic=-0.5_dp*dtweedie_dldphi(phi,mu,p,y)
call check_close('phi score',analytic,fd,2.0e-5_dp,failures)
phat=tweedie_phi_mle(y,mu,p,tweedie_method_series)
if(phat<=0.0_dp.or.phat>=2.0_dp)then
 print *, 'unexpected phi mle',phat
 failures=failures+1
end if
if(tweedie_loglik(y,mu,phat,p,tweedie_method_series)<ll1)then
 print *, 'phi mle did not improve likelihood'
 failures=failures+1
end if
aic=tweedie_aic(ll1,3)
call check_close('AIC',aic,-2.0_dp*ll1+6.0_dp,2.0e-14_dp,failures)
if(failures/=0)error stop 'test_likelihood failed'
print '(a)','test_likelihood: PASS'
contains
subroutine check_close(name,got,expected,tol,failures)
character(*),intent(in)::name
real(dp),intent(in)::got,expected,tol
integer,intent(inout)::failures
if(abs(got-expected)>tol*max(1.0_dp,abs(expected)))then
 print '(a,2es24.14)',trim(name)//' got/expected: ',got,expected
 failures=failures+1
end if
end subroutine
end program test_likelihood
