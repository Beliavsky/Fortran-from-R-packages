program test_profile
use tweedie
implicit none
real(dp)::x(6,1),y(6),powers(3),expected_mu
real(dp),allocatable::beta(:),mu(:)
logical::ok
integer::it,failures
type(tweedie_profile_result)::out
failures=0
x(:,1)=1.0_dp
y=[0.2_dp,0.5_dp,0.9_dp,1.1_dp,1.5_dp,2.0_dp]
expected_mu=sum(y)/real(size(y),dp)
call tweedie_glm_fit(x,y,1.5_dp,0.0_dp,beta,mu,ok,it)
if(.not.ok)then
print *,'GLM did not converge'
failures=failures+1
end if
call check_close('intercept-only fitted mean',mu(1),expected_mu,2.0e-9_dp,failures)
if(maxval(abs(mu-expected_mu))>2.0e-9_dp)then
 print *,'intercept-only fitted means differ'
 failures=failures+1
end if
powers=[1.3_dp,1.5_dp,1.7_dp]
call tweedie_profile_grid(x,y,powers,0.0_dp,out,method=tweedie_method_series)
if(any(.not.out%converged))then
print *,'profile GLM convergence failure'
failures=failures+1
end if
if(any(out%phi<=0.0_dp))then
print *,'profile nonpositive phi'
failures=failures+1
end if
if(any(out%loglik<=-huge(1.0_dp)/2.0_dp))then
print *,'profile invalid likelihood'
failures=failures+1
end if
if(failures/=0)error stop 'test_profile failed'
print '(a)','test_profile: PASS'
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
end program test_profile
