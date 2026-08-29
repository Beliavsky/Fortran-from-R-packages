program test_distribution
use tweedie
use r_compat, only: dpois, dgamma
implicit none
real(dp) :: x, d, p, q, ref
integer :: failures
failures = 0

! Public p=1 behavior: a scaled Poisson lattice, matching tweedie::dtweedie.
x = 2.0_dp
d = dtweedie(x, 4.0_dp, 2.0_dp, 1.0_dp)
ref = dpois(x / 2.0_dp, lambda=2.0_dp)
call check_close('p=1 density', d, ref, 2.0e-14_dp, failures)
p = ptweedie(x, 4.0_dp, 2.0_dp, 1.0_dp)
q = qtweedie(p, 4.0_dp, 2.0_dp, 1.0_dp)
if (q > x + 1.0e-12_dp) call fail('p=1 quantile relation', failures)

! Gamma special case.
x = 1.3_dp
d = dtweedie(x, 1.4_dp, 0.7_dp, 2.0_dp)
ref = dgamma(x, shape=1.0_dp/0.7_dp, rate=1.0_dp/(0.7_dp*1.4_dp))
call check_close('p=2 gamma density', d, ref, 2.0e-14_dp, failures)

! Inverse-Gaussian special case and p/q consistency.
x = 1.0_dp
d = dtweedie(x, 1.4_dp, 0.74_dp, 3.0_dp)
call check_close('p=3 IG density reference', d, 0.4388738851125460_dp, 5.0e-14_dp, failures)
p = ptweedie(x, 1.4_dp, 0.74_dp, 3.0_dp)
call check_close('p=3 IG cdf reference', p, 0.5294018089021599_dp, 5.0e-13_dp, failures)
q = qtweedie(p, 1.4_dp, 0.74_dp, 3.0_dp)
call check_close('p=3 p/q inverse', q, x, 2.0e-10_dp, failures)

! Compound-Poisson mass at zero.
ref = exp(-0.5_dp)
call check_close('p=1.5 zero mass', dtweedie(0.0_dp,1.0_dp,4.0_dp,1.5_dp), ref, &
   2.0e-14_dp, failures)
call check_close('p=1.5 cdf at zero', ptweedie(0.0_dp,1.0_dp,4.0_dp,1.5_dp), ref, &
   2.0e-14_dp, failures)

! Deviance and compound-Poisson conversion.
call check_close('Poisson deviance', tweedie_dev(2.0_dp,1.0_dp,1.0_dp), &
   2.0_dp*(2.0_dp*log(2.0_dp)-1.0_dp), 2.0e-14_dp, failures)
call check_close('lambda', tweedie_lambda(1.0_dp,4.0_dp,1.5_dp), 0.5_dp, 2.0e-14_dp, failures)

if (failures /= 0) error stop 'test_distribution failed'
print '(a)', 'test_distribution: PASS'
contains
subroutine check_close(name, got, expected, tol, failures)
character(*),intent(in)::name
real(dp),intent(in)::got,expected,tol
integer,intent(inout)::failures
if(abs(got-expected)>tol*max(1.0_dp,abs(expected)))then
 print '(a,2es24.14)', trim(name)//' got/expected: ',got,expected
 failures=failures+1
end if
end subroutine
subroutine fail(name, failures)
character(*),intent(in)::name
integer,intent(inout)::failures
print '(a)',trim(name)//': FAIL'
failures=failures+1
end subroutine
end program test_distribution
