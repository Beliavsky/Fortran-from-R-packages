program test_regression
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use tweedie
implicit none
real(dp)::v
integer::failures
failures=0
call check_density(0.5557088_dp,1.0_dp,0.5513086_dp,1.5_dp)
call check_cdf(0.01_dp,5.0_dp,0.01_dp,1.01_dp)
call check_cdf(0.095095_dp,1.0_dp,1.0_dp,4.5_dp)
call check_cdf(0.001_dp,0.01_dp,0.01_dp,6.0_dp)
call check_cdf(0.001_dp,5.0_dp,0.01_dp,6.0_dp)
call check_cdf(0.075075_dp,1.0_dp,1.0_dp,4.5_dp)
call check_cdf(0.001_dp,5.0_dp,10.0_dp,1.01_dp)
call check_cdf(0.001_dp,1.0_dp,0.01_dp,1.01_dp)
call check_cdf(0.3_dp,2.0_dp,1.0_dp,2.5_dp)
call check_cdf(0.3_dp,2.0_dp,1.0_dp,1.5_dp)
call check_cdf(0.3_dp,2.0_dp,1.0_dp,3.5_dp)
call check_cdf(7.709933e-308_dp,10.17691_dp,4.55_dp,1.98_dp)
call check_cdf(37470.91_dp,40597.51_dp,2161.194_dp,1.676_dp)
call check_cdf(50.0_dp,1.4_dp,0.74_dp,3.0_dp)
if(failures/=0)error stop 'test_regression failed'
print '(a)','test_regression: PASS'
contains
subroutine check_density(y,mu,phi,power)
real(dp),intent(in)::y,mu,phi,power
v=dtweedie_inversion(y,mu,phi,power)
if(.not.ieee_is_finite(v).or.v<0.0_dp)then
 print '(a,4es15.6)','bad density: ',y,mu,phi,power
 failures=failures+1
end if
end subroutine
subroutine check_cdf(q,mu,phi,power)
real(dp),intent(in)::q,mu,phi,power
v=ptweedie_inversion(q,mu,phi,power)
if(.not.ieee_is_finite(v).or.v<0.0_dp.or.v>1.0_dp)then
 print '(a,5es15.6)','bad cdf: ',q,mu,phi,power,v
 failures=failures+1
end if
end subroutine
end program test_regression
