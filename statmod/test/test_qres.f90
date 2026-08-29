program test_qres
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use statmod_qres
use r_compat, only: dp, set_seed_int
implicit none
real(dp),allocatable::r(:)
real(dp)::y(4),mu(4)
integer::fail
fail=0
y=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
mu=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
call qres_invgauss(y,mu,df=3,dispersion=0.5_dp,resid=r)
if(maxval(abs(r))>1e-14_dp)then
print *,'FAIL inverse gaussian qres center'
fail=fail+1
end if
call qres_default([1.0_dp,-2.0_dp],[1.0_dp,1.0_dp],df=2,dispersion=4.0_dp,resid=r)
if(maxval(abs(r-[0.5_dp,-1.0_dp]))>1e-14_dp)then
print *,'FAIL default qres'
fail=fail+1
end if
call set_seed_int(321)
call qres_pois([0.0_dp,1.0_dp,3.0_dp],[0.5_dp,1.0_dp,2.5_dp],r)
if(any(.not.ieee_is_finite(r)))then
print *,'FAIL poisson qres'
fail=fail+1
end if
call qres_gamma([0.5_dp,1.0_dp,3.0_dp],[1.0_dp,1.0_dp,2.0_dp],df=2,dispersion=0.5_dp,resid=r)
if(any(.not.ieee_is_finite(r)))then
print *,'FAIL gamma qres'
fail=fail+1
end if
call set_seed_int(654)
call qres_nbinom([0.0_dp,2.0_dp,5.0_dp],[1.0_dp,2.0_dp,4.0_dp],2.5_dp,r)
if(any(.not.ieee_is_finite(r)))then
print *,'FAIL NB qres'
fail=fail+1
end if
if(fail>0)error stop 'test_qres failed'
print '(a)','test_qres: PASS'
end program
