! Computational routines derived from statmod R/digamma.R and R/digammaf.R.
! Upstream license: GPL-2 | GPL-3. See LICENSE and NOTICE.md.
module statmod_special
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use r_compat, only: dp, r_digamma, r_trigamma
implicit none
private
public :: logmdigamma, cumulant_digamma, meanval_digamma, d2cumulant_digamma
public :: canonic_digamma, varfun_digamma, unitdeviance_digamma
contains

pure elemental function logmdigamma(x) result(ps)
real(dp), intent(in) :: x
real(dp) :: ps,z,t,tail
integer :: k
if(x<=0.0_dp) then
   ps=ieee_value(ps,ieee_quiet_nan)
   return
end if
z=x
if(z<5.0_dp) then
   z=x+5.0_dp
   t=1.0_dp/z**2
   tail=t*(-1.0_dp/12.0_dp+t*(1.0_dp/120.0_dp+t*(-1.0_dp/252.0_dp+t*(1.0_dp/240.0_dp+ &
      t*(-1.0_dp/132.0_dp+t*(691.0_dp/32760.0_dp+t*(-1.0_dp/12.0_dp+3617.0_dp*t/8160.0_dp)))))))
   ps=1.0_dp/(2.0_dp*z)-tail+log(x/z)
   do k=0,4
      ps=ps+1.0_dp/(x+real(k,dp))
   end do
else
   t=1.0_dp/z**2
   tail=t*(-1.0_dp/12.0_dp+t*(1.0_dp/120.0_dp+t*(-1.0_dp/252.0_dp+t*(1.0_dp/240.0_dp+ &
      t*(-1.0_dp/132.0_dp+t*(691.0_dp/32760.0_dp+t*(-1.0_dp/12.0_dp+3617.0_dp*t/8160.0_dp)))))))
   ps=1.0_dp/(2.0_dp*z)-tail
end if
end function logmdigamma

pure elemental function cumulant_digamma(theta) result(v)
real(dp),intent(in)::theta
real(dp)::v
v=2.0_dp*(theta*(log(-theta)-1.0_dp)+log_gamma(-theta))
end function cumulant_digamma

pure elemental function meanval_digamma(theta) result(v)
real(dp),intent(in)::theta
real(dp)::v
v=2.0_dp*logmdigamma(-theta)
end function meanval_digamma

pure elemental function d2cumulant_digamma(theta) result(v)
real(dp),intent(in)::theta
real(dp)::v
v=2.0_dp*(1.0_dp/theta+r_trigamma(-theta))
if(theta < -1.0e3_dp) v=(1.0_dp-1.0_dp/(3.0_dp*theta))/theta**2
end function d2cumulant_digamma

pure elemental function canonic_digamma(mu) result(theta)
real(dp),intent(in)::mu
real(dp)::theta,mlmt,mu1,varc,deriv
integer::i
mlmt=log(mu)
theta=-exp(-mlmt)
do i=1,3
   mu1=meanval_digamma(theta)
   varc=d2cumulant_digamma(theta)
   deriv=-varc/mu1*theta
   mlmt=mlmt-log(mu1/mu)/deriv
   theta=-exp(-mlmt)
end do
end function canonic_digamma

pure elemental function varfun_digamma(mu) result(v)
real(dp),intent(in)::mu
real(dp)::v,theta
theta=canonic_digamma(mu)
v=2.0_dp*(1.0_dp/theta+r_trigamma(-theta))
end function varfun_digamma

pure elemental function unitdeviance_digamma(y,mu) result(v)
real(dp),intent(in)::y,mu
real(dp)::v,thetay,theta
thetay=canonic_digamma(y)
theta=canonic_digamma(mu)
v=2.0_dp*(y*(thetay-theta)-(cumulant_digamma(thetay)-cumulant_digamma(theta)))
end function unitdeviance_digamma

end module statmod_special
