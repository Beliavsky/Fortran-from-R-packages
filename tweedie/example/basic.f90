program basic
use tweedie
implicit none
real(dp), parameter :: mu=1.0_dp, phi=1.0_dp, power=1.5_dp
real(dp) :: y, d, p
y=1.0_dp
d=dtweedie(y,mu,phi,power)
p=ptweedie(y,mu,phi,power)
print '(a,f8.3)','y       = ',y
print '(a,es14.6)','density = ',d
print '(a,f12.8)','cdf     = ',p
end program basic
