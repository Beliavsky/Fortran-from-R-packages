program basic
use truncnorm, only: dp, dtruncnorm, ptruncnorm, qtruncnorm, etruncnorm, vtruncnorm, rtruncnorm, set_seed_int
implicit none
real(dp) :: a,b,mu,sig,x,p
a=-1.0_dp
b=2.0_dp
mu=0.3_dp
sig=1.4_dp
x=0.2_dp
p=ptruncnorm(x,a,b,mu,sig)
print '(a,f12.8)', 'density  = ', dtruncnorm(x,a,b,mu,sig)
print '(a,f12.8)', 'cdf      = ', p
print '(a,f12.8)', 'quantile = ', qtruncnorm(p,a,b,mu,sig)
print '(a,f12.8)', 'mean     = ', etruncnorm(a,b,mu,sig)
print '(a,f12.8)', 'variance = ', vtruncnorm(a,b,mu,sig)
call set_seed_int(12345)
print '(a,f12.8)', 'draw     = ', rtruncnorm(a,b,mu,sig)
end program basic
