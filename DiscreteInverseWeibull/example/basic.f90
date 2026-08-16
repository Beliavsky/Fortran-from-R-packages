program basic
 use discrete_inverse_weibull
 implicit none
 real(dp)::x(100)
 type(diw_estimate)::e
 call set_rng_seed(123);call rdiweibull(x,0.5_dp,2.5_dp)
 e=estdiweibull(x,'H')
 print '(a,f10.6)','q    = ',e%q
 print '(a,f10.6)','beta = ',e%beta
 print '(a,f12.6)','negative log-likelihood = ',e%objective
end program
